import h5py
from astropy.io import fits
import numpy as np
import healpy as hp
import os
from astropy.table import Table, vstack
from joblib import Parallel, delayed
from collections import defaultdict
from scipy.spatial import cKDTree
from joblib import Parallel, delayed
import gc
import logging
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from memory_profiler import profile
from multiprocessing import Process, shared_memory
import multiprocessing
from multiprocessing.pool import ThreadPool as Pool
#from cluster_processing import process_cluster, get_name
from functools import partial
from multiprocessing import Lock
import tempfile
temp_dir = "/lustre/scratch/client/users/shuleic/temp_folder/"
tempfile.tempdir = temp_dir

## Initialize logging
#logging.basicConfig(level=logging.INFO)

def get_name(file_path, names, mask=None, file_format='fits', key='gold'):
    if file_format == 'fits':
        table_data = Table.read(file_path)  # Directly read the FITS file into an Astropy Table
        if mask is not None:
            table_data = table_data[mask]  # Apply the mask directly to the table

        data_names = {name: table_data[name] for name in names}  # Extract the relevant columns
    elif file_format == 'hdf5':
        with h5py.File(file_path, 'r') as f:
            data_names = {name: f['catalog/'][key][name][:][mask] if mask is not None else f['catalog/'][key][name][:] for name in names}
    else:
        raise ValueError(f"Unknown file format: {file_format}")
    return data_names

def init_logger(log_queue):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # Create a handler that logs to the Queue
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    handler.setStream(log_queue)
    logger.addHandler(handler)

def save_cluster_data_hdf5(cluster_id, data, output_folder):
    cluster_filename = os.path.join(output_folder, f"{cluster_id}.h5")
    temp_filename = os.path.join(output_folder, f"{cluster_id}_temp.h5")

    with h5py.File(temp_filename, 'w') as f:
        for category, data_list in data.items():
            category_group = f.create_group(category)
            for halo_data in data_list:
                halo_id = str(halo_data['haloid'])  # Use halo ID as the group name
                halo_group = category_group.create_group(halo_id)
                for col_name, col_data in halo_data['members'].items():
                    halo_group.create_dataset(col_name, data=col_data)
                halo_group.attrs['n_members'] = halo_data['n_members']
                halo_group.attrs['m200'] = halo_data['m200']
    os.rename(temp_filename, cluster_filename)

shared_vars_lock = Lock()

def process_cluster(cluster_id, redmapper_mem_data_path, tree_processes, shared_variables, output_folder):
    cluster_filename = os.path.join(output_folder, f"{cluster_id}.h5")
    # Check if the output file already exists
    if os.path.exists(cluster_filename):
        logging.info(f"Skipping cluster {cluster_id} as the output file already exists.")
        return  # Skip processing
    
    with shared_vars_lock:
        halo_coords_shared = shared_variables['halo_coords_shared']
        halo_r200_shared = shared_variables['halo_r200_shared']
        halo_ids_shared = shared_variables['halo_ids_shared']
        halo_masses_shared = shared_variables['halo_masses_shared']

    # 1. Filter members of the current cluster
    redmapper_mem_data = get_name(redmapper_mem_data_path, ['mem_match_id'])
    cluster_member_mask = redmapper_mem_data['mem_match_id'] == cluster_id
    cluster_member_indices = np.where(cluster_member_mask)[0]

#    logging.info("Creating cKDTree from galaxy positions...")
    redmapper_mem_data = get_name(redmapper_mem_data_path, ['px','py','pz'], cluster_member_indices)
    redmapper_mem_coords = np.vstack((redmapper_mem_data['px'], redmapper_mem_data['py'], redmapper_mem_data['pz'])).T
    del redmapper_mem_data
    gc.collect()

    galaxy_tree = cKDTree(redmapper_mem_coords)
    galaxy_indices = galaxy_tree.query_ball_point(halo_coords_shared, halo_r200_shared, workers=tree_processes)

    # # Create a mapping between original galaxy indices and their indices in redmapper_mem_coords
    # galaxy_mapping = {}
    
    # for original_idx, redmapper_mem_idx in enumerate(cluster_member_indices):
    #     galaxy_mapping[redmapper_mem_idx] = original_idx

    # Process the results
    all_associated_indices = []
    associated_halo_data = []
    
    for halo_id, halo_mass, member_indices in zip(halo_ids_shared, halo_masses_shared, galaxy_indices):
        if len(member_indices) > 0:
            member_data_combined = {}        
            for member_idx in member_indices:
                if member_idx < len(cluster_member_indices):  # Ensure member_idx is within bounds
                    original_idx = cluster_member_indices[member_idx]
                    all_associated_indices.append(original_idx)  # Save the original index
                    member_data = get_name(redmapper_mem_data_path, ['coadd_object_id', 'haloid', 'ra', 'dec', 'rhalo', 'r200', 'm200', 'px', 'py', 'pz', 'z'], mask=[original_idx])
    
                    # Aggregate data for each attribute, also ensuring single values are not wrapped in arrays
                    for key, value in member_data.items():
                        if key not in member_data_combined:
                            member_data_combined[key] = []
                        member_data_combined[key].append(value[0] if isinstance(value, (list, tuple, np.ndarray)) and len(value) == 1 else value)
                    
            associated_halo_data.append({
                'haloid': halo_id,
                'm200': halo_mass,
                'n_members': len(member_indices),
                'members': member_data_combined
            })
    
    for halo_data in associated_halo_data:
        for key, value in halo_data['members'].items():
            halo_data['members'][key] = np.array(value)
    
#     associated_halo_data_sorted = sorted(associated_halo_data, key=lambda x: x['n_members'], reverse=True)

#    logging.info("Finished processing assoicated data...")
        
    # Determine unassociated members
    unassociated_indices = [idx for idx in cluster_member_indices if idx not in all_associated_indices]
    unassociated_data = get_name(redmapper_mem_data_path, ['coadd_object_id', 'haloid', 'ra', 'dec', 'rhalo', 'r200', 'm200', 'px', 'py', 'pz','z'], mask=unassociated_indices)
    #    logging.info("Finished loading unassoicated data...")

    combined_key_str = np.array([f"{haloid}_{m200}" for haloid, m200 in zip(unassociated_data['haloid'], unassociated_data['m200'])], dtype='U25')
    unique_combined_key, counts = np.unique(combined_key_str, return_counts=True)

    unassociated_halo_data = []
    
    for combined_key, count in zip(unique_combined_key, counts):
        halo_id, halo_mass = map(float, combined_key.split('_'))
        halo_id = np.int64(halo_id)
        halo_mass = np.float32(halo_mass)
        halo_mask = (unassociated_data['haloid'] == halo_id) & (unassociated_data['m200'] == halo_mass)

        member_data = {
            key: value[halo_mask]
            for key, value in unassociated_data.items()
        }
        unassociated_halo_data.append({
            'haloid': halo_id,
            'm200': halo_mass,
            'n_members': count,
            'members': member_data
        })
    
#     unassociated_halo_data_sorted = sorted(unassociated_halo_data, key=lambda x: x['n_members'], reverse=True)
#     logging.info("Finished processing unassociated data...")
    
    # Sort the data for unassociated galaxies
#     unassociated_halo_data.sort(key=lambda x: x['n_members'], reverse=True)

#    logging.info("Check...")
    
    save_cluster_data_hdf5(cluster_id, {
        'associated': associated_halo_data,
        'unassociated': unassociated_halo_data
    }, output_folder)

if __name__ == "__main__":
    # Create a multiprocessing Queue to pass log messages
    log_queue = multiprocessing.Queue()

    # Create a logger process and initialize the logger
    logger_process = multiprocessing.Process(target=init_logger, args=(log_queue,))
    logger_process.start()

    # Set the logging level to INFO
    logging.basicConfig(level=logging.INFO)

    base_path = '/lustre/work/client/users/shuleic/Cardinalv3/'
    # suffix = '_lgt20'
    suffix = '_lgt05'
    
    halo_data_path = os.path.join(base_path, 'halo_data_all_new.fits')
    redmapper_mem_data_path = os.path.join(base_path, f'chunhao-redmapper{suffix}_mem_data_all_new.fits')
#     redmapper_centrals_data_path = os.path.join(base_path, 'chunhao-redmapper_centrals_data_all_new.fits')

    # Read the data from the saved FITS files
    halo_data = get_name(halo_data_path, ['m200', 'r200', 'haloid', 'px', 'py', 'pz'])
    halo_coords = np.vstack((halo_data['px'], halo_data['py'], halo_data['pz'])).T
    halo_r200 = halo_data['r200']
    halo_ids = halo_data['haloid']
    halo_masses = halo_data['m200']
    del halo_data
    gc.collect()

    # Create shared memory for halo_coords, halo_r200, and halo_ids
    halo_coords_shm = shared_memory.SharedMemory(create=True, size=halo_coords.nbytes)
    halo_r200_shm = shared_memory.SharedMemory(create=True, size=halo_r200.nbytes)
    halo_ids_shm = shared_memory.SharedMemory(create=True, size=halo_ids.nbytes)
    halo_masses_shm = shared_memory.SharedMemory(create=True, size=halo_masses.nbytes)

    # Copy the data to shared memory
    halo_coords_shared = np.ndarray(halo_coords.shape, dtype=halo_coords.dtype, buffer=halo_coords_shm.buf)
    halo_coords_shared[:] = halo_coords[:]
    halo_r200_shared = np.ndarray(halo_r200.shape, dtype=halo_r200.dtype, buffer=halo_r200_shm.buf)
    halo_r200_shared[:] = halo_r200[:]
    halo_ids_shared = np.ndarray(halo_ids.shape, dtype=halo_ids.dtype, buffer=halo_ids_shm.buf)
    halo_ids_shared[:] = halo_ids[:]
    halo_masses_shared = np.ndarray(halo_masses.shape, dtype=halo_masses.dtype, buffer=halo_masses_shm.buf)
    halo_masses_shared[:] = halo_masses[:]

    del halo_coords, halo_r200, halo_ids, halo_masses
#     gc.collect()
    
    logging.info("Processing clusters...")
    tree_processes = 4

    # Get the number of CPUs assigned to the current task
    assigned_cpus = os.getenv('SLURM_CPUS_PER_TASK')

    if assigned_cpus is not None:
        assigned_cpus = int(assigned_cpus)
        print(f'Assigned CPUs: {assigned_cpus}')
    else:
        print('Not running under SLURM or the variable is not set.')

    num_threads = assigned_cpus//2 # Adjust this if necessary

    redmapper_mem_data = get_name(redmapper_mem_data_path, ['mem_match_id'])
    unique_cluster_ids = np.unique(redmapper_mem_data['mem_match_id'])
    del redmapper_mem_data
#     gc.collect()

    shared_variables = {
        'halo_coords_shared': halo_coords_shared,
        'halo_r200_shared': halo_r200_shared,
        'halo_ids_shared': halo_ids_shared,
        'halo_masses_shared': halo_masses_shared
    }

    output_folder = os.path.join(base_path, f"cluster{suffix}_data")
    os.makedirs(output_folder, exist_ok=True)
    sorted_unique_cluster_ids = sorted(unique_cluster_ids)
    del unique_cluster_ids
    
    # Use the sorted list in the Pool
    with tqdm(total=len(sorted_unique_cluster_ids), desc="Processing clusters") as pbar, Pool(processes=num_threads) as pool:
        def update_progress(_):
            pbar.update(1)

        for cluster_id in sorted_unique_cluster_ids:
            pool.apply_async(
                process_cluster,
                args=(cluster_id, redmapper_mem_data_path, tree_processes, shared_variables, output_folder),
                callback=update_progress  # Corrected callback
            )

        # Wait for all tasks to complete
        pool.close()
        pool.join()

    logger_process.join()

    # Clean up shared memory when it's no longer needed
    halo_coords_shm.close()
    halo_coords_shm.unlink()
    halo_r200_shm.close()
    halo_r200_shm.unlink()
    halo_ids_shm.close()
    halo_ids_shm.unlink()
    halo_masses_shm.close()
    halo_masses_shm.unlink()
