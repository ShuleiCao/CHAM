import h5py
import numpy as np
import os
from astropy.table import Table, vstack
from scipy.spatial import cKDTree
import gc
import logging
from tqdm import tqdm
from multiprocessing import Manager, Pool, Lock
import argparse

import tempfile
temp_dir = "/lustre/scratch/client/users/shuleic/temp_folder/"
tempfile.tempdir = temp_dir

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

# shared_vars_lock = Lock()

def process_cluster(cluster_id, redmapper_mem_data_path, halo_coords, halo_r200, halo_ids, halo_masses, tree_processes, output_folder):
    cluster_filename = os.path.join(output_folder, f"{cluster_id}.h5")
    # Check if the output file already exists
    if os.path.exists(cluster_filename):
        logging.info(f"Skipping cluster {cluster_id} as the output file already exists.")
        return  # Skip processing
    
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
    galaxy_indices = galaxy_tree.query_ball_point(halo_coords, halo_r200, workers=tree_processes)

    # Process the results
    all_associated_indices = []
    associated_halo_data = []
    
    for halo_id, halo_mass, member_indices in zip(halo_ids, halo_masses, galaxy_indices):
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
    parser = argparse.ArgumentParser(description='Matching arguments.')
    parser.add_argument('--mismatched', action='store_true',
                        help='There exists mismatched or not')
    parser.add_argument('--masked', action='store_true',
                        help='Halos and galaxies masked or not')
    parser.add_argument('--base_path', type=str, default='/lustre/work/client/users/shuleic/Cardinalv3/',
                        help='Specify base path')
    parser.add_argument('--suffix', type=str, default='_lgt05', choices=['_lgt20', '_lgt05'], required=True, 
                        help='Suffix for redmapper richness cut')
    parser.add_argument('--halo_filename', type=str, default='_lgt05', choices=['_lgt20', '_lgt05'], required=True, 
                        help='Halo path filename')
    
    args = parser.parse_args()
    
    def wrapper(cluster_id):
        return process_cluster(cluster_id, redmapper_mem_data_path, halo_coords, halo_r200, halo_ids, halo_masses, tree_processes, output_folder)

    # Set the logging level to INFO
    logging.basicConfig(level=logging.INFO)

    base_path = args.base_path
    suffix = args.suffix
    
    halo_data_path = os.path.join(base_path, args.halo_filename)
    redmapper_mem_data_path = os.path.join(base_path, f'chunhao-redmapper{suffix}_mem_data_all.fits')

    redmapper_mem_data = get_name(redmapper_mem_data_path, ['mem_match_id'])
    unique_cluster_ids = np.unique(redmapper_mem_data['mem_match_id'])
    del redmapper_mem_data
#     gc.collect()

    output_folder = os.path.join(base_path, f"cluster{suffix}_data")
    os.makedirs(output_folder, exist_ok=True)
    sorted_unique_cluster_ids = sorted(unique_cluster_ids)
    del unique_cluster_ids
    tree_processes = 4
    num_processes = 64

    check_clusters_paths = [os.path.exists(os.path.join(output_folder, f'{cluster_id}.h5')) for cluster_id in sorted_unique_cluster_ids]

    all_exist = all(check_clusters_paths)

    if not all_exist or args.mismatched:
        # Read the data from the saved FITS files only if needed
        halo_data = get_name(halo_data_path, ['m200', 'r200', 'haloid', 'px', 'py', 'pz'])
        halo_coords = np.vstack((halo_data['px'], halo_data['py'], halo_data['pz'])).T
        halo_r200 = halo_data['r200']
        halo_ids = halo_data['haloid']
        halo_masses = halo_data['m200']
        del halo_data
        gc.collect()

    # Determine which clusters to process
    clusters_to_process = []

    if all_exist:
        # If all files exist and mismatched flag is set, process mismatched clusters
        if args.mismatched:
            mismatched_ids = np.load(os.path.join(base_path, f'mismatched_clusters{suffix}.npz'))['mismatched_ids']
            clusters_to_process.extend(mismatched_ids)

            # remove existing files for mismatched clusters
            for mismatched_id in mismatched_ids:
                mismatched_cluster_path = os.path.join(output_folder, f'{mismatched_id}.h5')
                if os.path.isfile(mismatched_cluster_path):
                    os.remove(mismatched_cluster_path)
        else:
            logging.info("All clusters exist. Exiting without processing.")
            exit(0)

    else:
        # If not all files exist, process the non-existing clusters
        existing_files = set(os.listdir(output_folder))  # Get all files in the output folder
        clusters_to_process = [cluster_id for cluster_id in sorted_unique_cluster_ids if f"{cluster_id}.h5" not in existing_files]

    # Proceed with processing the determined clusters
    if clusters_to_process:
        with tqdm(total=len(clusters_to_process), desc="Processing clusters") as pbar, Pool(processes=num_processes) as pool:
            for result in pool.imap(wrapper, clusters_to_process):
                pbar.update()  # Update progress for each completed task
    else:
        logging.info("No clusters to process.")

    # Use multiprocessing with Pool and tqdm
    # with tqdm(total=len(sorted_unique_cluster_ids), desc="Processing clusters") as pbar, Pool(processes=num_processes) as pool:
    #     for cluster_id in clusters_to_process:
    #         pool.apply_async(
    #             process_cluster,
    #             args=(cluster_id, redmapper_mem_data_path, tree_processes, shared_variables, output_folder),
    #             callback=lambda _: pbar.update(1)  # Update progress for each completed task
    #         )
    #     pool.close()
    #     pool.join()

