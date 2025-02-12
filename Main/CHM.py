import h5py
import numpy as np
import os
from astropy.table import Table
from scipy.spatial import cKDTree
import gc
import logging
from tqdm import tqdm
from multiprocessing import Manager, Pool
import tempfile
import argparse

parser = argparse.ArgumentParser(description='Matching arguments.')
parser.add_argument('--member_path', type=str, default='/path/to/your_cluster_member/member_data_lgt20.fits',
                    help='Specify cluster member path')
parser.add_argument('--lambda_cut_suffix', type=str, default='_lgt20', 
                    help='Suffix for cluster richness cut, default: _lgt20 (for richness cut 20)')
parser.add_argument('--halo_path', type=str, default='/path/to/your_halo/halo_data.fits', 
                    help='Specify halo path')
parser.add_argument('--output_loc', type=str, default='/path/to/your_output/', 
                    help='Specify folder location for outputs')
parser.add_argument('--keys', type=str, 
                    default='px,py,pz,haloid,m200,r200,mem_match_id,coadd_object_id,ra,dec,z', 
                    help=("Specify the keys as a comma-separated string (e.g., 'px,py,pz'). "
                          "Default: 'px,py,pz,haloid,m200,r200,mem_match_id,coadd_object_id,ra,dec,z' "
                          "(Includes true positions, halo ID, halo mass, halo radius, "
                          "cluster ID, member galaxy ID, RA, DEC, and redshift.) "
                          "Please pass the keys in the same order as the default, otherwise modify as you see fit."))
parser.add_argument('--temp_dir', type=str, default='/path/to/your_temp_folder/', 
                    help='Specify temporary folder for nuisance outputs')

args = parser.parse_args()
tempfile.tempdir = args.temp_dir

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

def save_cluster_data_hdf5(cluster_id, data, output_loc):
    cluster_filename = os.path.join(output_loc, f"{cluster_id}.h5")
    temp_filename = os.path.join(output_loc, f"{cluster_id}_temp.h5")

    with h5py.File(temp_filename, 'w') as f:
        for category, data_list in data.items():
            category_group = f.create_group(category)
            for halo_data in data_list:
                halo_id = str(halo_data[haloID_key])  # Use halo ID as the group name
                halo_group = category_group.create_group(halo_id)
                for col_name, col_data in halo_data['members'].items():
                    halo_group.create_dataset(col_name, data=col_data)
                halo_group.attrs['n_members'] = halo_data['n_members']
                halo_group.attrs[haloM_key] = halo_data[haloM_key]
    os.rename(temp_filename, cluster_filename)

def process_cluster(cluster_id, member_path, tree_processes, shared_variables, output_loc):
    cluster_filename = os.path.join(output_loc, f"{cluster_id}.h5")
    # Check if the output file already exists
    if os.path.exists(cluster_filename):
        logging.info(f"Skipping cluster {cluster_id} as the output file already exists.")
        return  # Skip processing
    
    halo_coords_shared = shared_variables['halo_coords']
    halo_radii_shared = shared_variables['halo_radii']
    halo_ids_shared = shared_variables['halo_ids']
    halo_masses_shared = shared_variables['halo_masses']

    # 1. Filter members of the current cluster
    member_clusterIDs = get_name(member_path, [clusterID_key])
    cluster_member_mask = member_clusterIDs[clusterID_key] == cluster_id
    cluster_member_indices = np.where(cluster_member_mask)[0]

#    logging.info("Creating cKDTree from galaxy positions...")
    cluster_member_coords = get_name(member_path, coords_keys, cluster_member_indices)
    cluster_member_coords = np.vstack((cluster_member_coords[coords_keys[0]], 
                                       cluster_member_coords[coords_keys[1]], 
                                       cluster_member_coords[coords_keys[2]])).T

    galaxy_tree = cKDTree(cluster_member_coords)
    galaxy_indices = galaxy_tree.query_ball_point(halo_coords_shared, halo_radii_shared, workers=tree_processes)

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
                    member_data = get_name(member_path, keys, mask=[original_idx])
    
                    # Aggregate data for each attribute, also ensuring single values are not wrapped in arrays
                    for key, value in member_data.items():
                        if key not in member_data_combined:
                            member_data_combined[key] = []
                        member_data_combined[key].append(value[0] if isinstance(value, (list, tuple, np.ndarray)) and len(value) == 1 else value)
                    
            associated_halo_data.append({
                haloID_key: halo_id,
                haloM_key: halo_mass,
                'n_members': len(member_indices),
                'members': member_data_combined
            })

    for halo_data in associated_halo_data:
        for key, value in halo_data['members'].items():
            halo_data['members'][key] = np.array(value)
    
#    logging.info("Finished processing assoicated data...")
        
    # Determine unassociated members
    unassociated_indices = [idx for idx in cluster_member_indices if idx not in all_associated_indices]
    unassociated_data = get_name(member_path, keys, mask=unassociated_indices)
    #    logging.info("Finished loading unassoicated data...")

    combined_key_str = np.array([f"{haloid}_{m200}" for haloid, m200 in zip(unassociated_data[haloID_key], 
                                                                            unassociated_data[haloM_key])], 
                                                                            dtype='U25')
    unique_combined_key, counts = np.unique(combined_key_str, return_counts=True)

    unassociated_galaxy_data = []
    
    for combined_key, count in zip(unique_combined_key, counts):
        halo_id, halo_mass = map(float, combined_key.split('_'))
        halo_id = np.int64(halo_id)
        halo_mass = np.float32(halo_mass)
        halo_mask = (unassociated_data[haloID_key] == halo_id) & (unassociated_data[haloM_key] == halo_mass)

        member_data = {
            key: value[halo_mask]
            for key, value in unassociated_data.items()
        }
        unassociated_galaxy_data.append({
            haloID_key: halo_id,
            haloM_key: halo_mass,
            'n_members': count,
            'members': member_data
        })
    
#     logging.info("Finished processing unassociated data...")
    
    save_cluster_data_hdf5(cluster_id, {
        'associated': associated_halo_data,
        'unassociated': unassociated_galaxy_data
    }, output_loc)

if __name__ == "__main__":

    def wrapper(cluster_id):
        return process_cluster(cluster_id, member_path, tree_processes, shared_variables, output_loc)

    # Set the logging level to INFO
    logging.basicConfig(level=logging.INFO)

    suffix = args.lambda_cut_suffix
    keys = args.keys.split(',')
    coords_keys = keys[:3] # Positions
    halo_keys = keys[:6] # Positions, ID, mass, and radius
    clusterID_key = keys[6] # Cluster ID key
    haloID_key = halo_keys[3]
    haloM_key = halo_keys[4]

    halo_path = args.halo_path
    member_path = args.member_path

    member_clusterIDs = get_name(member_path, [clusterID_key])
    unique_cluster_ids = np.unique(member_clusterIDs[clusterID_key])
    del member_clusterIDs
#     gc.collect()

    output_loc = args.output_loc
    os.makedirs(output_loc, exist_ok=True)
    sorted_unique_cluster_ids = sorted(unique_cluster_ids)
    del unique_cluster_ids
    tree_processes = 4

    # Get the number of CPUs assigned to the current task
    assigned_cpus = os.getenv('SLURM_CPUS_PER_TASK')

    if assigned_cpus is not None:
        assigned_cpus = int(assigned_cpus)
        logging.info(f'Assigned CPUs: {assigned_cpus}')
    else:
        logging.info('Not running under SLURM or the variable is not set.')

    num_processes = assigned_cpus//2 # Adjust this if necessary

    check_clusters_paths = [os.path.exists(os.path.join(output_loc, f'{cluster_id}.h5')) for cluster_id in sorted_unique_cluster_ids]

    all_exist = all(check_clusters_paths)

    if not all_exist:
        # Read the data from the saved FITS files only if needed
        halo_data = get_name(halo_path, halo_keys)
        halo_coords = np.vstack((halo_data[halo_keys[0]], halo_data[halo_keys[1]], halo_data[halo_keys[2]])).T
        halo_ids = halo_data[haloID_key]
        halo_masses = halo_data[haloM_key]
        halo_radii = halo_data[halo_keys[5]]
        del halo_data
        gc.collect()

        # Create a Manager for shared variables
        manager = Manager()
        shared_variables = manager.dict({
            'halo_coords': halo_coords,
            'halo_radii': halo_radii,
            'halo_ids': halo_ids,
            'halo_masses': halo_masses
        })

        del halo_coords, halo_radii, halo_ids, halo_masses
        gc.collect()  # Optional: Garbage collection can help with memory management

    # Determine which clusters to process
    clusters_to_process = []

    if all_exist:
        logging.info("All clusters exist. Exiting without processing.")
        exit(0)
    else:
        # If not all files exist, process the non-existing clusters
        existing_files = set(os.listdir(output_loc))  # Get all files in the output folder
        clusters_to_process = [cluster_id for cluster_id in sorted_unique_cluster_ids if f"{cluster_id}.h5" not in existing_files]

    # Proceed with processing the determined clusters
    if clusters_to_process:
        with tqdm(total=len(clusters_to_process), desc="Processing clusters") as pbar, Pool(processes=num_processes) as pool:
            for result in pool.imap(wrapper, clusters_to_process):
                pbar.update()  # Update progress for each completed task
    else:
        logging.info("No clusters to process.")

