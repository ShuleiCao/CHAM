import os
import h5py
from joblib import Parallel, delayed
import logging
import tempfile
import argparse

parser = argparse.ArgumentParser(description='Matching arguments.')
parser.add_argument('--lambda_cut_suffix', type=str, default='_lgt20', 
                    help='Suffix for cluster richness cut, default: _lgt20 (for richness cut 20)')
parser.add_argument('--output_loc', type=str, default='/path/to/your_output/', 
                    help='Specify folder location for outputs')
parser.add_argument('--temp_dir', type=str, default='/path/to/your_temp_folder/', 
                    help='Specify temporary folder for nuisance outputs')

args = parser.parse_args()
tempfile.tempdir = args.temp_dir

# Configure logging to capture errors and warnings
logging.basicConfig(
    filename='consolidate_errors.log',
    level=logging.INFO,  # Set to INFO to capture detailed logs
    format='%(asctime)s:%(levelname)s:%(message)s'
)

def process_file(filename, input_folder):
    cluster_id = filename.split('.')[0]
    try:
        with h5py.File(os.path.join(input_folder, filename), 'r') as cluster_file:

            # Process other categories
            categories_data = {}
            for category in cluster_file.keys():
                categories_data[category] = {}
                for halo_id in cluster_file[category].keys():
                    halo_group = cluster_file[category][halo_id]
                    halo_data = {}
                    for dataset in halo_group.keys():
                        data = halo_group[dataset][:]
                        halo_data[dataset] = data  # Preserve original dtype
                    # Extract attributes
                    halo_attrs = {}
                    for attr in halo_group.attrs:
                        halo_attrs[attr] = halo_group.attrs[attr]  # Preserve original dtype
                    halo_data['attrs'] = halo_attrs
                    categories_data[category][halo_id] = halo_data
        return (cluster_id, categories_data)
    except Exception as e:
        logging.error(f"Error processing file {filename}: {e}")
        return None  # Indicate failure for this file

def write_to_hdf5(main_file, cluster_id, categories_data):
    try:
        if cluster_id in main_file:
            logging.warning(f"Cluster {cluster_id} already exists in the consolidated file. Skipping.")
            return

        cluster_group = main_file.create_group(cluster_id)
        logging.info(f"Created group for cluster '{cluster_id}'.")

        for category, category_data in categories_data.items():
            category_group = cluster_group.create_group(category)
            logging.info(f"  Created category group '{category}' for cluster '{cluster_id}'.")
            for halo_id, halo_data in category_data.items():
                halo_group = category_group.create_group(halo_id)
                logging.info(f"    Created halo group '{halo_id}' in category '{category}'.")
                for dataset, data in halo_data.items():
                    if dataset == 'attrs':
                        continue  # Attributes are handled separately
                    halo_group.create_dataset(dataset, data=data)
                    logging.info(f"      Written dataset '{dataset}' with data shape {data.shape}.")
                # Write attributes
                for attr, value in halo_data['attrs'].items():
                    halo_group.attrs[attr] = value
                    logging.info(f"      Written attribute '{attr}' with value '{value}'.")
    except Exception as e:
        logging.error(f"Error writing cluster {cluster_id} to HDF5: {e}")

def consolidate_hdf5_files_with_centrals_parallel(input_folder, output_file, num_workers=4, batch_size=1000, dry_run=False, test_size=100):
    # Get the list of files to process
    files = [filename for filename in os.listdir(input_folder) if filename.endswith('.h5')]

    # Check if there are any files to process
    if not files:
        print("No files found in input folder.")
        return

    if dry_run:
        files = files[:test_size]
        print(f"Dry run: Processing {len(files)} out of {len(files)} files.")
    else:
        print(f"Processing {len(files)} files using {num_workers} workers in batches of {batch_size} ...")

    # Initialize the output HDF5 file
    if not dry_run:
        with h5py.File(output_file, 'w') as main_file:
            pass  # Create or clear the file

    # Process in batches to manage memory usage
    for i in range(0, len(files), batch_size):
        batch_files = files[i:i + batch_size]
        print(f"Processing batch {i//batch_size + 1} with {len(batch_files)} files...")
        
        # Parallel processing of the batch
        results = Parallel(n_jobs=num_workers, verbose=0)(
            delayed(process_file)(filename, input_folder) for filename in batch_files
        )

        # Open the main HDF5 file in append mode
        with h5py.File(output_file, 'a') as main_file:
            for res in results:
                if res is None:
                    continue  # Skip files that failed to process
                cluster_id, categories_data = res
                write_to_hdf5(main_file, cluster_id, categories_data)

        if dry_run:
            print("Dry run complete! No data has been written.")
            break

    if not dry_run:
        print("Consolidation complete!")

if __name__ == "__main__":
    # Initialize SLURM settings and paths as before
    assigned_cpus = os.getenv('SLURM_CPUS_PER_TASK')
    slurm_partition = os.getenv('SLURM_JOB_PARTITION')

    if assigned_cpus is not None:
        assigned_cpus = int(assigned_cpus)
        print(f'Assigned CPUs: {assigned_cpus}')
    else:
        print('Not running under SLURM or the variable is not set.')
        assigned_cpus = 4  # Default to 4 if not set

    job_num = assigned_cpus
    print(f'Job number is {job_num}')

    suffix = args.lambda_cut_suffix
    
    output_file = os.path.join(args.output_loc, f'sorted_clusters{suffix}.h5')
    
    input_folder = os.path.join(args.output_loc, f'clusters{suffix}')
    
    # Perform a dry run first
    # consolidate_hdf5_files_with_centrals_parallel(
    #     input_folder=input_folder, 
    #     output_file=output_file, 
    #     num_workers=job_num, 
    #     batch_size=1000, 
    #     dry_run=True,  # Set to True for a dry run
    #     test_size=10  # Number of files to test
    # )
    
    # If dry run is successful, run the full consolidation
    consolidate_hdf5_files_with_centrals_parallel(
        input_folder=input_folder, 
        output_file=output_file, 
        num_workers=job_num, 
        batch_size=1000, 
        dry_run=False, 
        # test_size=100
    )
