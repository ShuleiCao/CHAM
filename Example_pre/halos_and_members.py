import h5py
from astropy.io import fits
import numpy as np
import healpy as hp
import os
from astropy.table import Table, vstack
from joblib import Parallel, delayed
import gc
from tqdm import tqdm
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

import tempfile
import argparse

parser = argparse.ArgumentParser(description='Matching arguments.')
parser.add_argument('--mock_path', type=str, default='/path/to/your_mock/Cardinal-3_v2.0_Y6a_gold.h5',
                    help='Specify mock path (without redshifts)')
parser.add_argument('--redshift_path', type=str, default='/path/to/your_mock/Cardinal-3_v2.0_Y6a_bpz.h5',
                    help='Specify mock path for redshifts')
parser.add_argument('--member_path', type=str, default='/path/to/your_cluster_member/original_member_data_lgt20.fit',
                    help='Specify cluster member path')
parser.add_argument('--lambda_cut_suffix', type=str, default='_lgt20', 
                    help='Suffix for cluster richness cut, default: _lgt20 (for richness cut 20)')
parser.add_argument('--output_loc', type=str, default='/path/to/your_output/', 
                    help='Specify folder location for outputs')
parser.add_argument('--keys', type=str, 
                    default='px,py,pz,haloid,m200,r200,mem_match_id,coadd_object_id,id,ra,dec,z', 
                    help=("Specify the keys as a comma-separated string (e.g., 'px,py,pz'). "
                          "Default: 'px,py,pz,haloid,m200,r200,rhalo,coadd_object_id,id,mem_match_id,ra,dec,z' "
                          "(Includes true positions, halo ID, halo mass, halo radius, "
                          "member galaxy ID, member galaxy ID, cluster ID, RA, DEC, and redshift.) "
                          "Please pass the keys in the same order as the default, otherwise modify as you see fit."))
parser.add_argument('--temp_dir', type=str, default='/path/to/your_temp_folder/', 
                    help='Specify temporary folder for nuisance outputs')

args = parser.parse_args()
tempfile.tempdir = args.temp_dir


def correct_endianness(data):
    for col in data.keys():
        if data[col].dtype.byteorder == '>':
            data[col] = data[col].byteswap().newbyteorder()

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

# Function to compute pixel ID from RA, Dec
def compute_pixel_id(ra, dec, nside=8):
    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)
    return hp.ang2pix(nside, theta, phi)

# Load specific columns from a dataset by pixel ID
def load_data_for_pixel(file_path, pixel_id, columns, pixel_ids_all, file_format='fits', key=None):
    relevant_mask = pixel_ids_all == pixel_id
    try:
        if file_format == 'hdf5':
            with h5py.File(file_path, 'r') as f:
                return {col: f['catalog/'][key][col][:][relevant_mask] for col in columns}
        elif file_format == 'fits':
            with fits.open(file_path) as hdul:
                data_table = hdul[1].data
                data = {col: data_table[col][relevant_mask] for col in columns}
                return data
    except TypeError:
        # Debugging information
        print(f"File Path: {file_path}")
        print(f"Pixel ID: {pixel_id}")
        print(f"Columns: {columns}")
        print(f"File Format: {file_format}")
        print(f"Key: {key}")
        raise

def get_relevant_neighboring_pixels(pixel_id, mock_pixel_ids, nside=8):
    """Return a list of neighboring pixel IDs for a given pixel from the mock dataset."""
    theta, phi = hp.pix2ang(nside, pixel_id)
    neighbors = hp.get_all_neighbours(nside, theta, phi)
    neighbors = neighbors[~np.isnan(neighbors)].astype(int)
    return [pid for pid in neighbors if pid in mock_pixel_ids]

def load_data_for_pixel_and_neighbors(mock_data_path, pixel_id, columns, mock_pixel_ids, file_format='hdf5', key='gold'):
    # Start with the main pixel's data
    data = load_data_for_pixel(mock_data_path, pixel_id, columns, mock_pixel_ids, file_format, key)
    
    # Get the relevant neighboring pixels for the given pixel from the mock dataset
    pixels_to_load = get_relevant_neighboring_pixels(pixel_id, mock_pixel_ids)
    
    # Ensure consistent order for processing neighboring pixels
    pixels_to_load.sort()
    
    for pid in pixels_to_load:
        pixel_data = load_data_for_pixel(mock_data_path, pid, columns, mock_pixel_ids, file_format, key)
        for col in columns:
            data[col] = np.concatenate([data[col], pixel_data[col]])
    return data

nside = 8

suffix = args.lambda_cut_suffix
mock_path = args.mock_path
member_path = args.member_path
redshift_path = args.redshift_path
keys = args.keys
output_loc = args.output_loc

mock_pixels_path = os.path.dirname(mock_path, f'mock_pixels_nside{nside}.npz')
if os.path.isfile(mock_pixels_path):
    mock_pixels = np.load(mock_pixels_path)['mock_pixels']
else:
    mock_radec = get_name(mock_path, keys[-3:-1], file_format='hdf5', key='gold')
    mock_pixels = compute_pixel_id(mock_radec[keys[-3]], mock_radec[keys[-2]])
    del mock_radec
    np.savez(mock_pixels_path, mock_pixels=mock_pixels)

member_pixels_path = os.path.join(output_loc, f'member{suffix}_pixels_nside{nside}.npz')
if os.path.isfile(member_pixels_path):
    member_pixels = np.load(os.path.join(output_loc, f'member{suffix}_pixels_nside{nside}.npz'))['member_pixels']
else:
    member_radec = get_name(member_path, keys[-3:-1], file_format='fits')
    member_pixels = compute_pixel_id(member_radec[keys[-3]], member_radec[keys[-2]])
    del member_radec
    np.savez(member_pixels_path, member_pixels=member_pixels)

unique_mock_pixels = np.unique(mock_pixels)
unique_member_pixels = np.unique(member_pixels)

# Save halo properties
def process_halo_pixel(pixel_id):
    print(f"Processing pixel {pixel_id}")

    mock_data_pixel = load_data_for_pixel(mock_path, pixel_id, keys[:8] + keys[-3:-1], mock_pixels, file_format='hdf5', key='gold')

    # Load cosmological redshift (redshift_cos) from bpz (need to be modified for your own usage!!!)
    z_values = load_data_for_pixel_and_neighbors(redshift_path, pixel_id, ['redshift_cos'], mock_pixels, file_format='hdf5', key='bpz')
    mock_data_pixel[keys[-1]] = z_values['redshift_cos']
    del z_values

    # Filter halo data where m200 > 0 and rhalo = 0
    halo_mask = (mock_data_pixel['rhalo'] == 0) & (mock_data_pixel['m200'] > 0)
    halo_data_pixel_table = Table({col: list(mock_data_pixel[col][halo_mask]) for col in mock_data_pixel.keys()})

    del halo_mask, mock_data_pixel
    gc.collect()
    
    return halo_data_pixel_table

results = Parallel(n_jobs=-1, verbose=1)(delayed(process_halo_pixel)(pixel_id) 
                                         for pixel_id in tqdm(unique_mock_pixels, 
                                                              desc="Processing pixels"))

# Convert to a single Table
halo_data = vstack(results)
halo_path = os.path.join(output_loc, 'halo_data.fits')
halo_data.write(halo_path, overwrite=True)
del halo_data, results
gc.collect()

# Save cluster member properties
def process_cluster_pixel(pixel_id):
    print(f"Processing pixel {pixel_id}")

    mock_data_pixel = load_data_for_pixel_and_neighbors(mock_path, pixel_id, keys[:8], mock_pixels, file_format='hdf5', key='gold') # suppose 'coadd_object_id' is used
    member_data_pixel = load_data_for_pixel(member_path, pixel_id, keys[-5:-1], member_pixels, file_format='fits') # suppose 'id' is used

    # Load cosmological redshift (redshift_cos) from bpz (need to be modified for your own usage!!!)
    z_values = load_data_for_pixel_and_neighbors(redshift_path, pixel_id, ['redshift_cos'], mock_pixels, file_format='hdf5', key='bpz')
    mock_data_pixel[keys[-1]] = z_values['redshift_cos']
    del z_values
    
    correct_endianness(mock_data_pixel)
    correct_endianness(member_data_pixel)

    # Convert dictionaries to pandas DataFrames
    mock_data_pixel = pd.DataFrame(mock_data_pixel)
    member_data_pixel = pd.DataFrame(member_data_pixel)
    
    # Merge DataFrames on 'coadd_object_id' and 'id' (if same then change 'left_on' and 'right_on' into 'on')
    merged_data = pd.merge(mock_data_pixel, member_data_pixel, left_on='coadd_object_id', right_on='id', how='inner')
    merged_data = merged_data.drop(columns='id') # keep only 'coadd_object_id'

    # Convert the merged DataFrame back to an astropy Table
    matched_member_data_pixel_table = Table.from_pandas(merged_data)
    del mock_data_pixel, member_data_pixel, merged_data
    
#     print(f"Length of matched_member_data_pixel_table: {len(matched_member_data_pixel_table)}")
    return matched_member_data_pixel_table

results = Parallel(n_jobs=-1, verbose=1)(delayed(process_cluster_pixel)(pixel_id) 
                                         for pixel_id in tqdm(unique_member_pixels, 
                                                              desc="Processing pixels"))

member_data = vstack(results)
matched_member_path = os.path.join(output_loc, f'member_data{suffix}.fits')
member_data.write(matched_member_path, overwrite=True)
