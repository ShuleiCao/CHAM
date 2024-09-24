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

def get_relevant_neighboring_pixels(pixel_id, gold_pixel_ids, nside=8):
    """Return a list of neighboring pixel IDs for a given pixel from the gold dataset."""
    theta, phi = hp.pix2ang(nside, pixel_id)
    neighbors = hp.get_all_neighbours(nside, theta, phi)
    neighbors = neighbors[~np.isnan(neighbors)].astype(int)
    return [pid for pid in neighbors if pid in gold_pixel_ids]

def load_data_for_pixel_and_neighbors(gold_data_path, pixel_id, columns, gold_pixel_ids, file_format='hdf5', key='gold'):
    # Start with the main pixel's data
    data = load_data_for_pixel(gold_data_path, pixel_id, columns, gold_pixel_ids, file_format, key)
    
    # Get the relevant neighboring pixels for the given pixel from the gold dataset
    pixels_to_load = get_relevant_neighboring_pixels(pixel_id, gold_pixel_ids)
    
    # Ensure consistent order for processing neighboring pixels
    pixels_to_load.sort()
    
    for pid in pixels_to_load:
        pixel_data = load_data_for_pixel(gold_data_path, pid, columns, gold_pixel_ids, file_format, key)
        for col in columns:
            data[col] = np.concatenate([data[col], pixel_data[col]])
    return data


# suffix = '_lgt20'
suffix = '_lgt05'

base_path = '/lustre/work/client/users/shuleic/Cardinalv3/'
gold_path = os.path.join(base_path, 'Cardinal-3_v2.0_Y6a_gold.h5')
redmapper_mem_path = os.path.join(base_path, 'redmapper_v4_v8_v51_y6_v7/run/', f'Cardinal-3Y6a_v2.0_run_run_redmapper_v0.8.1{suffix}_vl02_catalog_members.fit')
redmapper_path = os.path.join(base_path, 'redmapper_v4_v8_v51_y6_v7/run/', f'Cardinal-3Y6a_v2.0_run_run_redmapper_v0.8.1{suffix}_vl02_catalog.fit')
bpz_path = os.path.join(base_path, 'Cardinal-3_v2.0_Y6a_bpz.h5')

gold_pixels_path = os.path.join(base_path, 'gold_pixels_nside8.npz')
if os.path.isfile(gold_pixels_path):
    gold_pixels = np.load(gold_pixels_path)['gold_pixels']
else:
    gold_radec = get_name(gold_path, ['ra', 'dec'], file_format='hdf5', key='gold')
    gold_pixels = compute_pixel_id(gold_radec['ra'], gold_radec['dec'])
    del gold_radec
    np.savez(gold_pixels_path, gold_pixels=gold_pixels)

redmapper_mem_pixels_path = os.path.join(base_path, f'chunhao-redmapper{suffix}_mem_pixels_nside8.npz')
if os.path.isfile(redmapper_mem_pixels_path):
    redmapper_mem_pixels = np.load(os.path.join(base_path, f'chunhao-redmapper{suffix}_mem_pixels_nside8.npz'))['redmapper_mem_pixels']
else:
    redmapper_mem_radec = get_name(redmapper_mem_path, ['ra', 'dec'], file_format='fits')
    redmapper_mem_pixels = compute_pixel_id(redmapper_mem_radec['ra'], redmapper_mem_radec['dec'])
    del redmapper_mem_radec
    np.savez(redmapper_mem_pixels_path, redmapper_mem_pixels=redmapper_mem_pixels)

unique_gold_pixels = np.unique(gold_pixels)
unique_redmapper_mem_pixels = np.unique(redmapper_mem_pixels)

def correct_endianness(data):
    for col in data.keys():
        if data[col].dtype.byteorder == '>':
            data[col] = data[col].byteswap().newbyteorder()

def process_redmapper_pixel(pixel_id):
    print(f"Processing pixel {pixel_id}")

    # Load data into pandas DataFrames
    gold_data_pixel = load_data_for_pixel_and_neighbors(gold_path, pixel_id, ['coadd_object_id', 'haloid', 'rhalo', 'r200', 'm200', 'px', 'py', 'pz'], gold_pixels, file_format='hdf5', key='gold')
    redmapper_mem_data_pixel = load_data_for_pixel(redmapper_mem_path, pixel_id, ['ra', 'dec', 'id', 'mem_match_id'], redmapper_mem_pixels, file_format='fits')
    z_values = load_data_for_pixel_and_neighbors(bpz_path, pixel_id, ['redshift_cos'], gold_pixels, file_format='hdf5', key='bpz')
    gold_data_pixel['z'] = z_values['redshift_cos']
    del z_values
    
    correct_endianness(gold_data_pixel)
    correct_endianness(redmapper_mem_data_pixel)

    # Convert dictionaries to pandas DataFrames
    gold_data_pixel = pd.DataFrame(gold_data_pixel)
    redmapper_mem_data_pixel = pd.DataFrame(redmapper_mem_data_pixel)
    
    # Merge DataFrames on 'coadd_object_id' and 'id'
    merged_data = pd.merge(gold_data_pixel, redmapper_mem_data_pixel, left_on='coadd_object_id', right_on='id', how='inner')
    merged_data = merged_data.drop(columns='id')

    # Convert the merged DataFrame back to an astropy Table
    matched_redmapper_mem_data_pixel_table = Table.from_pandas(merged_data)
    del gold_data_pixel, redmapper_mem_data_pixel, merged_data
    
#     print(f"Length of matched_redmapper_mem_data_pixel_table: {len(matched_redmapper_mem_data_pixel_table)}")
    return matched_redmapper_mem_data_pixel_table

results = Parallel(n_jobs=-1, verbose=10)(delayed(process_redmapper_pixel)(pixel_id) for pixel_id in tqdm(unique_redmapper_mem_pixels, desc="Processing pixels"))
redmapper_mem_data = vstack(results)
matched_member_path = os.path.join(base_path,f'chunhao-redmapper{suffix}_mem_data_all_new.fits')
redmapper_mem_data.write(matched_member_path, overwrite=True)
