# Cluster-Halo Matching (CHM)
Users can directly provide (and modify as you like) **command-line arguments** for cluster member path, halo path, output path, and keys (**in a specific order**) in the 'Main/CHM.py' using `argparse`.
## Inputs
* Halos: halo's 3D positions, ID, mass, and radius in `fits` or `hdf5` format (See example for saving Cardinal halos in `fits` format)
* Cluster members: members' 3D positions, cluster ID, and galaxy ID in `fits` or `hdf5` format (See example)
## Outputs
Detailed descriptions can be found at: https://github.com/ShuleiCao/CHM/wiki

## Dependencies
The following modules are required:
* numpy
* astropy
* h5py
* scipy
* multiprocessing
* joblib
* tqdm (optional, can be removed)
* pandas (optional)
* healpy (optional) 
