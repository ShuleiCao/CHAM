# Cluster-Halo Matching (CHM)
## Inputs
* Halos: halo ID, mass, radius, and 3D true positions in `fits` or `hdf5` format (See example for preparing Cardinal halos in `fits` format)
* Cluster members: cluster ID, and galaxy ID and 3D positions in `fits` or `hdf5` format (See example for preparing Cardinal redMaPPer members in `fits` format)
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
