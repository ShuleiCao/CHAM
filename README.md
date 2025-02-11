# Cluster-Halo Matching (CHM)
Input (halo ID, mass, radius, and 3D true positions saved in `fits` or `hdf5` format) with redMaPPer members (with cluster ID, and galaxy ID and 3D positions saved in `fits` or `hdf5` format).

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

## Outputs explanation
1. Clusters are stored as `hdf5` files, one file per cluster, named by cluster ID (e.g., `1.hdf5` for a cluster with ID 1). Each file contains `associated` and `unassociated` categories, where `associated` has members associated with halos and `unassociated` has members not belonging to any halo. The `associated` category organizes data by halo ID.  Each halo's subgroup contains datasets of galaxy properties for associated members and attributes for the number of members and halo mass. The `unassociated` catagory has similar structure but with halo ID assigned to galaxies in the original (`gold` for Cardinal) catalog.
2. The clusters (along with `centrals` from cluster catalog) will be saved into one `hdf5` file in different catagories organized by cluster ID.
