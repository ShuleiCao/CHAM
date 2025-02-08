## Workflow for Matching Halos with redMaPPer Members

1. **Save Host Halos (Optional)**:  
   If host (main) halos (matched with `gold` and `bpz` catalogs for Cardinal) are not available, run `host_halo_save.py` to create and save them in `fits` format.

2. **Save redMaPPer Members (Optional)**:  
   Run `redmapper_member_save.py` to process and save redMaPPer members (also matched with `gold` and `bpz` catalogs for Cardinal) in `fits` format.

3. **Match Halos with redMaPPer Members**:  
  Run `match_redmapper_with_gold_process` to match the halos (with halo ID, mass, radius, and 3D true positions saved in `fits` or `hdf5` format) with redMaPPer members (with cluster ID, and galaxy ID and 3D positions saved in `fits` or `hdf5` format). Clusters are stored as `hdf5` files, one file per cluster, named by cluster ID (e.g., `1.hdf5` for a cluster with ID 1). Each file contains `associated` and `unassociated` categories, where `associated` has members associated with halos and `unassociated` has members not belonging to any halo.

4. **Consolidate**:  
   Consolidate all clusters into one `hdf5` file for post analyses.

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
