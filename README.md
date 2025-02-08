### Workflow for Matching Halos with redMaPPer Members

1. **Save Host Halos (Optional)**:  
   If host (main) halos (matched with `gold` and `bpz` catalogs for Cardinal) are not available, run `host_halo_save.py` to create and save them in `fits` format.

2. **Save redMaPPer Members (Optional)**:  
   Run `redmapper_member_save.py` to process and save redMaPPer members (also matched with `gold` and `bpz` catalogs for Cardinal) in `fits` format.

3. **Match Halos with redMaPPer Members**:  
  Run `match_redmapper_with_gold_process` to match the halos ([`m200`, `r200`, `haloid`, `px`, `py`, `pz`] saved in `fits` or `hdf5` format) with redMaPPer members (in `fits` or `hdf5` format). Clusters with `associated` halos

4. **Consolidate**:  
   Consolidate all clusters into one data file for post analyses.
