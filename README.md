### Workflow for Matching Halos with redMaPPer Members

1. **Save Host Halos (Optional)**:  
   If host halos (matched with `gold` and `bpz` catalogs for Cardinal) are not available, run `host_halo_save.py` to create and save them.

2. **Save redMaPPer Members (Optional)**:  
   Run `redmapper_member_save.py` to process and save redMaPPer members (also matched with `gold` and `bpz` catalogs for Cardinal).

3. **Match Halos with redMaPPer Members**:  
  Run `match_redmapper_with_gold_process` to match the halos with redMaPPer members.

4. **Consolidate**:  
   Consolidate all clusters into one data file for post analyses.
