### Workflow for Matching Halos with redMaPPer Members

1. **Save Host Halos (Optional)**:  
   If host halos (matched with `gold` and `bpz` catalogs) are not available, run `host_halo_save.py` to create and save them.

2. **Save redMaPPer Members (Optional)**:  
   Run `redmapper_member_save.py` to process and save redMaPPer members (also matched with `gold` and `bpz` catalogs).

3. **Match Halos with redMaPPer Members**:  
   After saving both, choose one of `match_redmapper_with_gold` files (`_process` is better) to match the halos with redMaPPer members, depending on available memory.
