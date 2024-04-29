#### 3 contacts

rt_process_vs.py read --input_db output.db --hydrogen_bond :GLN:19:,:GLY:65:,:GLY:66:  --log H_19_65_66_log.txt --bookmark_name H_19_65_66 --output_all_poses
mkdir H_19_65_66
rt_process_vs.py read --input_db output.db --bookmark_name H_19_65_66 --export_sdf_path H_19_65_66

#### 4 contacts

rt_process_vs.py read --input_db output.db --hydrogen_bond :GLN:19:,:SER:61,:GLY:65:,:GLY:66:  --log H_19_61_65_66_log.txt --bookmark_name H_19_61_65_66 --output_all_poses
mkdir H_19_61_65_66
rt_process_vs.py read --input_db output.db --bookmark_name H_19_61_65_66 --export_sdf_path H_19_61_65_66
