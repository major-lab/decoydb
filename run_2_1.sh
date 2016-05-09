
python 2_1_mcflashfold_mask_generator.py \
--digested_data ../workdir_decoy/digested_data.pk \
--mcfold_cmd "../FlashFold/bin/flashfold -ft 1000 -s {seq} -n {name} -m {mask} -tables ../FlashFold/tables > ../workdir_decoy/flashfold_raw/{accession}"
