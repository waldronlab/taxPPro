
## Step 1 - Get prunde and reference trees, and trait table in the right format.
format_tree_and_trait_table.py -t ../../extdata/LTP_all_08_2023.newick -i trait.tab

## Step 1b - It seems that single quotes are added. Remove them.
sed -i -e "s/'//g" formatted/reference_tree.newick
sed -i -e "s/'//g" formatted/pruned_tree.newick

## Perform ancestral state reconstruction 
ancestral_state_reconstruction.py -i formatted/trait_table.tab -t formatted/pruned_tree.newick -o asr_counts.tab -c asr_ci.tab

## Create predictions
predict_traits.py -a -i formatted/trait_table.tab -t formatted/reference_tree.newick -r asr_counts.tab -c asr_ci.tab -o predict_traits.tab
