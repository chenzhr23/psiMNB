python PUSscan_build.py -training_file human_PUS_MNB_input_k-mer_overall.txt -model_name overall -to_predict Day0_common_anno_group_redundance_mix.txt -output_dir /public/home/chenzr/PSI_Seq_brainCell/A1-A12-totalRNA-result/psiFinder_ANN_res/PUSscan_test 

python PUSscan_predict.py -model_file overall_multinomialnb_model.pkl -to_predict Day0_common_anno_group_redundance_mix.txt -output_dir /public/home/chenzr/PSI_Seq_brainCell/A1-A12-totalRNA-result/psiFinder_ANN_res/PUSscan_test 
