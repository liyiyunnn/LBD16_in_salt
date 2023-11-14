# Author: Yiyun Li
# Affiliation: Wageningen University & Research

for i in {1..100}; do python /Network_inference_GRNBOOST2.py
"/LBD_100run/all_layer_salt/run${i}.txt" "/LBD_100run/all_layer_salt/run_anno${i}.txt"
