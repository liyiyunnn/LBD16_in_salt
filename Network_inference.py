## adapted from http://arboreto.readthedocs.io.

import sys
sys.path.append('/Library/anaconda/anaconda3/lib/python3.7/site-packages/')

# import python modules
import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

if __name__ == '__main__':
    # load the data
    ex_matrix = pd.read_csv("/network-grn/input_saltdataset_ste.txt", sep='\t')
    tf_names = load_tf_names("/network-grn/input_genelist.txt")

    # infer the gene regulatory network
    network = grnboost2(expression_data=ex_matrix)
                        tf_names=tf_names)
    network.to_csv('/network-grn/output_salt_grn_ste.txt', sep='\t', index=False, header=False)
    network.head()
