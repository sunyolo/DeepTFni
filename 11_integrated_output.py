import argparse
import numpy as np
import pandas as pd
from pandas import Series
import scipy.sparse as sp
import os

parser = argparse.ArgumentParser()
parser.add_argument("-sample", "--sample_type", type=str, help="sample models used")
args = parser.parse_args()
#print(args)

def integrated_output(args):
    sample_type = args.sample_type
    print('processing type ' + sample_type)
    k_fold_num = 10
    replicates_num = 20

    base_dir = os.getcwd()
    summary_dir = base_dir + "/{}/train_info/result_k_{}_summary".format(sample_type,str(k_fold_num))
    if not os.path.exists(summary_dir):
        os.makedirs(summary_dir)

    tf = pd.read_csv(base_dir + "/Adjacency_matrix/{}.TF.list.txt".format(sample_type), header=None)
    len_ = tf.shape[0]

    result_tmp = np.zeros((len_, len_))
    for i in range(1, replicates_num+1):

        result_dir = base_dir + "/{}/train_info/result_k_{}_r_{}".format(sample_type,str(k_fold_num),str(i))

        result_for_curr_replicate = np.array(pd.read_csv(result_dir +
                         "/adj_matrix_predicted_binary_final.txt", sep="\t", header=None).values[1:,1:].astype(int))

        result_tmp = result_tmp + result_for_curr_replicate

    shared_TF_link_array = (result_tmp >= replicates_num - 2) + 0
    tmp = pd.concat([tf.T, pd.DataFrame(shared_TF_link_array)], axis=0)
    TF_list_new = pd.concat([pd.DataFrame(['TF']), tf], axis=0)
    shared_TF_link_array_with_TF = pd.concat([TF_list_new, tmp], axis=1)

    shared_TF_link_array_with_TF.to_csv(summary_dir + '/adj_matrix_predicted_binary_final_integrated_from_20_replicates.txt', float_format='%d',
                                  sep='\t', index=False, header=None)

    np.savetxt(summary_dir + "/adj_matrix_predicted_sum_final_integrated_from_20_replicates.txt",
               result_tmp, fmt='%d', delimiter='\t')

if __name__ == '__main__':
    integrated_output(args)
