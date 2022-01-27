import argparse
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse import coo_matrix

parser = argparse.ArgumentParser()
parser.add_argument("-sample", "--sample_type", type=str, help="sample models used")
args = parser.parse_args()

def weighted_2_binary(args):
    sample_type = args.sample_type
    print('processing type ' + sample_type)
    k_fold_num = 10
    replicates_num = 20

    base_dir = os.getcwd()
    tf = pd.read_csv(base_dir + "/Adjacency_matrix/{}.TF.list.txt".format(sample_type), header=None)
    len_ = tf.shape[0]

    for i in range(1, replicates_num+1):
        result_tmp = np.zeros((len_,len_))
        test_accuracy_for_given_replicate = 0
        data_dir = "/{}/train_info/data_k_{}_r_{}".format(sample_type,str(k_fold_num), str(i))
        result_dir = "/{}/train_info/result_k_{}_r_{}".format(sample_type,str(k_fold_num), str(i))

        for j in range(1, k_fold_num+1):
            coo = np.array(pd.read_csv(base_dir + result_dir +
                             "/adj_matrix_predicted_weighted_using_RP_feature_"+str(j)+".txt", sep=" ", header=None))
            coo_new = coo_matrix((coo[:, 2], (coo[:, 0] - 1, coo[:, 1] - 1)), shape=(len_, len_))
            matrix = np.array(coo_new.todense())

            neg_ = np.array(pd.read_csv(base_dir + data_dir +
                                        "/test_negative_edge_coordinate_k_"+str(j)+".txt", sep="\t", header=None))
            pos_ = np.array(pd.read_csv(base_dir + data_dir +
                                        "/test_positive_edge_coordinate_k_"+str(j)+".txt", sep="\t", header=None))

            neg_list = matrix[neg_[:, 0], neg_[:, 1]]
            pos_list = matrix[pos_[:, 0], pos_[:, 1]]
            test_list = np.hstack((neg_list, pos_list))
            thresh = np.median(test_list)

            matrix_binary = (matrix > thresh) + 0
            matrix_binary = pd.DataFrame(matrix_binary)

            result_tmp = result_tmp + matrix_binary

            final_tmp = pd.concat([tf.T, matrix_binary], axis=0)
            tf_new = pd.concat([pd.DataFrame(['gene']), tf], axis=0)
            final_matrix = pd.concat([tf_new, final_tmp], axis=1)

            final_matrix.to_csv(base_dir + result_dir +
                                "/adj_matrix_predicted_binary_k_" + str(j) + ".txt",
                                float_format='%d', sep='\t', index=False, header=None)
            print("     {} KFold = {} :  replicate index = {}, k fold index = {} prediction determined".format(
                sample_type,str(k_fold_num), str(i), str(j)))

        final_result_tmp = pd.concat([tf.T,  pd.DataFrame(result_tmp)], axis=0)
        final_result = pd.concat([tf_new,final_result_tmp], axis=1)
        final_result.to_csv(base_dir + result_dir +
                            "/adj_matrix_predicted_binary_sum.txt",
                            float_format='%d', sep='\t', index=False, header=None)

        result_tmp_2 = (result_tmp >= int(k_fold_num*0.6) ) + 0
        final_result_tmp_2 = pd.concat([tf.T, pd.DataFrame(result_tmp_2)], axis=0)
        final_result_2 = pd.concat([tf_new, final_result_tmp_2], axis=1)
        final_result_2.to_csv(base_dir + result_dir +
                            "/adj_matrix_predicted_binary_final.txt",
                            float_format='%d', sep='\t', index=False, header=None)
        print("{} KFold = {} : replicate index = {} FINAL prediction determined".format(
                sample_type,str(k_fold_num), str(i)))

if __name__ == '__main__':
    weighted_2_binary(args)
