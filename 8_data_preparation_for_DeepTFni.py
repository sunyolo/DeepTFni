#-*- coding: UTF-8 -*-  
from __future__ import division
from __future__ import print_function
import argparse
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.model_selection import KFold
import os
from utils import sparse_to_tuple

def ismember(a, b, tol=5):
    rows_close = np.all(np.round(a - b[:, None], tol) == 0, axis=-1)
    return np.any(rows_close)

def findByRow(mat, row):
    return np.where((mat == row).all(1))[0]

parser = argparse.ArgumentParser()
parser.add_argument("-sample","--sample_type", type=str, help="models used")
args = parser.parse_args()
sample_type = args.sample_type

replicates_num = 20
k_fold_num = 10
kf = KFold(n_splits=k_fold_num)

base_dir = os.getcwd()
input_adj_file = base_dir+"/Adjacency_matrix/{}.txt".format(sample_type)

# Read adjacency matrix
df_ = pd.read_csv(input_adj_file, sep="\t",index_col=0)
df = df_.values
adj = sp.csr_matrix(df.astype(int))

# Remove diagonal elements
adj = adj - sp.dia_matrix((adj.diagonal()[np.newaxis, :], [0]), shape=adj.shape)
adj.eliminate_zeros()
# Check that diag is zero:
assert np.diag(adj.todense()).sum() == 0

adj_triu = sp.triu(adj)
adj_tuple = sparse_to_tuple(adj_triu)
edges = adj_tuple[0]
edges_all = sparse_to_tuple(adj)[0]
a = adj.todense().A
b = np.ones(adj.shape)
b = b.astype(np.int64)
adj_xor = a ^ b
edges_all_false = sparse_to_tuple(sp.csr_matrix(adj_xor))[0]

dir_train = base_dir + '/' + sample_type + "/train_info"
if not os.path.exists(dir_train):
    os.mkdir(dir_train)

for m in range (1,replicates_num+1):

    print("processing {} cell for replicate {}".format(sample_type, str(m)))
    out_dir = base_dir + "/{}/train_info/data_k_10_r_".format(sample_type) +str(m)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    np.savetxt(out_dir + '/all_positive_edges.txt', edges_all, fmt='%d', delimiter='\t')
    np.savetxt(out_dir + '/all_negative_edges.txt', edges_all_false, fmt='%d', delimiter='\t')

    all_edge_idx = np.arange(0, edges.shape[0])
    np.random.shuffle(all_edge_idx)

    i = 0
    for train, test in kf.split(all_edge_idx):
        i+=1
        test_edges_idx = all_edge_idx[test]
        test_edges = edges[test_edges_idx]

        test_edges_false = []
        while len(test_edges_false) < len(test_edges):
            idx_i = np.random.randint(0, adj.shape[0])
            idx_j = np.random.randint(0, adj.shape[0])
            if idx_i >= idx_j:
                continue
            if ismember([idx_i, idx_j], edges_all):
                continue
            if test_edges_false:
                if ismember([idx_j, idx_i], np.array(test_edges_false)):
                    continue
                if ismember([idx_i, idx_j], np.array(test_edges_false)):
                    continue
            test_edges_false.append([idx_i, idx_j])

        test_edges_false = np.array(test_edges_false)
        np.savetxt(out_dir + '/test_positive_edge_coordinate_k_' + str(i) + '.txt', test_edges, fmt='%d',
                   delimiter='\t')
        np.savetxt(out_dir + '/test_negative_edge_coordinate_k_' + str(i) + '.txt', test_edges_false, fmt='%d',
                   delimiter='\t')

        data = np.ones(test_edges.shape[0])
        adj_test = sp.csr_matrix((data, (test_edges[:, 0], test_edges[:, 1])), shape=adj.shape)
        adj_test = adj_test + adj_test.T
        adj_test_false = sp.csr_matrix((data, (test_edges_false[:, 0], test_edges_false[:, 1])), shape=adj.shape)
        adj_test_false = adj_test_false + adj_test_false.T
        adj_test_postive_negative = adj_test + adj_test_false

        np.savetxt(out_dir + '/adj_matrix_for_test_postive_negative_edges_k_' + str(i) + '.txt',
                   adj_test_postive_negative.todense().A, fmt='%d',
                   delimiter='\t')
        df_ = pd.read_csv(input_adj_file, sep="\t", index_col=0)
        df = df_.values
        adj_curr = sp.csr_matrix(df.astype(int))
        
    ###################################################################
        # Remove diagonal elements
        adj_curr = adj_curr - sp.dia_matrix((adj_curr.diagonal()[np.newaxis, :], [0]),
                                  shape=adj_curr.shape)
        adj_curr.eliminate_zeros()
        # Check that diag is zero:
        assert np.diag(adj_curr.todense()).sum() == 0

        adj_triu_curr = sp.triu(adj_curr)
        adj_tuple_curr = sparse_to_tuple(adj_triu_curr)
        edges_curr = adj_tuple_curr[0]
        a = adj_curr.todense().A
        b = np.ones(adj_curr.shape)
        b = b.astype(np.int64)
        adj_xor_curr = a ^ b
        edges_false_curr = sparse_to_tuple(sp.triu(sp.csr_matrix(adj_xor_curr)))[0]

        test_pos_idx_to_delete=[]
        test_neg_idx_to_delete = []
        for k in range (test_edges.shape[0]):
            tp1 = findByRow(edges_curr, test_edges[k])
            tp2 = findByRow(edges_false_curr, test_edges_false[k])
            if (tp1):
                test_pos_idx_to_delete.append(tp1)
            if (tp2):
                test_neg_idx_to_delete.append(tp2)

        new_edges_curr = np.delete(edges_curr, test_pos_idx_to_delete,0)
        new_edges_false_curr = np.delete(edges_false_curr, test_neg_idx_to_delete,0)

        rest_edge_idx = np.arange(0, new_edges_curr.shape[0])
        np.random.shuffle(rest_edge_idx)
        train_num = len(rest_edge_idx) * 4 // 5

        train_idx = rest_edge_idx[0:train_num]
        val_idx =  rest_edge_idx[train_num:]
        train_edges =  new_edges_curr[train_idx]
        val_edges = new_edges_curr[val_idx]

        rest_false_edge_idx = np.arange(0, new_edges_false_curr.shape[0])
        np.random.shuffle(rest_false_edge_idx)

        rest_false_edge_idx = pd.DataFrame(rest_false_edge_idx)
        train_val_false_idx = np.array(rest_false_edge_idx.sample(n=len(rest_edge_idx),
                                                           frac=None, replace=False, axis=0))
        train_val_false_idx = train_val_false_idx.reshape(-1,)

        train_false_idx =  train_val_false_idx[0:train_num]
        val_false_idx = train_val_false_idx[train_num:]
        train_false_edges = new_edges_false_curr[train_false_idx]
        val_false_edges = new_edges_false_curr[val_false_idx]

        data = np.ones(train_edges.shape[0])
        adj_train = sp.csr_matrix((data, (train_edges[:, 0], train_edges[:, 1])), shape=adj.shape)
        adj_train = adj_train + adj_train.T

        data = np.ones(val_edges.shape[0])
        adj_val = sp.csr_matrix((data, (val_edges[:, 0], val_edges[:, 1])), shape=adj.shape)
        adj_val = adj_val + adj_val.T

        #test_edges_false = np.array(test_edges_false)
        np.savetxt(out_dir + '/adj_train_k_{}.txt'.format(str(i)),
                   adj_train.todense().A, fmt='%d',delimiter='\t')
        np.savetxt(out_dir + '/adj_val_k_{}.txt'.format(str(i)),
                   adj_val.todense().A, fmt='%d', delimiter='\t')
        np.savetxt(out_dir + '/train_positive_edge_coordinate_k_{}.txt'.format(str(i)),
                   train_edges,fmt='%d',delimiter='\t')
        np.savetxt(out_dir + '/train_negative_edge_coordinate_k_{}.txt'.format(str(i)),
                   train_false_edges,fmt='%d',delimiter='\t')
        np.savetxt(out_dir + '/val_positive_edge_coordinate_k_{}.txt'.format(str(i)),
                   val_edges, fmt='%d', delimiter='\t')
        np.savetxt(out_dir + '/val_negative_edge_coordinate_k_{}.txt'.format(str(i)),
                   val_false_edges, fmt='%d', delimiter='\t')
