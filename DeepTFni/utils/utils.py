import pickle as pkl

import networkx as nx
import numpy as np
import scipy.sparse as sp
import torch
import copy
from sklearn.metrics import roc_auc_score, average_precision_score,roc_curve,auc

def parse_index_file(filename):
    index = []
    for line in open(filename):
        index.append(int(line.strip()))
    return index


def sparse_to_tuple(sparse_mx):
    if not sp.isspmatrix_coo(sparse_mx):
        sparse_mx = sparse_mx.tocoo()
    coords = np.vstack((sparse_mx.row, sparse_mx.col)).transpose()
    values = sparse_mx.data
    shape = sparse_mx.shape
    return coords, values, shape


def mask_test_edges(adj):
    # Remove diagonal elements
    adj = adj - sp.dia_matrix((adj.diagonal()[np.newaxis, :], [0]), shape=adj.shape)
    adj_dense = adj.todense()
    adj.eliminate_zeros()
    # Check that diag is zero:
    assert np.diag(adj.todense()).sum() == 0

    adj_triu = sp.triu(adj)
    adj_tuple = sparse_to_tuple(adj_triu)
    edges = adj_tuple[0]
    edges_all = sparse_to_tuple(adj)[0]
    num_test = int(np.floor(edges.shape[0] * 1 / 20.))
    num_val = int(np.floor(edges.shape[0] * 1 / 3.))

    all_edge_idx = list(range(edges.shape[0]))
    np.random.shuffle(all_edge_idx)
    val_edge_idx = all_edge_idx[:num_val]
    test_edge_idx = all_edge_idx[num_val:(num_val + num_test)]
    test_edges = edges[test_edge_idx]
    val_edges = edges[val_edge_idx]
    train_edges = np.delete(edges, np.hstack([test_edge_idx, val_edge_idx]), axis=0)



    def ismember(a, b, tol=5):
        rows_close = np.all(np.round(a - b[:, None], tol) == 0, axis=-1)
        return np.any(rows_close)

    test_edges_false = []
    while len(test_edges_false) < len(test_edges):
        idx_i = np.random.randint(0, adj.shape[0])
        idx_j = np.random.randint(0, adj.shape[0])
        if idx_i == idx_j:
            continue
        if ismember([idx_i, idx_j], edges_all):
            continue
        if test_edges_false:
            if ismember([idx_j, idx_i], np.array(test_edges_false)):
                continue
            if ismember([idx_i, idx_j], np.array(test_edges_false)):
                continue
        test_edges_false.append([idx_i, idx_j])

    val_edges_false = []
    while len(val_edges_false) < len(val_edges):
        idx_i = np.random.randint(0, adj.shape[0])
        idx_j = np.random.randint(0, adj.shape[0])
        if idx_i > idx_j:  #added by sunyu
            continue
        if idx_i == idx_j:
            continue
        if ismember([idx_i, idx_j], train_edges):
            continue
        if ismember([idx_j, idx_i], train_edges):
            continue
        if ismember([idx_i, idx_j], test_edges):
            continue
        if ismember([idx_j, idx_i], test_edges):
            continue
        if ismember([idx_i, idx_j], val_edges):
            continue
        if ismember([idx_j, idx_i], val_edges):
            continue
        if val_edges_false:
            if ismember([idx_j, idx_i], np.array(val_edges_false)):
                continue
            #if ismember([idx_i, idx_j], np.array(val_edges_false)): #
                #continue
        val_edges_false.append([idx_i, idx_j])
    
    assert ~ismember(test_edges_false, edges_all)
    #assert ~ismember(val_edges_false, edges_all) 
    assert ~ismember(val_edges, train_edges)
    assert ~ismember(test_edges, train_edges)
    assert ~ismember(val_edges, test_edges)

    data = np.ones(train_edges.shape[0])
    # Re-build adj matrix
    adj_train = sp.csr_matrix((data, (train_edges[:, 0], train_edges[:, 1])), shape=adj.shape)
    adj_train = adj_train + adj_train.T

    return adj_train, train_edges, val_edges, val_edges_false, test_edges, test_edges_false


def my_mask_test_edges(adj):

    column_sum = np.sum(adj, axis=0)
    row_sum = np.sum(adj, axis=1)
    row_i = (row_sum[:] != 0)
    column_i = (column_sum[:] != 0)
    row_index = np.argwhere(row_i == True)
    column_index = np.argwhere(column_i == True)

    edges = len(row_index)
    num_val = int(np.floor(edges * 1 / 10.))


    np.random.shuffle(row_index)
    val_edge_row_idx = row_index[:num_val]
    train_edge_row_idx = row_index[num_val:]
    val_edge_row_adj = np.zeros(adj.shape)
    val_edge_row_adj[val_edge_row_idx] = adj[val_edge_row_idx]
    train_edge_row_adj = np.zeros(adj.shape)
    train_edge_row_adj[train_edge_row_idx] = adj[train_edge_row_idx]


    np.random.shuffle(column_index)
    val_edge_column_idx = column_index[:num_val]
    train_edge_column_idx = column_index[num_val:]
    val_edge_column_adj = np.zeros(adj.shape)
    val_edge_column_adj[val_edge_column_idx] = adj[val_edge_column_idx]
    train_edge_column_adj = np.zeros(adj.shape)
    train_edge_column_adj[:, train_edge_column_idx] = adj[:, train_edge_column_idx]

    #  merge row_adj and column_adj
    val_adj = val_edge_row_adj + val_edge_column_adj
    train_adj = train_edge_row_adj + train_edge_column_adj

    val_adj = val_adj + val_adj.T
    train_adj = train_adj + train_adj.T

    val_adj[val_adj > 0] = 1
    train_adj[train_adj > 0] = 1

    train_adj = sp.csr_matrix(train_adj)
    val_adj = sp.csr_matrix(val_adj)
    return train_adj, train_edge_row_idx, train_edge_column_idx, val_adj, val_edge_row_idx, val_edge_column_idx



def preprocess_graph(adj):
    adj = sp.coo_matrix(adj)
    adj_ = adj + sp.eye(adj.shape[0])
    rowsum = np.array(adj_.sum(1))
    degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -0.5).flatten())
    adj_normalized = adj_.dot(degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt).tocoo()
    # return sparse_to_tuple(adj_normalized)
    return sparse_mx_to_torch_sparse_tensor(adj_normalized)


def sparse_mx_to_torch_sparse_tensor(sparse_mx):
    sparse_mx = sparse_mx.tocoo().astype(np.float32)
    indices = torch.from_numpy(
        np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.int64))
    values = torch.from_numpy(sparse_mx.data)
    shape = torch.Size(sparse_mx.shape)
    return torch.sparse.FloatTensor(indices, values, shape)


def get_roc_score(emb, adj_orig, edges_pos, edges_neg):
    def sigmoid(x):
        return 1 / (1 + np.exp(-x))

    # Predict on test set of edges
    adj_rec = np.dot(emb, emb.T)
    preds = []
    pos = []
    for e in edges_pos:
        preds.append(sigmoid(adj_rec[e[0], e[1]]))
        pos.append(adj_orig[e[0], e[1]])

    preds_neg = []
    neg = []
    for e in edges_neg:
        preds_neg.append(sigmoid(adj_rec[e[0], e[1]]))
        neg.append(adj_orig[e[0], e[1]])

    preds_all = np.hstack([preds, preds_neg])

    labels_all = np.hstack([np.ones(len(edges_pos)), np.zeros(len(edges_neg))])
    fpr, tpr, thresholds = roc_curve(labels_all, preds_all)
    roc_auc = auc(fpr, tpr)
    roc_score = roc_auc_score(labels_all, preds_all)
    ap_score = average_precision_score(labels_all, preds_all)

    return fpr, tpr,thresholds, roc_auc, roc_score, ap_score, preds_all

def retain_edges(all, part):
    all = all.detach().numpy()
    tmp = copy.deepcopy(all)
    all_new = []
    for i in range(len(part)):
        all_new.append(tmp[part[i, 0], part[i, 1]])
    all_new = torch.from_numpy(np.array(all_new))
    return all_new

def remove_edges(all, part):
    all = all.detach().numpy()
    tmp = copy.deepcopy(all)
    for i in range(len(part)):
        tmp[part[i,0],part[i,1]]=np.nan
    all_new = tmp[~np.isnan(tmp)]
    all_new = torch.from_numpy(np.array(all_new))
    return all_new
