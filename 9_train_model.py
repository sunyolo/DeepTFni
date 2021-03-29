from __future__ import division
from __future__ import print_function

from __future__ import division
from __future__ import print_function

import argparse
import time
import logging
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from torch import optim

import os
import inspect
from model import GCNModelVAE
# print (inspect.getsourcelines(GCNModelVAE))
from optimizer import loss_function
from torch.optim.lr_scheduler import StepLR, ReduceLROnPlateau
from utils import load_data, mask_test_edges, preprocess_graph, get_roc_score
from pytorchtools import EarlyStopping
import matplotlib
matplotlib.use('Agg')

parser = argparse.ArgumentParser()  
parser.add_argument("-sample", "--sample_type", type=str, help="sample models used")
parser.add_argument('--model', type=str, default='gcn_vae', help="models used") 
parser.add_argument('--seed', type=int, default=42, help='Random seed.')
parser.add_argument('--epochs', type=int, default=2000, help='Number of epochs to train.')
parser.add_argument('--hidden1', type=int, default=32, help='Number of units in hidden layer 1.')
parser.add_argument('--hidden2', type=int, default=16, help='Number of units in hidden layer 2.')
parser.add_argument('--lr', type=float, default=0.0001, help='Initial learning rate.')
parser.add_argument('--dropout', type=float, default=0., help='Dropout rate (1 - keep probability).')
parser.add_argument('--feature_type', type=str, default='RP', help='type of feature(TOY RP).')

args = parser.parse_args()
print(args)

def gae_for(args):

    print("start DeepTFni forward")
    feature_type = args.feature_type
    epochs = 2000

    sample_type = args.sample_type
    print('processing type ' + sample_type )

    k_fold_num = 10
    replicate_num = 20
    base_dir = os.getcwd()

    ##
    df_ = pd.read_csv(base_dir + "/Adjacency_matrix/{}.txt".format(sample_type), sep="\t")
    df = df_.values[:,1:]
    adj = sp.csr_matrix(df.astype(int))
    print(type(adj), adj.shape, adj.data.shape)

    ##
    df = pd.read_csv(base_dir + "/{}/10_{}_TF_rp_score.txt".format(sample_type,sample_type), sep="\t")
    df = df.values
    features = torch.FloatTensor(np.array(df))
    print(type(features), features.shape, features.data.shape)
    n_nodes, feat_dim = features.shape

    ##
    adj_orig = adj
    adj_orig = adj_orig - sp.dia_matrix((adj_orig.diagonal()[np.newaxis, :], [0]), shape=adj_orig.shape)
    adj_orig.eliminate_zeros()

    for m in range(1, replicate_num + 1):

        dir_train_info = base_dir + "/{}/train_info/data_k_{}_r_{}".format(sample_type,str(k_fold_num),str(m))
        df = pd.read_csv(dir_train_info + '/all_positive_edges.txt', sep="\t", header=None)
        positive_edges_all = np.array(df.values)
        df = pd.read_csv(dir_train_info + '/all_negative_edges.txt', sep="\t", header=None)
        negative_edges_all = np.array(df.values)

        dir_train_result = base_dir + "/{}/train_info/result_k_{}_r_{}".format(sample_type,str(k_fold_num),str(m))
        if not os.path.exists(dir_train_result):
            os.makedirs(dir_train_result)
        #########################################################################################
        for k_th in range(1, k_fold_num + 1):
            df = pd.read_csv(dir_train_info + '/adj_train_k_{}.txt'.format(str(k_th)), sep="\t", header=None)
            adj_train = sp.csr_matrix(df.values)
            df = pd.read_csv(dir_train_info + '/adj_val_k_{}.txt'.format(str(k_th)), sep="\t", header=None)
            adj_val = sp.csr_matrix(df.values)
            df = pd.read_csv(dir_train_info + '/val_positive_edge_coordinate_k_{}.txt'.format(str(k_th)), sep="\t",
                             header=None)
            val_edges = np.array(df.values)
            df = pd.read_csv(dir_train_info + '/val_negative_edge_coordinate_k_{}.txt'.format(str(k_th)), sep="\t",
                             header=None)
            val_edges_false = np.array(df.values)
            df = pd.read_csv(dir_train_info + '/test_positive_edge_coordinate_k_{}.txt'.format(str(k_th)), sep="\t",
                             header=None)
            test_edges = np.array(df.values)
            df = pd.read_csv(dir_train_info + '/test_negative_edge_coordinate_k_{}.txt'.format(str(k_th)), sep="\t",
                             header=None)
            test_edges_false = np.array(df.values)

            adj = adj_train
            adj_norm_train = preprocess_graph(adj)
            adj_label_train = adj_train + sp.eye(adj_train.shape[0])
            adj_label_train = torch.FloatTensor(adj_label_train.toarray())
            pos_weight_train = torch.Tensor([float(adj.shape[0] * adj.shape[0] - adj.sum()) / adj.sum()])
            norm_train = adj.shape[0] * adj.shape[0] / float((adj.shape[0] * adj.shape[0] - adj.sum()) * 2)

            adj = adj_val
            adj_norm_val = preprocess_graph(adj)
            adj_label_val = adj_val + sp.eye(adj_val.shape[0])
            adj_label_val = torch.FloatTensor(adj_label_val.toarray())
            pos_weight_val = torch.Tensor([float(adj.shape[0] * adj.shape[0] - adj.sum()) / adj.sum()])
            norm_val = adj.shape[0] * adj.shape[0] / float((adj.shape[0] * adj.shape[0] - adj.sum()) * 2)


            model = GCNModelVAE(feat_dim, args.hidden1, args.hidden2, args.dropout)
            print("for the " + str(k_th) + " model training")  # print(model)
            optimizer = optim.Adam(model.parameters(), lr=args.lr)
            scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=1000,
                                          verbose=True, threshold=1e-4, threshold_mode='rel',
                                          cooldown=0, min_lr=0, eps=1e-8)

            early_stopping = EarlyStopping(patience=100, verbose=True,
                                           path=dir_train_result + '/checkpoint_k_{}.pt'.format(str(k_th)))

            ####################################################################################################
            hidden_emb = None
            train_loss_record = []
            test_auc = []
            test_ap = []
            test_accuracy_for_record = []
            label_test = np.hstack([np.ones(len(test_edges)), np.zeros(len(test_edges_false))])

            for epoch in range(epochs):

                t = time.time()
                model.train()
                optimizer.zero_grad()
                recovered_train, mu_train, logvar_train = model(features, adj_norm_train)  #
                train_loss = loss_function(preds=recovered_train, labels=adj_label_train,
                                           mu=mu_train, logvar=logvar_train, n_nodes=n_nodes,
                                           norm=norm_train, pos_weight=pos_weight_train)
                train_loss.backward()
                cur_loss = train_loss.item()
                optimizer.step()

                recovered_val, mu_val, logvar_val = model(features, adj_norm_val)
                val_loss = loss_function(preds=recovered_val, labels=adj_label_val,
                                         mu=mu_val, logvar=logvar_val, n_nodes=n_nodes,
                                         norm=norm_val, pos_weight=pos_weight_val)

                scheduler.step(val_loss)
                early_stopping(val_loss, model)

                if early_stopping.early_stop:
                    print("EARLY STOPPING")
                    break

                hidden_emb = mu_train.data.numpy()
                fpr_curr, tpr_curr, thresholds, auc_curr, roc_curr, ap_curr, preds_curr = get_roc_score(hidden_emb,
                                                                                                        adj_orig,
                                                                                                        test_edges,
                                                                                                        test_edges_false)
                train_loss_record.append(cur_loss)
                test_auc.append(auc_curr)
                test_ap.append(ap_curr)

                if ((epoch + 1) % 20 == 0):
                    thresh_ = np.median(preds_curr)
                    preds_curr_binary = (preds_curr >= thresh_) + 0
                    temp = preds_curr_binary + label_test
                    test_accuracy_curr = (sum(temp == 2) + sum(temp == 0)) / len(temp)
                    test_accuracy_for_record.append(test_accuracy_curr)
                    print("\n##   Epoch:", '%04d' % (epoch + 1),
                          "test_accuracy=", "{:.5f}\n".format(test_accuracy_curr)
                          )

                print("Epoch:", '%04d' % (epoch + 1), "train_loss=", "{:.5f}".format(cur_loss),
                      "lr=", "{:.5f}".format(optimizer.param_groups[0]['lr']),
                      "time=", "{:.5f}".format(time.time() - t)
                      )

            model.load_state_dict(torch.load(dir_train_result + '/checkpoint_k_{}.pt'.format(str(k_th))))
            recovered_train, mu_train, logvar_train = model(features, adj_norm_train)
            hidden_emb = mu_train.data.numpy()
            print("Optimization Finished for the " + str(k_th) + "th model!")

            np.savetxt(dir_train_result + '/train_loss_using_' + feature_type + '_feature_' + str(k_th) + '.txt',
                       train_loss_record, fmt='%.6f', delimiter='\t')
            np.savetxt(dir_train_result + '/test_accuracy_using_' + feature_type + '_feature_' + str(k_th) + '.txt',
                       test_accuracy_for_record, fmt='%.6f', delimiter='\t')

            fpr_curr_all, tpr_curr_all, thresholds, auc_curr_all, roc_curr_all, ap_curr_for_edges_all, preds_curr_for_edges_all = get_roc_score(
                hidden_emb, adj_orig, positive_edges_all, negative_edges_all)

            edges_all = np.vstack((positive_edges_all, negative_edges_all))
            edges_all_with_score = np.vstack(
                (edges_all[:, 0] + 1, edges_all[:, 1] + 1, preds_curr_for_edges_all)).T
            np.savetxt(
                dir_train_result + '/adj_matrix_predicted_weighted_using_' + feature_type + '_feature_' + str(
                    k_th) + '.txt',
                edges_all_with_score, fmt='%d %d %.6f', delimiter='\t')

print("Finally finished!")

# """

if __name__ == '__main__':
    gae_for(args)