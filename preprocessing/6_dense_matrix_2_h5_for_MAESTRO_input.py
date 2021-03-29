import pandas as pd
import numpy as np
import os
import h5py
import numpy
import numpy as np
import scipy.io
import gzip
import scipy.sparse as sp_sparse
import pandas as pd
from scipy.sparse import coo_matrix

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-sample", "--sample_type",type=str,help="the name of input")
parser.add_argument("-input_file", "--input_dir_file",type=str,help="the path of input")
parser.add_argument("-output_dir", "--output_directory",type=str,help="the path of h5 file output")

args = parser.parse_args()

def mtxtoh5_parser(subparsers):
    """

    Add main function mtx-to-h5 argument parsers.

    """

    workflow = subparsers.add_parser("mtx-to-h5",

                                     help="Convert 10X mtx format matrix to HDF5 format.")

    group_input = workflow.add_argument_group("Input files arguments")

    group_input.add_argument("--type", dest="datatype", default="Peak",

                             choices=["Peak", "Gene"],

                             help="Type of the count matrix (Peak for scATAC-seq and Gene for scRNA-seq). DEFAULT: Peak.")

    group_input.add_argument("--matrix", dest="matrix", default="matrix.mtx",

                             help="Location of .mtx formatted matrix file. DEFAULT: matrix.mtx.")

    group_input.add_argument("--feature", dest="feature", default="features.tsv",

                             help="Location of feature file (required for the format of 'mtx'). Features correspond to row indices of count matrix. "

                                  "If the type is Peak, please provide the peak bed file with 3 columns. "

                                  "If the type is Gene, each row should include gene name. DEFAULT: features.tsv.")

    group_input.add_argument("--gene-column", dest="gene_column", default=2, type=int,

                             help="If the type is 'Peak', please specify which column of the feature file to use for gene names. DEFAULT: 2.")

    group_input.add_argument("--barcode", dest="barcode", default="barcodes.tsv",

                             help="Location of barcode file. Cell barcodes correspond to column indices of count matrix. DEFAULT: barcodes.tsv. ")

    group_input.add_argument("--species", dest="species", default="GRCh38",

                             choices=["GRCh38", "GRCm38"], type=str,

                             help="Species (GRCh38 for human and GRCm38 for mouse). DEFAULT: GRCh38.")

    group_output = workflow.add_argument_group("Output arguments")

    group_output.add_argument("-d", "--directory", dest="directory", default="MAESTRO",

                              help="Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")

    group_output.add_argument("--outprefix", dest="outprefix", default="MAESTRO",

                              help="Prefix of output files. DEFAULT: MAESTRO.")


def write_10X_h5(filename, matrix, features, barcodes, genome='GRCh38', datatype='Peak'):
    """Write 10X HDF5 files, support both gene expression and peaks."""

    f = h5py.File(filename, 'w')

    if datatype == 'Peak':

        M = sp_sparse.csc_matrix(matrix, dtype=numpy.int8)

    else:

        M = sp_sparse.csc_matrix(matrix, dtype=numpy.float32)

    B = numpy.array(barcodes, dtype='|S200')

    P = numpy.array(features, dtype='|S100')

    GM = numpy.array([genome] * len(features), dtype='|S10')

    FT = numpy.array([datatype] * len(features), dtype='|S100')

    AT = numpy.array(['genome'], dtype='|S10')

    mat = f.create_group('matrix')

    mat.create_dataset('barcodes', data=B)

    mat.create_dataset('data', data=M.data)

    mat.create_dataset('indices', data=M.indices)

    mat.create_dataset('indptr', data=M.indptr)

    mat.create_dataset('shape', data=M.shape)

    fet = mat.create_group('features')

    fet.create_dataset('_all_tag_keys', data=AT)

    fet.create_dataset('feature_type', data=FT)

    fet.create_dataset('genome', data=GM)

    fet.create_dataset('id', data=P)

    fet.create_dataset('name', data=P)

    f.close()



def read_10X_mtx(matrix_file, feature_file, barcode_file, datatype, gene_column=2):
    """Convert 10x mtx as matrix."""

    print('read_10X_mtx function info')
    
    print('\tprocessing feature_file')
    if feature_file.split('.')[-1] == 'gz' or feature_file.split('.')[-1] == 'gzip':

        feature_in = gzip.open(feature_file, "r")

    else:

        feature_in = open(feature_file, "r")

    features = feature_in.readlines()

    if datatype == "Peak":

        features = ["_".join(feature.strip().split("\t")[0:3]) for feature in features]

    else:

        if type(features[0]) == str:
            features = [feature.strip().split("\t")[gene_column - 1] for feature in features]

        if type(features[0]) == bytes:
            features = [feature.decode().strip().split("\t")[gene_column - 1] for feature in features]
    print('\texample of feature[0] : '+ str(features[0]))
    print('\n')

    print('\tprocessing barcode_file')
    if barcode_file.split('.')[-1] == 'gz' or barcode_file.split('.')[-1] == 'gzip':

        barcode_in = gzip.open(barcode_file, "r")

    else:

        barcode_in = open(barcode_file, "r")

    barcodes = barcode_in.readlines()

    if type(barcodes[0]) == str:
        barcodes = [barcode.strip().split("\t")[0] for barcode in barcodes]

    if type(barcodes[0]) == bytes:
        barcodes = [barcode.decode().strip().split("\t")[0] for barcode in barcodes]
    print('\texample of barcode[0] : '+ str(barcodes[0]))
    print('\n')

    print('\tprocessing matrix_file')
    matrix = scipy.io.mmread(matrix_file)

    matrix = sp_sparse.csc_matrix(matrix, dtype=numpy.float32)

    return {"matrix": matrix, "features": features, "barcodes": barcodes}


def mtx_2_h5(directory, outprefix, matrix_file, feature_file, barcode_file, gene_column=2, genome='GRCh38',
             datatype='Peak'):
    """Convert 10x mtx format matrix to HDF5."""

    try:

        os.makedirs(directory)

    except OSError:

        # either directory exists (then we can ignore) or it will fail in the

        # next step.

        pass

    if datatype == "Peak":

        filename = os.path.join(directory, outprefix + "_peak_count.h5")

    else:

        filename = os.path.join(directory, outprefix + "_gene_count.h5")

    print('mtx_2_h5 function info')
    print('\ttarget_file : ' +filename)
    print('\tmatrix_file : ' + matrix_file)
    print('\tfeature_file : ' + feature_file)
    print('\tbarcode_file : ' + barcode_file)
    print('\n')

    matrix_dict = read_10X_mtx(matrix_file=matrix_file, feature_file=feature_file, barcode_file=barcode_file, datatype=datatype, gene_column=gene_column)

    write_10X_h5(filename=filename, matrix=matrix_dict["matrix"],features=matrix_dict["features"], barcodes=matrix_dict["barcodes"], genome=genome, datatype=datatype)


sample_type = args.sample_type
input_dir_file = args.input_dir_file

data = pd.read_csv(input_dir_file, sep=',')
data_sparse = coo_matrix(data)

coordinate_with_count = np.array([data_sparse.row+1, data_sparse.col+1, data_sparse.data])
Num = len(data_sparse.data)
P_num = np.size(data.values,0)
C_num = np.size(data.values,1)

input_directory = "./"	
dir_1 = input_directory + sample_type
if not(os.path.exists(dir_1)):
	os.mkdir(dir_1)
mtx_file_dir = dir_1 + "/mtx_file"
if not(os.path.exists(mtx_file_dir)):
	os.mkdir(mtx_file_dir)

np.savetxt(mtx_file_dir + "/feature.txt", data.axes[0], fmt='%s')
np.savetxt(mtx_file_dir +  "/barcodes.txt", data.axes[1], fmt='%s')
print()

with open(mtx_file_dir + "/" + sample_type + ".mtx", 'a') as f:
    f.write("%%MatrixMarket matrix coordinate integer general\n")
    f.write("%\n")
    f.write(str(P_num) + " " + str(C_num) + " " + str(Num) + "\n")
    np.savetxt(f, coordinate_with_count.transpose(), fmt='%d', delimiter=' ')

mtx_file_dir = "./"+ sample_type + "/mtx_file"
matrix_file = os.path.join(mtx_file_dir + "/" + sample_type + ".mtx")
feature_file = os.path.join(mtx_file_dir + "/feature.txt")
barcode_file = os.path.join(mtx_file_dir + "/barcodes.txt")

output_directory = args.output_directory
outprefix = "8_"+sample_type	

mtx_2_h5(directory=output_directory, outprefix=outprefix, matrix_file=matrix_file, feature_file=feature_file, barcode_file=barcode_file, gene_column=2, genome='GRCh38',
         datatype='Peak')


