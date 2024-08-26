import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.linalg import cholesky
import argparse

def read_file(file_path):
    with open(file_path, 'r') as file:
        return [int(line.strip()) for line in file]

def plot_sparsity_pattern(csr):
    plt.spy(csr, markersize=1)
    plt.title("Sparsity Pattern")
    plt.show()

def check_symmetric(csr):
    mat = (csr - csr.transpose())
    for row in range(mat.shape[0] - 1):
        begin = mat.indptr[row]
        end = mat.indptr[row+1]
        for i in range(begin, end):
            print("(" + str(row) + "," + str( mat.indices[i]) + "," + str(mat.data[i]) + ")")
    exit(0)

def main(row_ptr_file, columns_file, values_file):

    row_ptr = read_file(row_ptr_file)
    columns = read_file(columns_file)
    values = np.loadtxt(values_file)

    num_rows = len(row_ptr) - 1
    csr = csr_matrix((values, columns, row_ptr), shape=(num_rows, num_rows))

    if check_symmetric(csr):
        print("The matrix is symmetric.")
    else:
        print("The matrix is not symmetric.")

    plot_sparsity_pattern(csr)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot the sparsity pattern of a CSR matrix from files, and check symmetry.")
    parser.add_argument("row_ptr_file", help="Path to the row pointer file")
    parser.add_argument("columns_file", help="Path to the columns file")
    parser.add_argument("values_file", help="Path to the values file")

    args = parser.parse_args()

    main(args.row_ptr_file, args.columns_file, args.values_file)
