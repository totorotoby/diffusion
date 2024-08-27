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

def is_symmetric_csr(matrix: csr_matrix, tol=1e-10) -> bool:
    
    if matrix.shape[0] != matrix.shape[1]:
        return False

    for i in range(matrix.shape[0]):
        row_start = matrix.indptr[i]
        row_end = matrix.indptr[i + 1]
        for idx in range(row_start, row_end):
            j = matrix.indices[idx]
            a_ij = matrix.data[idx]
            
            row_j_start = matrix.indptr[j]
            row_j_end = matrix.indptr[j + 1]
            found_symmetric = False
            for k in range(row_j_start, row_j_end):
                if matrix.indices[k] == i:
                    a_ji = matrix.data[k]
                    if np.abs(a_ij - a_ji) < tol:
                        found_symmetric = True
                        break
                    else:
                        print(i, j, a_ij, a_ji)
                    
            if not found_symmetric:
                return False
    
    return True

def main(row_ptr_file, columns_file, values_file):

    row_ptr = read_file(row_ptr_file)
    columns = read_file(columns_file)
    values = np.loadtxt(values_file)

    num_rows = len(row_ptr) - 1
    csr = csr_matrix((values, columns, row_ptr), shape=(num_rows, num_rows))

    if is_symmetric_csr(csr):
        print("The matrix is symmetric.")
    else:
        print("The matrix is not symmetric.")

    eigenvalues = eigsh(csr, k=1, which='SA', return_eigenvectors=False)
    pos_def = np.all(eigenvalues > 0)
    print(eigenvalues)
    print(pos_def)
    
        
    plot_sparsity_pattern(csr)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot the sparsity pattern of a CSR matrix from files, and check symmetry.")
    parser.add_argument("row_ptr_file", help="Path to the row pointer file")
    parser.add_argument("columns_file", help="Path to the columns file")
    parser.add_argument("values_file", help="Path to the values file")

    args = parser.parse_args()

    main(args.row_ptr_file, args.columns_file, args.values_file)
