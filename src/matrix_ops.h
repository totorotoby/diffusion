#include <stdbool.h>

int csr_add(csr *A, csr *B, double scalar, csr **result);
int spmv(csr *A, double *y, double *x);
bool is_symmetric_csr(csr *A);
void csr_to_csc(int num_rows, int num_cols,
		const int *row_ptr, const int *cols, const double *vals,
                int **col_ptr, int **row_indices, double **values);
