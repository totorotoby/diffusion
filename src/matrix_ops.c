/* matrix_ops.c -- matrix operations, solvers, LU, spmv, etc.
   *
   * Written on Thursday, 15 August 2024.
  */
  
#include <stdio.h>
#include "matrix.h"
#include "matrix_ops.h"


// A must have more nnz than B
int csr_subtract(csr *A, csr *B, csr **result) {

  *result = malloc(sizeof(csr));
  
  (*result)->n = A->n;
  (*result)->m = A->m;

  // temp storage to figure out nnz of (*result), and get vals/cols
  int *temp_vals = (int *)malloc(A->nnz * sizeof(int));
  int *temp_cols = (int *)malloc(A->nnz * sizeof(int));
  int *temp_row_ptr = (int *)malloc((A->nnz + 1) * sizeof(int));
  temp_row_ptr[0] = 0;
  int temp_count = 0;


#pragma omp parallel
  {
#pragma omp for
    for (int row = 0; row < (*result)->n; row++) {
      int row_start_A = A->row_ptr[row];
      int row_end_A = A->row_ptr[row + 1];
      int row_start_B = B->row_ptr[row];
      int row_end_B = B->row_ptr[row + 1];

      int local_count = 0;

      for (int i = row_start_A; i < row_end_A; i++) {
	int col = A->cols[i];
	int A_val = A->vals[i];
	int B_val = 0;

	for (int j = row_start_B; j < row_end_B; j++) {
	  if (B->cols[j] == col) {
	    B_val = B->vals[j];
	    break;
	  }
	}

	int result_val = A_val - B_val;
	if (result_val != 0) {
	  temp_vals[temp_count + local_count] = result_val;
	  temp_cols[temp_count + local_count] = col;
	  local_count++;
	}
      }

      temp_row_ptr[row + 1] = temp_count + local_count;
      temp_count += local_count;
    }
  }


  printf("nnz in LHS: %d\n", temp_count);
  (*result)->nnz = temp_count;
  (*result)->vals = (double *)malloc((*result)->nnz * sizeof(double));
  (*result)->cols = (int *)malloc((*result)->nnz * sizeof(int));
  (*result)->row_ptr = (int *)malloc(((*result)->n + 1) * sizeof(int));

  for (int i = 0; i < (*result)->nnz; ++i) {
    (*result)->vals[i] = temp_vals[i];
    (*result)->cols[i] = temp_cols[i];
  }

  for (int i = 0; i <= (*result)->n; ++i) {
    (*result)->row_ptr[i] = temp_row_ptr[i];
  }

  // Free temporary arrays
  free(temp_vals);
  free(temp_cols);
  free(temp_row_ptr);

  return 0;
}

int spmv(csr *A, double *y, double *x)
{


}
