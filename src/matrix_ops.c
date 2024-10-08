/* matrix_ops.c -- matrix operations, solvers, LU, spmv, etc.
 *
 * Written on Thursday, 15 August 2024.
 */

#include "matrix.h"
#include "matrix_ops.h"
#include <stdio.h>

// A must have more nnz than B
int csr_add(csr *A, csr *B, double scalar, csr **result) {

  *result = malloc(sizeof(csr));
  (*result)->n = A->n;
  (*result)->m = A->m;

  // temp storage to figure out nnz of (*result), and get vals/cols
  double *temp_vals = (double *)malloc(A->nnz * sizeof(double));
  int *temp_cols = (int *)malloc(A->nnz * sizeof(int));
  int *temp_row_ptr = (int *)malloc((A->nnz + 1) * sizeof(int));
  temp_row_ptr[0] = 0;
  int temp_count = 0;
  int local_count = 0;

  // #pragma omp parallel for firstprivate(local_count)
  for (int row = 0; row < (*result)->n; row++) {

    int row_start_A = A->row_ptr[row];
    int row_end_A = A->row_ptr[row + 1];
    int row_start_B = B->row_ptr[row];
    int row_end_B = B->row_ptr[row + 1];

    local_count = 0;

    for (int i = row_start_A; i < row_end_A; i++) {
      int col = A->cols[i];
      double A_val = A->vals[i];
      double B_val = 0;

      for (int j = row_start_B; j < row_end_B; j++) {
        if (B->cols[j] == col) {
          B_val = B->vals[j];
          break;
        }
      }

      double result_val = A_val + scalar * B_val;
      if (result_val != 0.0) {
        temp_vals[temp_count + local_count] = result_val;
        temp_cols[temp_count + local_count] = col;
        local_count++;
      }
    }

    temp_row_ptr[row + 1] = temp_count + local_count;
#pragma omp atomic
    temp_count += local_count;
  }

  (*result)->nnz = temp_count;
  (*result)->vals = (double *)malloc((*result)->nnz * sizeof(double));
  (*result)->cols = (int *)malloc((*result)->nnz * sizeof(int));
  (*result)->row_ptr = (int *)malloc(((*result)->n + 1) * sizeof(int));

#pragma omp parallel for
  for (int i = 0; i < (*result)->nnz; ++i) {
    (*result)->vals[i] = temp_vals[i];
    (*result)->cols[i] = temp_cols[i];
  }

  for (int i = 0; i <= (*result)->n; ++i) {
    (*result)->row_ptr[i] = temp_row_ptr[i];
    // printf("%d\n", temp_row_ptr[i]);
  }

  // Free temporary arrays
  free(temp_vals);
  free(temp_cols);
  free(temp_row_ptr);

  return 0;
}

int spmv(csr *A, double *y, double *x) {

  int rows = A->n;
  double *vals = A->vals;
  int *cols = A->cols;
  int *row_ptr = A->row_ptr;

#pragma omp parallel for
  for (int row = 0; row < rows; row++) {

    int begin = row_ptr[row];
    int end = row_ptr[row + 1];
    double sum = 0;

    for (int i = begin; i < end; i++) {
      sum += vals[i] * x[cols[i]];
    }
    y[row] = sum;
  }
  return 0;
}

bool is_symmetric_csr(csr *A) {

  int n = A->n;
  int *row_ptr = A->row_ptr;
  int *cols = A->cols;
  double *vals = A->vals;

  bool sym = true;
  for (int row = 0; row < n; row++) {

    int begin = row_ptr[row];
    int end = row_ptr[row + 1];
    for (int i = begin; i < end; i++) {

      int col = cols[i];
      double val = vals[i];

      int sbegin = row_ptr[col];
      int send = row_ptr[col + 1];
      int found = 0;

      for (int j = sbegin; j < send; j++) {

        int scol = cols[j];
        double sval = vals[j];

        if (scol == row) {
          if (sval == val) {
            found = 1;
          } else {
            sym = false;
            //printf("miss match vals %f and %f at %d %d\n", sval, val, row, col);
          }
        }
      }
      if (found == 0) {
        sym = false;
        //printf("bad sparsity at %d %d\n", row, col);
      }
    }
  }
  return sym;
}


void csr_to_csc(int num_rows, int num_cols,
		const int *row_ptr, const int *cols, const double *vals,
                int **col_ptr, int **row_indices, double **values) {
 
    *col_ptr = (int *)calloc(num_cols + 1, sizeof(int));
 
     for (int i = 0; i < num_rows; i++) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            (*col_ptr)[cols[j] + 1]++;
        }
    }

     for (int i = 1; i <= num_cols; i++) {
        (*col_ptr)[i] += (*col_ptr)[i - 1];
    }

    int nnz = (*col_ptr)[num_cols];
    *row_indices = (int *)malloc(nnz * sizeof(int));
    *values = (double *)malloc(nnz * sizeof(double));
    int *next_index = (int *)malloc((num_cols + 1) * sizeof(int));
    
    for (int i = 0; i <= num_cols; i++) {
        next_index[i] = (*col_ptr)[i];
    }

    for (int i = 0; i < num_rows; i++) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            int col = cols[j];
            int dest = next_index[col];
            (*row_indices)[dest] = i;
            (*values)[dest] = vals[j];
            next_index[col]++;
        }
    }

    free(next_index);
}
