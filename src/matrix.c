/*matrix.c -- sparse matrix operations
 *
 * Written on Monday, 29 July 2024.
 */

#include "matrix.h"
#include <stdio.h>
#include <string.h>

int allocate_csr_novals(csr **matrix, int nnz, int n, int m) {

  *matrix = (csr *)malloc(sizeof(csr));
  (*matrix)->nnz = nnz;
  (*matrix)->n = n;
  (*matrix)->m = m;
  (*matrix)->cols = (int *)malloc(nnz * sizeof(int));
  (*matrix)->row_ptr = (int *)calloc((n + 1), sizeof(int));

  return 0;
}

int allocate_csr(csr **matrix, int nnz, int n, int m) {

  *matrix = (csr *)malloc(sizeof(csr));
  (*matrix)->nnz = nnz;
  (*matrix)->n = n;
  (*matrix)->m = m;
  (*matrix)->cols = (int *)malloc(nnz * sizeof(int));
  (*matrix)->vals = (double *)malloc(nnz * sizeof(double));
  (*matrix)->row_ptr = (int *)calloc((n + 1), sizeof(int));

  return 0;
}

int allocate_csr_row_ptr(csr **matrix, int *row_ptr, int nnz, int n, int m) {

  *matrix = (csr *)malloc(sizeof(csr));
  (*matrix)->nnz = nnz;
  (*matrix)->n = n;
  (*matrix)->m = m;
  (*matrix)->cols = (int *)malloc(nnz * sizeof(int));
  (*matrix)->vals = (double *)malloc(nnz * sizeof(double));
  (*matrix)->row_ptr = row_ptr;

  return 0;
}

int free_csr_novals(csr *matrix) {

  free(matrix->cols);
  free(matrix->row_ptr);
  free(matrix);

  return 0;
}

int free_csr(csr *matrix) {

  free(matrix->cols);
  free(matrix->vals);
  free(matrix->row_ptr);
  free(matrix);

  return 0;
}

int write_csr_novals(const csr *matrix, const char *f1, const char *f2) {

  int n = matrix->n;
  int nnz = matrix->nnz;
  int *cols = matrix->cols;
  int *row_ptr = matrix->row_ptr;

  FILE *fout = fopen(f1, "w");

  for (int row = 0; row < n + 1; row++) {
    fprintf(fout, "%d\n", row_ptr[row]);
  }
  fclose(fout);

  fout = fopen(f2, "w");
  for (int i = 0; i < nnz; i++) {
    fprintf(fout, "%d\n", cols[i]);
  }
  fclose(fout);

  return 0;
}

int write_csr(const csr *matrix, const char *f1, const char *f2,
              const char *f3) {

  int n = matrix->n;
  int nnz = matrix->nnz;
  int *cols = matrix->cols;
  int *row_ptr = matrix->row_ptr;
  double *vals = matrix->vals;

  FILE *fout = fopen(f1, "w");

  for (int row = 0; row < n + 1; row++) {
    fprintf(fout, "%d\n", row_ptr[row]);
  }
  fclose(fout);

  fout = fopen(f2, "w");
  for (int i = 0; i < nnz; i++) {
    fprintf(fout, "%d\n", cols[i]);
  }
  fclose(fout);

  fout = fopen(f3, "w");
  for (int i = 0; i < nnz; i++) {
    fprintf(fout, "%f\n", vals[i]);
  }
  fclose(fout);

  return 0;
}

int write_int_array(const int *array, const int n, const char *filename) {

  FILE *fout = fopen(filename, "w");

  for (int row = 0; row < n; row++) {
    fprintf(fout, "%d\n", array[row]);
  }
  fclose(fout);
}

int write_EToV_array(const int *array, const int n, const char *filename) {

  FILE *fout = fopen(filename, "w");

  for (int e = 0; e < n; e++) {
    fprintf(fout, "%d %d %d\n", array[e * 3], array[e * 3 + 1],
            array[e * 3 + 2]);
  }
  fclose(fout);
}

int write_double_array(const double *array, const int n, const char *filename) {

  FILE *fout = fopen(filename, "w");

  for (int row = 0; row < n; row++) {
    fprintf(fout, "%f\n", array[row]);
  }
  fclose(fout);
}
