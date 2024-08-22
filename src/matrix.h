#ifndef MATRIX_HEADER
#define MATRIX_HEADER

#include <stdlib.h>

typedef struct csr {

  int n;
  int m;
  int nnz;
  double *vals;
  int *cols;
  int *row_ptr;
  
} csr;

int allocate_csr_novals(csr **matrix, const int nnz, const int n, int m);
int allocate_csr(csr **matrix, const int nnz, const int n, int m);
int allocate_csr_row_ptr(csr **matrix, int *row_ptr, const int nnz, const int n, int m);
int free_csr_novals(csr *matrix);
int free_csr(csr *matrix);
int write_csr_novals(const csr *matrix, const char *f1, const char *f2);
int write_csr(const csr *matrix, const char *f1, const char *f2, const char *f3);
int write_int_array(const int *array, const int n, const char *filename);
int write_double_array(const double *array, const int n, const char *filename);
int write_EToV_array(const int *array, const int n, const char *filename);

#endif
