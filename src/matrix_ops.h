#include <stdbool.h>

int csr_add(csr *A, csr *B, double scalar, csr **result);
int spmv(csr *A, double *y, double *x);
bool is_symmetric_csr(csr *A);
