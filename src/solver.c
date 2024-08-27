/*solver.c --
 *
 * Written on Monday, 26 August 2024.
 */

#include "solver.h"
#include "matrix.h"
#include "matrix_ops.h"
#include <stdio.h>

int linear_solver(double *u, double *u_prev, double *f, csr *A, csr *RHS,
                  int tsteps) {

  int num_nodes = A->n;
  double *b = (double *)malloc(num_nodes * sizeof(double));

  for (int t = 0; t < tsteps; t++) {

    spmv(RHS, b, u_prev);
  }

  free(b);
}
