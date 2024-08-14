#include "matrix.h"

int assemble_matrices(csr **S, csr **M, const int *EToV, const double *vx, const double *vy,
		      const int num_nodes, const int num_elements);
double local_jacobian(double p1x, double p1y,
		   double p2x, double p2y,
		   double p3x, double p3y);
int psi0(double *psi, double r, double s);
int psi1(double *psi, double r, double s);
int psi2(double *psi, double r, double s);
