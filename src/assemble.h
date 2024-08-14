int assemble_matices(const int *EToV, const double *vx, const double *vy,
		     const int num_nodes, const int num_elements);
int local_jacobian(double *J, double *p1, double *p2, double *p3);
int psi0(double *psi, double r, double s);
int psi1(double *psi, double r, double s);
int psi2(double *psi, double r, double s);
