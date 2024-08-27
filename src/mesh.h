#include "matrix.h"

int min(const int a, const int b);
int max(const int a, const int b);
int get_mesh(const char *filename, int *num_nodes, int *num_elements,
             double **vx, double **vy, int **EToV);
int get_boundary_nodes(const int *EToV, const int *neighbors,
                       int **boundary_nodes, int num_nodes, int num_elements);
int get_mesh_degrees(const int *EToV, int **degree, int *num_edges,
                     const int num_elements, const int num_nodes);
int tri_mesh_csr_adj_matrix(csr **adj_matrix, int **boundary_nodes,
                            int *num_boundary_nodes, const int *EToV,
                            const int num_elements, const int num_edges,
                            int num_nodes);
int degree_comp(const void *v1, const void *v2, const void *degree);
int cuthill_mckee(const csr *adj_matrix, const int *degree, int **EToV,
                  double **vx, double **vy, int **boundary_nodes,
                  const int num_boundary_nodes, const int num_nodes,
                  const int num_elements);
