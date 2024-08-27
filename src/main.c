/*main.c -- 2D heat diffusion FEM solver
 *
 * Written on Friday, 26 July 2024.
  */

#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "mesh.h"
#include "assemble.h"
#include "matrix_ops.h"
#include "params.h"
#include "solver.h"

int main(int argc, char *argv[])
{

  if (argc != 2)
    {
      printf("Usage:\n");
      printf("\t ./diffusion <mesh_file>\n");
      return 0;
    }

  char *mesh_filename = argv[1];

  double *vx;          // length = num_nodes                -  x coordinates by nodes
  double *vy;          // length = num_nodes                -  y coordinates by nodes
  int *EToV;           // length = 3 * num_elements         -  vertices by element
  int *degree;         // length = num_nodes                -  degree - 1 (# neighboring elements) of node by node
  int *boundary_nodes; // length = # boundary nodes         -  indices of nodes on the boundary
  csr *adj_matrix;     // length = [num_nodes x num_nodes]  -  adjacency matrix of mesh
  csr *M;              // length = [num_nodes x num_nodes]  -  global mass matrix
  csr *S;              // length = [num_nodes x num_nodes]  -  global stiffness matrix
  csr *A;            // length = [num_nodes x num_nodes]  -  LHS with time discretization
  csr *RHS;            // length = [num_nodes x num_nodes]  - RHS with time discretization
  
  int num_nodes;
  int num_boundary_nodes;
  int num_elements;
  int num_edges;

  int t_step = .01;
  
  printf("Mesh file: %s\n", mesh_filename);
  int error = get_mesh(mesh_filename, &num_nodes, &num_elements, &vx, &vy, &EToV);
  
  error = get_mesh_degrees(EToV, &degree, &num_edges, num_elements, num_nodes);
  printf("Got Mesh\n\n");
  printf("# nodes: %d\n", num_nodes);
  printf("# elements: %d\n", num_elements);
  printf("# edges: %d\n", num_edges);

  error = tri_mesh_csr_adj_matrix(&adj_matrix, &boundary_nodes, &num_boundary_nodes, EToV, num_elements, num_edges, num_nodes);
  error = cuthill_mckee(adj_matrix, degree, &EToV, &vx, &vy, &boundary_nodes, num_boundary_nodes, num_nodes, num_elements);
  free(degree);
  printf("\nFound bandwidth reduced node ordering\n\n");
  free_csr_novals(adj_matrix);  
  
  assemble_matrices(&S, &M, 1.0, 1.0, EToV, vx, vy, boundary_nodes, num_nodes, num_boundary_nodes, num_elements);

  
  //getting lhs for crank-nicholson
  csr_add(M, S, -t_step/2, &A);
  csr_add(M, S, t_step/2, &RHS);

  printf("nnz in A: %d\n", A->nnz);
  printf("nnz in RHS: %d\n", A->nnz);

  if (is_symmetric_csr(A))
      printf("A is symmetric\n");
  else
    printf("A is not symmetric\n");
  if (is_symmetric_csr(M))
      printf("M is symmetric\n");
  else
    printf("M is not symmetric\n");
  if (is_symmetric_csr(S))
      printf("S is symmetric\n");
  else
    printf("S is not symmetric\n");
  if (is_symmetric_csr(RHS))
      printf("RHS is symmetric\n");
  else
    printf("RHS is not symmetric\n");

  
  
      

  write_csr(A, "plotting/A_row_ptr", "plotting/A_cols", "plotting/A_vals");
  
  free_csr(M);
  free_csr(S);

  // solution a previous timestep
  double *u_prev = (double*) malloc(sizeof(double) * num_nodes);
  double *u = (double*) malloc(sizeof(double) * num_nodes);

  // forcing vector
  double *f = (double*) calloc(sizeof(double), num_nodes);
  //double *f_prev = (double*) malloc(sizeof(double) * num_nodes);
  
  // set initial conditions
  for (int node = 0; node < num_nodes; node++)
    u_prev[node] = gaussian(vx[node], vy[node], .5, .5, .1, .1);

  linear_solver(u, u_prev, f, A, RHS, 400);
  
  
  free(f);
  free_csr(RHS);
  free_csr(A);
  free(u_prev);
  free(u);
  free(vx);
  free(vy);
  free(boundary_nodes);
  free(EToV);
  
  return error;
}
