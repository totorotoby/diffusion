/* stiffness_assemble.c -- assembling stiffness matrix using
 *  gaussian quad on reference element
 *
 * Written on Wednesday, 31 July 2024.
  */
  
#include <stdio.h>
#include <stdbool.h>
#include "assemble.h"
#include "params.h"
#include "mesh.h"
#include "matrix.h"
#include "linkedlist.h"

#define DIM 2
#define NUM_ABSCISSAS 4
#define NUM_BASIS 3

/*
  linear lagrangian basis horribly low order.
*/

double local_jacobian(double p1x, double p1y,
		   double p2x, double p2y,
		   double p3x, double p3y)
{

  return (((p2x-p1x) * (p3y-p1y)) -
	  ((p2y-p1y) * (p3x-p1x)))/4;
}

int assemble_matrices(csr **S, csr **M, const int *EToV, const double *vx, const double *vy,
		      const int num_nodes, const int num_elements)
{

  double quad_weights[NUM_ABSCISSAS] = {-9.0/8, 25.0/24, 25.0/24, 25.0/24};

  // quad eval points [ABS x DIM + 1]
  double abss[NUM_ABSCISSAS][DIM + 1] = {{-1.0/3, -1.0/3, 1.0},
					 {1.0/5, -3.0/5, 1.0},
					 {-3.0/5, 1.0/5, 1.0},
					 {-3.0/5, -3.0/5, 1.0}};

  // NOTE: Might need to reverse the order so that nodes and basis functions line up counter clockwise
  // linear basis function coefficents [BASIS x DIM + 1]
  double psi[NUM_BASIS][DIM + 1] = {{-1.0/2, -1.0/2, 0.0},
                                    {1.0/2, 0.0, 1.0/2},
                                    {0.0, 1.0/2, 1.0/2}};
  
  double psix[NUM_BASIS][DIM + 1] = {{0.0, -1.0/2, -1.0/2},
                                    {0.0, 0.0, 1.0/2},
                                    {0.0, 1.0/2, 0}};

  double psiy[NUM_BASIS][DIM + 1] = {{-1.0/2, 0.0, -1.0/2},
                                    {1.0/2, 0.0, 0.0},
                                    {0.0, 0.0, 1.0/2}};
  
  double local_mass[NUM_BASIS][NUM_BASIS];
  double local_stiffness[NUM_BASIS][NUM_BASIS];

  //-----------------------------------------------------------------------------
  // Assembling Global Matrices
  // creating linked list and row_ptr

  int *row_ptr = (int *) calloc((num_nodes + 1), sizeof(int));
  Node **column_lists = (Node **) malloc(num_nodes * sizeof(Node *));
  for (int row = 0; row < num_nodes; row++)
    {
      column_lists[row] = NULL;
    }

  for (int e = 0; e < num_elements; e++)
    {
      
      double p1x = vx[EToV[3*e]];
      double p1y = vy[EToV[3*e]];
      double p2x = vx[EToV[3*e + 1]];
      double p2y = vy[EToV[3*e + 1]];
      double p3x = vx[EToV[3*e + 2]];
      double p3y = vy[EToV[3*e + 2]];

      double J = local_jacobian(p1x, p1y,
				p2x, p2y,
				p3x, p3y);

      for (int v = 0; v < NUM_BASIS; v++) 
	{
	  int row = EToV[NUM_BASIS * e + v];
	  for (int u = 0; u < NUM_BASIS; u++)
	    {
	      int col = EToV[NUM_BASIS*e + u];

	      // generate local mass matrix
	      double mval = 0;
	      double sval = 0;
	      // loop over quadtrature
	      for (int k = 0; k < NUM_ABSCISSAS; k++)
		{
		  double fi = 0;
		  double fj = 0;
		  double fis = 0;
		  double fjs = 0;
		  double x = phys_coord(abss[k][0], abss[k][1],
					p1x, p2x, p3x);
		  double y = phys_coord(abss[k][0], abss[k][1],
					p1y, p2y, p3y);
		  double s = speed(x, y);
		  double t = diffusivity(x,y);
		  // loop over function eval
		  for (int l = 0; l < DIM + 1; l++)
		    {
		      fi += psi[v][l] * abss[k][l];
		      fj += psi[u][l] * abss[k][l];
		      fis += psix[v][l] * abss[k][l] + psiy[v][l] * abss[k][l];
		      fjs += psix[u][l] * abss[k][l] + psiy[u][l] * abss[k][l];
		    }
		  mval += fi * s * fj * quad_weights[k];
		  sval += fis * t * fjs * quad_weights[k];
		}
	      if (!is_col_in_list(column_lists[row], col)) {
		row_ptr[row + 1]++;
		add_to_list(&column_lists[row], col, J * mval, sval);
	      }
	      else {
		add_to_val(column_lists[row], col, J * mval, sval);
	      }
	    }
        }
    }

  // could prefix sum this
  for (int row = 0; row < num_nodes; row++)
    row_ptr[row + 1] += row_ptr[row];
  int nnz = row_ptr[num_nodes];
  printf("nnz in mass matrix: %d\n\n", nnz);
  
  // right now the mass and stiffness matrices have the same sparsity pattern
  allocate_csr_row_ptr(M, row_ptr, nnz, num_nodes, num_nodes);
  allocate_csr_row_ptr(S, row_ptr, nnz, num_nodes, num_nodes);

  int nz = 0;
  for (int row = 0; row < num_nodes; row++)
    {
      Node *current = column_lists[row];
      while (current != NULL)
	{
	  int col = current -> col;
	  (*M) -> cols[nz] = col;
	  (*S) -> cols[nz] = col;
	  (*M) -> vals[nz] = current->mval;
	  (*S) -> vals[nz] = current->sval;
	  
	  nz++;
	  current = current->next;
	}
      free_list(column_lists[row]);
    }

  write_csr(*M, "mass_row_ptr", "mass_cols", "mass_vals");
  
  free(column_lists);

  return 0;
}
