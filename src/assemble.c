/* stiffness_assemble.c -- assembling stiffness matrix using
 *  gaussian quad on reference element
 *
 * Written on Wednesday, 31 July 2024.
  */
  
#include <stdio.h>
#include "assemble.h"
#include "params.h"
#include "mesh.h"
#include "matrix.h"

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
				    
  double local_mass[NUM_BASIS][NUM_BASIS];
  double local_stiffness[NUM_BASIS][NUM_BASIS];

  // ---------------------------------------------------------------------------
  // Getting Local Matrices
  // not optimized at all...
  // loop over basis functions
  for (int i = 0; i < NUM_BASIS; i++)
    {
      // loop over basis functions
      for (int j = 0; j < NUM_BASIS; j++)
	{
	  double mval = 0;
	  // loop over quadtrature
	  for (int k = 0; k < NUM_ABSCISSAS; k++)
	    {
	      double fi = 0;
	      double fj = 0;
	      double s = speed(abss[k][0], abss[k][1]);
	      // loop over function eval
	      for (int l = 0; l < DIM + 1; l++)
		{
		  fi += psi[i][l] * abss[k][l];
		  fj += psi[j][l] * abss[k][l];
		}
	      mval += fi * s * fj * quad_weights[k];
	    }
	  local_mass[i][j] = mval;
	}
    }

  printf("unscaled local mass matrix: \n");
  for (int i = 0; i < NUM_BASIS; i++)
    {
      for (int j = 0; j < NUM_BASIS; j++)
	{
	  printf("%f ", local_mass[i][j]);
	}
      printf("\n");
    }

  //-----------------------------------------------------------------------------
  // Assembling Global Matrices
  // getting row_ptr
  int *row_ptr = (int *) malloc((num_nodes + 1) * sizeof(int));
  for (int e = 0; e < num_elements; e++)
    {
      for (int v = 0; v < NUM_BASIS; v++)
	{
	  int row = EToV[3*e + v];
	  row_ptr[row]++;
	}
    }
  // could prefix sum this
  for (int row = 0; row < num_nodes; row++)
    row_ptr[row + 1] += row_ptr[row];
  
  int nnz = row_ptr[num_nodes];
  (*M) = (csr *) malloc(sizeof(csr));
  (*M) -> row_ptr = row_ptr;
  (*M) -> cols = (int *) malloc(nnz * sizeof(int));
  (*M) -> vals = (double *) malloc(nnz * sizeof(double));

  // loop through elements scale local matrices by jacobian and place in correct global system

  int *index = (int *) malloc(num_nodes * sizeof(int));
  for (int row = 0; row < num_nodes; row++)
      index[row] = row_ptr[row];

  for (int e = 0 ; e < num_elements ; e++)
    {

      double J = local_jacobian(vx[EToV[3*e]], vy[EToV[3*e]],
				vx[EToV[3*e + 1]], vy[EToV[3*e + 1]],
				vx[EToV[3*e + 2]], vy[EToV[3*e + 2]]);
      
      for (int v = 0; v < NUM_BASIS; v++)
	{
	  int row = EToV[3*e + v];
	  for (int u = 0; u < NUM_BASIS; u++)
	    {

	      int col = EToV[3*e + u];
	      int found = 0;

	      for (int pos = row_ptr[row]; pos < index[row]; pos++)
		{
		  if ((*M) -> cols[pos] == col)
		    {
		      (*M) -> vals[pos] += J * local_mass[u][v];
		      found = 1;
		      break;
		    }
		}

	      if (!found)
		{
		  fprintf(stderr, "index[%d]: %d", row, index[row]);
		  (*M) -> cols[index[row]] = col;
		  (*M) -> vals[index[row]] = J * local_mass[u][v];
		  index[row]++;
		}
	    }
	}
    }

  free(index);

  return 0;
}
