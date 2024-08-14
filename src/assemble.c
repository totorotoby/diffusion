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

#define NV_PER_ELEM 3
#define DIM 2
#define NUM_ABSCISSAS 4
#define NUM_BASIS 3

/*
  linear lagrangian basis horribly low order.
*/

double psi0(double r, double s)
{
  return .5*s + .5;
}

double psi0r(double r, double s)
{
  return 0.0;
}

double psi0s(double r, double s)
{
  return 0.5;
}

double psi1(double r, double s)
{
  return .5*r + .5;
}

double psi1r(double r, double s)
{
  return .5;
}

double psi1(double r, double s)
{
  return 0;
}

double psi2(double r, double s)
{
  return -.5*r - .5*s;
}

double psi2r(double r, double s)
{
  return -.5;
}

double psi2s(double r, double s)
{
  return -.5;
}

int local_jacobian(double *J,
		   double *p1x, double *p1y,
		   double *p2x, double *p2y,
		   double *p3x, double *p3y)
{
  *J = (((p2x-p1x) * (p3y-p1y)) -
	((p2y-p1y) * (p3x-p1x)))/4;
  return 0;
}

int assemble_matrices(csr **S, csr **M, const int *EToV, const double *vx, const double *vy,
		      const int num_nodes, const int num_elements)
{

  double quad_weights[NUM_ABSCISSAS] = {9/32, 25/96, 25/96, 25/96};
  double absx[NUM_ABSCISSAS] = {1/3, 3/5, 1/5, 1/5};
  double absy[NUM_ABSCISSAS] = {1/3, 1/5, 3/5, 1/5};
  double (*basis[NUM_BASIS]) (double r, double s) = {psi0, psi1, psi2};
  double (*basisr[NUM_BASIS]) (double r, double s) = {psi0r, psi1r, psi2r};
  double (*basiss[NUM_BASIS]) (double r, double s) = {psi0s, psi1s, psi2s};
									 
  double local_mass[NUM_BASIS][NUM_BASIS];
  double local_stiffness[NUM_BASIS][NUM_BASIS];

  // construct unscaled local mass, and stiffness matrix
  // probably more efficent way to do this with MM and unrolling or something
  for (int i = 0; i < NUM_BASIS ; i++)
    {
      for (int j = 0;  < NUM_BASIS ; j++)
      {
	double mval = 0;
	double rval = 0;
	double sval = 0;
	for (int k = 0 ; k < NUM_ABSCISSAS ; k++)
	  {
	    double ax = absx[k];
	    double ay = absy[k];
	    val += quad_weights[k] * basis[i](ax, ay) * basis[j](ax, ay) * speed(ax, ay);
	    rval += quad_weights[k] * basisr[i](ax, ay) * basisr[j](ax, ay) * diffusivity(ax, ay);
	    sval += quad_weights[k] * basiss[i](ax, ay) * basiss[j](ax, ay) * diffusivity(ax, ay);
	  }
	local_mass[i][j] = val;
	local_stiffness[i][j] = rval + sval;
      }
    }

  // loop through elements scale local matrices by jacobian and place in correct global system
  for (int e = 0 ; e < num_elements ; e++)
    {

      double J;
      local_jacobian(&J,
		     vx[EToV[e]], vy[EToV[e]],
		     vx[EToV[e + 1]], vy[EToV[e + 1]],
		     vx[EToV[e + 2]], vy[EToV[e + 2]]);
      
    }
}
