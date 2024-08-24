/* params.c -- Physical Params for diffusion 
 *
 * Written on Sunday,  4 August 2024.
  */
  
#include <stdio.h>
#include <math.h>
#include "params.h"

double diffusivity(double x, double y)
{
  return 1.0;
  //return x*x + y*y;
}

double phys_coord(double r, double s,
		  double p1, double p2, double p3)
{
  return ((r+s)/2) * p1 + ((r+1)/2) * p2 + ((s+1)/2) * p3;
}

double specific_yield(double x, double y)
{
  return 1.0;
  //return x + y;
}

double gaussian(double x, double y,
		double mu_x, double mu_y,
		double sigma_x, double sigma_y) {
  
  double coef = 1.0 / (2.0 * M_PI * sigma_x * sigma_y);
  double exp_term = exp(-0.5 * ((pow((x - mu_x) / sigma_x, 2.0)) + (pow((y - mu_y) / sigma_y, 2.0))));
  return coef * exp_term;
}

