/* params.c -- Physical Params for diffusion 
 *
 * Written on Sunday,  4 August 2024.
  */
  
  #include <stdio.h>
  #include "params.h"
  
double diffusivity(double x, double y)
{
  return x*x + y*y;
}

double phys_coord(double r, double s,
		  double p1, double p2, double p3)
{
  return ((r+s)/2) * p1 + ((r+1)/2) * p2 + ((s+1)/2) * p3;
}

double speed(double x, double y)
{
  return x + y;
}
