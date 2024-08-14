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

double speed(double x, double y)
{
  return x + y;
}
