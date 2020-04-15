#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <math.h>
# include "myfunction.h"
//static const double t1 = 0.15;
//static const double t2 = 2.18;
//static const double t3 = 0.725;
//static const double t4 = 0.45;

/* Exponential rejection sampling (a,inf) */
double ers_a_inf(double a,const gsl_rng * r) {
  //SAMPLER_DEBUG("ers_a_inf", a, R_PosInf);
  const double ainv = 1.0 / a;
  double x, rho;
  do {
    //x = rexp(ainv) + a; /* rexp works with 1/lambda */
    x=gsl_ran_exponential (r, ainv)+a;
  rho = exp(-0.5 * pow((x - a), 2));
  } while (gsl_ran_flat (r, 0,1) > rho);
  return x;
}

/* Normal rejection sampling (a,inf) 
double nrs_a_inf(double a,const gsl_rng * r) {
  //SAMPLER_DEBUG("nrs_a_inf", a, R_PosInf);
  double x = -DBL_MAX;
  while (x < a) {
    x = gsl_ran_ugaussian(r);
  }
  return x;
}
*/

/*
 double r_lefttruncnorm(double a, double mean, double sd,const gsl_rng * r) {
  const double alpha = (a - mean) / sd;
  if (alpha < t4) {
    return mean + sd * nrs_a_inf(alpha,r);
  } else {
    return mean + sd * ers_a_inf(alpha,r);
  }
}
double r_righttruncnorm(double b, double mean, double sd,const gsl_rng * r) {
  const double beta = (b - mean) / sd;
  // Exploit symmetry: 
  return mean - sd * r_lefttruncnorm(-beta, 0.0, 1.0,r);
}
*/






  
  

