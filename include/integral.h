#ifndef __INTEGRAL_H__
#define __INTEGRAL_H__

#include "basis.h"
#include "gslextra.h"

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>

typedef struct gaussian_chain{
    double R[3];
    int a[3];
    
    double exponent;

    double coefficient;

    gaussian_chain * NEXT;

}gaussian_chain;

gaussian_chain * gaussian_chain_calloc();

void gaussian_chain_free(gaussian_chain * HEAD);

double gaussian_chain_get(gaussian_chain * HEAD,double x, double y, double z);
double orbital_get(orbital * orbital,double x, double y, double z);

double SIntegral(double ra[3], double rb[3], int ax, int ay, int az, int bx, int by, int bz, double alpha,double beta);

double gaussian_chain_SIntegral(gaussian_chain * a, gaussian_chain * b);
double gaussian_chain_full_SIntegral(gaussian_chain * a_HEAD, gaussian_chain * b_HEAD);
void single_electron_transform(gaussian_chain * HEAD, orbital * a);

double orbital_SIntegral(orbital * a,orbital * b);

void orbital_S_matrix(gsl_matrix * dest, orbital * HEAD, int length);

#endif