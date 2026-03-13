/**
 * @file genz.h
 * @brief Genz test integrand family for numerical integration testing.
 *
 * Declares the seven Genz test functions, their exact integrals, and
 * supporting utilities (error-function tables, tuple enumeration).
 *
 * @author  Knut Petras
 */

#ifndef GENZ_H
#define GENZ_H

#include "smolpack.h"

/* ------------------------------------------------------------------ */
/*                         Global parameters                          */
/* ------------------------------------------------------------------ */

/** @brief Coefficient vector (randomised per test run). */
extern double c[maxdim];

/** @brief Shift vector (randomised per test run). */
extern double w[maxdim];

/** @brief Precomputed constant sqrt(pi). */
extern double wp;

/** @brief Error-function Taylor coefficients for @c my_erfc. */
extern double derf[8][20];

/* ------------------------------------------------------------------ */
/*                        Genz integrands                             */
/* ------------------------------------------------------------------ */

/**
 * @brief Oscillatory test integrand: cos(2*pi*w[0] + sum(c[i]*x[i])).
 */
double f1(int dim, double x[]);

/**
 * @brief Product-peak test integrand: 1 / prod((1/c[i])^2+(x[i]-w[i])^2).
 */
double f2(int dim, double x[]);

/**
 * @brief Corner-peak test integrand: 1 / (1+sum(c[i]*x[i]))^(dim+1).
 */
double f3(int dim, double x[]);

/**
 * @brief Gaussian test integrand: exp(-sum((c[i]*(x[i]-w[i]))^2)).
 */
double f4(int dim, double x[]);

/**
 * @brief Continuous test integrand: exp(-sum(c[i]*|x[i]-w[i]|)).
 */
double f5(int dim, double x[]);

/**
 * @brief Discontinuous test integrand (zero outside a box).
 */
double f6(int dim, double x[]);

/**
 * @brief Exponential-sum test integrand: exp(sum(x[i])).
 */
double f7(int dim, double x[]);

/* ------------------------------------------------------------------ */
/*                       Support utilities                            */
/* ------------------------------------------------------------------ */

/**
 * @brief Return the exact integral of Genz function @p fnum.
 *
 * @param fnum  Function index (1..7).
 * @param dim   Spatial dimension.
 * @return      Exact integral over [0,1]^dim.
 */
double integral(int fnum, int dim);

/**
 * @brief Initialise Taylor-coefficient table for @c my_erfc.
 */
void init_erf(void);

/**
 * @brief Evaluate the complementary error function erfc(x).
 *
 * @param x  Argument.
 * @return   erfc(x).
 */
double my_erfc(double x);

/**
 * @brief Print the current date and time as a timestamp.
 */
void timestamp(void);

/**
 * @brief Compute the next element of a tuple space.
 *
 * @param m1    Minimum entry value.
 * @param m2    Maximum entry value.
 * @param n     Number of components.
 * @param rank  In/out rank counter (set to 0 on first call).
 * @param x     In/out tuple vector.
 */
void tuple_next(int m1, int m2, int n, int *rank, int x[]);

#endif /* GENZ_H */
