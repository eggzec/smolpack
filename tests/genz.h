/**
 * @file genz.h
 * @brief Genz test integrand family for numerical integration testing.
 *
 * Declares the seven Genz test functions, their exact integrals, and
 * supporting utilities (error-function tables, tuple enumeration).
 *
 * The test suite is drawn from:
 *
 *   Alan Genz, "A Package for Testing Multiple Integration
 *   Subroutines", in <em>Numerical Integration: Recent Developments,
 *   Software and Applications</em>, edited by Patrick Keast and
 *   Graeme Fairweather, Reidel, 1987, pp. 337–340.
 *
 * In all descriptions below, @e d denotes the spatial dimension,
 * @b c is a randomly chosen scaling vector, @b w is a randomly
 * chosen displacement/shift vector, and the domain of integration
 * is @f$[0,1]^d@f$.
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
 * @brief Genz #1 — Oscillatory.
 *
 * @f[
 *   f_1(\mathbf{x}) = \cos\!\bigl(2\pi\,w_0 +
 *   \textstyle\sum_{i=0}^{d-1} c_i\,x_i\bigr)
 * @f]
 */
double f1(int dim, double x[]);

/**
 * @brief Genz #2 — Product Peak.
 *
 * @f[
 *   f_2(\mathbf{x}) = \prod_{i=0}^{d-1}
 *   \frac{1}{c_i^{-2} + (x_i - w_i)^2}
 * @f]
 */
double f2(int dim, double x[]);

/**
 * @brief Genz #3 — Corner Peak.
 *
 * @f[
 *   f_3(\mathbf{x}) = \frac{1}{\bigl(1 +
 *   \textstyle\sum_{i=0}^{d-1} c_i\,x_i\bigr)^{d+1}}
 * @f]
 */
double f3(int dim, double x[]);

/**
 * @brief Genz #4 — Gaussian.
 *
 * @f[
 *   f_4(\mathbf{x}) = \exp\!\Bigl(-\textstyle\sum_{i=0}^{d-1}
 *   c_i^2\,(x_i - w_i)^2\Bigr)
 * @f]
 */
double f4(int dim, double x[]);

/**
 * @brief Genz #5 — Continuous.
 *
 * @f[
 *   f_5(\mathbf{x}) = \exp\!\Bigl(-\textstyle\sum_{i=0}^{d-1}
 *   c_i\,|x_i - w_i|\Bigr)
 * @f]
 */
double f5(int dim, double x[]);

/**
 * @brief Genz #6 — Discontinuous.
 *
 * @f[
 *   f_6(\mathbf{x}) =
 *   \begin{cases}
 *     \exp\!\bigl(\textstyle\sum_{i=0}^{d-1} c_i\,x_i\bigr)
 *       & \text{if } x_0 \le w_0 \text{ and } x_1 \le w_1,\\
 *     0 & \text{otherwise.}
 *   \end{cases}
 * @f]
 */
double f6(int dim, double x[]);

/**
 * @brief Genz #7 — Exponential Sum.
 *
 * @f[
 *   f_7(\mathbf{x}) = \exp\!\Bigl(
 *   \textstyle\sum_{i=0}^{d-1} x_i\Bigr)
 * @f]
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
