/**
 * @file smolpack.h
 * @brief SMOLPACK — Smolyak Cubature for Multidimensional Integration
 *
 * SMOLPACK provides efficient numerical integration (cubature) over the
 * unit hypercube [0,1]^d using Smolyak's algorithm with Clenshaw-Curtis
 * quadrature rules. Two variants are available:
 *
 *   - **Delayed Clenshaw-Curtis** (`int_smolyak`): Uses a specially
 *     constructed basic sequence that produces fewer function evaluations
 *     for a given accuracy level.
 *
 *   - **Standard Clenshaw-Curtis** (`cc_int_smolyak`): Uses the classical
 *     Clenshaw-Curtis nested sequence as the underlying 1D rule.
 *
 * Both methods construct sparse grids via the Smolyak combination technique,
 * achieving polynomial exactness with substantially fewer nodes than a full
 * tensor-product grid, especially in high dimensions.
 *
 * @par Usage Example
 * @code
 *     #include "smolpack.h"
 *
 *     double my_func(int dim, double x[]) {
 *         double sum = 0.0;
 *         for (int i = 0; i < dim; i++) sum += x[i];
 *         return exp(sum);
 *     }
 *
 *     int main(void) {
 *         int dim = 3, q = dim + 2, print_stats = 1;
 *         double result = int_smolyak(dim, q, my_func, print_stats);
 *         printf("Integral = %e\n", result);
 *         return 0;
 *     }
 * @endcode
 *
 * @author  Knut Petras
 * @date    2007
 *
 * @par References
 *   -# Alan Genz,
 *      "A Package for Testing Multiple Integration Subroutines",
 *      in <em>Numerical Integration: Recent Developments, Software
 *      and Applications</em>, edited by Patrick Keast and Graeme
 *      Fairweather, Reidel, 1987, pp. 337–340,
 *      ISBN 9027725144.
 *   -# Erich Novak, Klaus Ritter,
 *      "High Dimensional Integration of Smooth Functions over Cubes",
 *      <em>Numerische Mathematik</em>, Vol. 75, No. 1,
 *      November 1996, pp. 79–97.
 *   -# Erich Novak, Klaus Ritter,
 *      "Simple Cubature Formulas with High Polynomial Exactness",
 *      <em>Constructive Approximation</em>, Vol. 15, No. 4,
 *      December 1999, pp. 499–522.
 *   -# Knut Petras,
 *      "Fast Calculation of Coefficients in the Smolyak Algorithm",
 *      <em>Numerical Algorithms</em>, Vol. 26, No. 2,
 *      February 2001, pp. 93–109.
 *   -# Knut Petras,
 *      "Smolyak Cubature of Given Polynomial Degree with Few Nodes
 *       for Increasing Dimension",
 *      <em>Numerische Mathematik</em>, Vol. 93, No. 4,
 *      February 2003, pp. 729–753.
 *   -# Sergey Smolyak,
 *      "Quadrature and Interpolation Formulas for Tensor Products
 *       of Certain Classes of Functions",
 *      <em>Doklady Akademii Nauk SSSR</em>, Vol. 4, 1963,
 *      pp. 240–243.
 *
 * @note This implementation is C11 compliant.
 */

#ifndef SMOLPACK_H
#define SMOLPACK_H

/**
 * @brief Maximum spatial dimension supported by SMOLPACK.
 *
 * All arrays indexed by spatial dimension are sized to accommodate
 * up to @c maxdim components. Attempting to integrate in more than
 * @c maxdim dimensions results in undefined behavior.
 */
#define maxdim 40

/**
 * @brief Global counter for function evaluations.
 *
 * This variable is reset to zero at the start of each integration call
 * (`int_smolyak` or `cc_int_smolyak`) and is incremented by the user-
 * supplied integrand function on each evaluation. After the integration
 * completes, @c count holds the total number of integrand evaluations
 * performed.
 *
 * @note The user's integrand callback is expected to execute @c count++
 *       on each invocation for accurate function-call statistics.
 */
extern int count;

/**
 * @brief Approximate a multidimensional integral using the delayed
 *        Clenshaw-Curtis Smolyak algorithm.
 *
 * Constructs a sparse grid from a "delayed" Clenshaw-Curtis basic
 * sequence and applies the Smolyak combination technique to approximate
 *
 * @f[
 *   \int_{[0,1]^d} f(\mathbf{x})\,d\mathbf{x}
 * @f]
 *
 * The delayed sequence typically requires fewer function evaluations
 * than the standard Clenshaw-Curtis sequence for the same level of
 * polynomial exactness.
 *
 * @param[in] dim          Spatial dimension of the integrand.
 *                         Must satisfy @c 1 <= dim < maxdim.
 * @param[in] qq           Quadrature level parameter. The "number of
 *                         stages" is @c k = qq - dim. Higher values of
 *                         @c k yield more accurate approximations.
 *                         Must satisfy @c qq - dim < 48.
 * @param[in] ff           Pointer to the integrand function. The function
 *                         receives the dimension @c d and an array
 *                         @c x[0..d-1] with each component in [0,1].
 * @param[in] print_stats  If nonzero, print function-call and weight-
 *                         evaluation counts to stdout after integration.
 *
 * @return  The approximate value of the integral.
 *
 * @see cc_int_smolyak
 */
extern double int_smolyak(int dim, int qq,
                          double (*ff)(int, double x[]),
                          int print_stats);

/**
 * @brief Approximate a multidimensional integral using the standard
 *        Clenshaw-Curtis Smolyak algorithm.
 *
 * Constructs a sparse grid from the standard Clenshaw-Curtis nested
 * sequence and applies the Smolyak combination technique to approximate
 *
 * @f[
 *   \int_{[0,1]^d} f(\mathbf{x})\,d\mathbf{x}
 * @f]
 *
 * This variant uses the classical Clenshaw-Curtis rule hierarchy
 * (1, 3, 5, 9, 17, 33, 65, … nodes) and may require more function
 * evaluations than `int_smolyak` at the same level but can be more
 * accurate for smooth integrands.
 *
 * @param[in] dim          Spatial dimension of the integrand.
 *                         Must satisfy @c 1 <= dim < maxdim.
 * @param[in] qq           Quadrature level parameter. The "number of
 *                         stages" is @c k = qq - dim. Higher values of
 *                         @c k yield more accurate approximations.
 *                         Must satisfy @c qq - dim < 48.
 * @param[in] ff           Pointer to the integrand function. The function
 *                         receives the dimension @c d and an array
 *                         @c x[0..d-1] with each component in [0,1].
 * @param[in] print_stats  If nonzero, print function-call and weight-
 *                         evaluation counts to stdout after integration.
 *
 * @return  The approximate value of the integral.
 *
 * @see int_smolyak
 */
extern double cc_int_smolyak(int dim, int qq,
                             double (*ff)(int, double x[]),
                             int print_stats);

#endif /* SMOLPACK_H */
