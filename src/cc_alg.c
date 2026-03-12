/**
 * @file cc_alg.c
 * @brief Public API and core algorithm for standard CC Smolyak cubature.
 *
 * Contains the public @c cc_int_smolyak entry point plus the internal
 * @c cc_formula, @c cc_eval, and @c cc_fsum recursions.
 *
 * @author  Knut Petras
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "smolpack.h"
#include "cc_internal.h"

/******************************************************************************/
/**
 * @brief Approximate a multidimensional integral using standard CC Smolyak.
 * @date  25 April 2007
 *
 * @param[in] dim          Spatial dimension (1 <= dim < maxdim).
 * @param[in] qq           Cubature parameter (qq - dim < CC_GESFN).
 * @param[in] ff           Integrand function pointer.
 * @param[in] print_stats  Non-zero to print function-call / weight statistics.
 * @return                 Approximated value of the integral.
 */
double cc_int_smolyak(int dim, int qq, double (*ff)(int, double xx[]),
                      int print_stats) {
    /*
     * Make the parameters global.
     */
    cc_d = dim;
    cc_q = qq;
    cc_f = ff;
    /*
     * Initialisation.
     */
    cc_wcount = 0;
    count = 0;
    cc_quafo = 0.0;

    cc_setup(dim);
    /*
     * Call the Smolyak algorithm.  Q-DIM = K, the number of stages.
     */
    cc_formula(1, qq - dim);
    /*
     * Free the allocated memory.
     */
    cc_frei(cc_root);
    /*
     * Print statistics, if desired.
     */
    if (print_stats) {
        printf("\n");
        printf("  Number of function calls =  %d\n", count);
        printf("  Weight evaluations =        %d\n", cc_wcount);
    }

    return cc_quafo;
}

/******************************************************************************/
/**
 * @brief Carry out the Smolyak combination algorithm.
 * @date  25 April 2007
 *
 * @param k  Current dimension index.
 * @param l  Remaining index sum to distribute across dimensions.
 */
void cc_formula(int k, int l) {
    int i;

    if (k == cc_d + 1) {
        cc_quafo = cc_quafo + cc_eval(0);
    } else {
        for (i = 0; i <= l; i++) {
            if (cc_sw[i] < CC_FN) {
                cc_indices[k] = cc_sw[i];
                cc_formula(k + 1, l - i);
            }
        }
    }
    return;
}

/******************************************************************************/
/**
 * @brief Evaluate a product formula using tree-cached coefficients.
 * @date  25 April 2007
 *
 * @param k  Recursion depth (0 = entry, d+1 = leaf evaluation).
 * @return   Accumulated product-formula value.
 */
double cc_eval(int k) {
    int i;

    if (k == 0) {
        cc_summe = 0.0;
        cc_eval(1);
    } else if (k == cc_d + 1) {
        /*
         * Summation corresponding to one coefficient.
         */
        cc_summe = cc_summe + cc_coeff() * cc_fsum(0);
    } else {
        /*
         * Choice of the nodes.
         */
        for (i = 0; i <= cc_n[cc_indices[k]]; i++) {
            cc_argind[k] = i;
            cc_eval(k + 1);
        }
    }
    return cc_summe;
}

/******************************************************************************/
/**
 * @brief Compute unweighted sum of integrand values at symmetric nodes.
 * @date  26 April 2007
 *
 * Evaluates the integrand at nodes where the same weight applies,
 * exploiting the symmetry f(x) + f(1-x).
 *
 * @param k  Recursion depth (0 = entry, d+1 = leaf evaluation).
 * @return   Accumulated function sum.
 */
double cc_fsum(int k) {
    if (k == 0) {
        cc_ftotal = 0.0;
        cc_fsum(1);
    } else if (k == cc_d + 1) {
        cc_ftotal = cc_ftotal + (*cc_f)(cc_d, cc_x);
    } else {
        if (cc_indices[k] == 0) {
            cc_x[k - 1] = 0.5;
            cc_fsum(k + 1);
        } else {
            /*
             * Use symmetry to get both X and -X.
             */
            cc_x[k - 1] = cc_xnu[cc_indices[k]][2 * cc_argind[k] + 1];
            cc_fsum(k + 1);
            cc_x[k - 1] = 1.0 - cc_x[k - 1];
            cc_fsum(k + 1);
        }
    }
    return cc_ftotal;
}
