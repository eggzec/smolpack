/**
 * @file smolyak_alg.c
 * @brief Smolyak combination algorithm and public API for the delayed
 *        Clenshaw-Curtis subsystem.
 *
 * Contains the public entry point @c int_smolyak() together with the
 * recursive formula dispatcher, product-formula evaluator, and the
 * symmetric function-value summation.
 */

#include <stdio.h>
#include <stdlib.h>

#include "smolyak_internal.h"

/* ------------------------------------------------------------------ */
/*                   Public API implementation                        */
/* ------------------------------------------------------------------ */

/**
 * @brief Approximate a multidimensional integral (delayed Clenshaw-Curtis).
<<<<<<< HEAD
<<<<<<< HEAD
 * @date  28 April 2007
=======
>>>>>>> f26fbd8 (Implement Smolyak quadrature subsystem with Clenshaw-Curtis integration)
=======
 * @date  28 April 2007
>>>>>>> 62f2bf2 (Refactor Smolyak algorithms and add Genz test suite)
 *
 * Constructs a sparse grid from a "delayed" Clenshaw-Curtis basic
 * sequence and applies the Smolyak combination technique.
 *
 * @param[in] dim          Spatial dimension (1 <= dim < maxdim).
 * @param[in] qq           Level parameter; k = qq - dim stages.
 * @param[in] ff           Integrand callback.
 * @param[in] print_stats  Nonzero to print evaluation statistics.
 * @return Approximated integral value.
 */
double int_smolyak(int dim, int qq, double (*ff)(int, double xx[]),
                   int print_stats) {
    sm_d = dim;
    sm_q = qq;
    sm_f = ff;

    sm_wcount = 0;
    count = 0;
    sm_quafo = 0.0;

    sm_init(dim);

    sm_formula(1, qq - dim);

    sm_frei(sm_root);

    if (print_stats) {
        printf("\n");
        printf("  Number of function calls =  %d\n", count);
        printf("  Weight evaluations =        %d\n", sm_wcount);
    }

    return sm_quafo;
}

/* ------------------------------------------------------------------ */
/*                     Internal algorithms                            */
/* ------------------------------------------------------------------ */

/**
 * @brief Carry out the Smolyak combination algorithm.
<<<<<<< HEAD
<<<<<<< HEAD
 * @date  25 April 2007
=======
>>>>>>> f26fbd8 (Implement Smolyak quadrature subsystem with Clenshaw-Curtis integration)
=======
 * @date  25 April 2007
>>>>>>> 62f2bf2 (Refactor Smolyak algorithms and add Genz test suite)
 *
 * Recursively distributes the level-sum budget across dimensions.
 * When all dimensions have been assigned a formula index the
 * product-formula evaluation @c sm_eval() is invoked.
 *
 * @param[in] k  Current dimension index (1-based).
 * @param[in] l  Remaining index-sum budget.
 */
void sm_formula(int k, int l) {
    int i;

    if (k == sm_d + 1) {
        sm_quafo = sm_quafo + sm_eval(0);
    } else {
        for (i = 0; i <= l; i++) {
            if (sm_sw[i] < SM_FN) {
                sm_indices[k] = sm_sw[i];
                sm_formula(k + 1, l - i);
            }
        }
    }
}

/**
 * @brief Calculate the value of a product formula.
<<<<<<< HEAD
<<<<<<< HEAD
 * @date  25 April 2007
=======
>>>>>>> f26fbd8 (Implement Smolyak quadrature subsystem with Clenshaw-Curtis integration)
=======
 * @date  25 April 2007
>>>>>>> 62f2bf2 (Refactor Smolyak algorithms and add Genz test suite)
 *
 * Iterates over all node combinations for the current formula
 * assignment, multiplying the weight coefficient by the symmetrised
 * function sum.
 *
 * @param[in] k  Recursion depth (0 = entry, d+1 = leaf).
 * @return Accumulated product-formula value.
 */
double sm_eval(int k) {
<<<<<<< HEAD
<<<<<<< HEAD
=======
    double dummy;
>>>>>>> f26fbd8 (Implement Smolyak quadrature subsystem with Clenshaw-Curtis integration)
=======
>>>>>>> 62f2bf2 (Refactor Smolyak algorithms and add Genz test suite)
    int i;

    if (k == 0) {
        sm_summe = 0.0;
<<<<<<< HEAD
<<<<<<< HEAD
        sm_eval(1);
=======
        dummy = sm_eval(1);
>>>>>>> f26fbd8 (Implement Smolyak quadrature subsystem with Clenshaw-Curtis integration)
=======
        sm_eval(1);
>>>>>>> 62f2bf2 (Refactor Smolyak algorithms and add Genz test suite)
    } else if (k == sm_d + 1) {
        sm_summe = sm_summe + sm_coeff() * sm_fsum(0);
    } else {
        for (i = 0; i <= sm_n[sm_indices[k]]; i++) {
            sm_argind[k] = i;
<<<<<<< HEAD
<<<<<<< HEAD
            sm_eval(k + 1);
=======
            dummy = sm_eval(k + 1);
>>>>>>> f26fbd8 (Implement Smolyak quadrature subsystem with Clenshaw-Curtis integration)
=======
            sm_eval(k + 1);
>>>>>>> 62f2bf2 (Refactor Smolyak algorithms and add Genz test suite)
        }
    }
    return sm_summe;
}

/**
 * @brief Compute symmetric sums of integrand values.
<<<<<<< HEAD
<<<<<<< HEAD
 * @date  26 April 2007
=======
>>>>>>> f26fbd8 (Implement Smolyak quadrature subsystem with Clenshaw-Curtis integration)
=======
 * @date  26 April 2007
>>>>>>> 62f2bf2 (Refactor Smolyak algorithms and add Genz test suite)
 *
 * Evaluates the integrand at nodes for which the same weight
 * applies, exploiting symmetry: f(x) and f(1-x) share a weight.
 *
 * @param[in] k  Recursion depth (0 = entry, d+1 = leaf).
 * @return Accumulated function sum.
 */
double sm_fsum(int k) {
<<<<<<< HEAD
<<<<<<< HEAD

    if (k == 0) {
        sm_ftotal = 0.0;
        sm_fsum(1);
=======
    double dummy;

    if (k == 0) {
        sm_ftotal = 0.0;
        dummy = sm_fsum(1);
>>>>>>> f26fbd8 (Implement Smolyak quadrature subsystem with Clenshaw-Curtis integration)
=======

    if (k == 0) {
        sm_ftotal = 0.0;
        sm_fsum(1);
>>>>>>> 62f2bf2 (Refactor Smolyak algorithms and add Genz test suite)
    } else if (k == sm_d + 1) {
        sm_ftotal = sm_ftotal + (*sm_f)(sm_d, sm_x);
    } else {
        if (sm_indices[k] == 0) {
            sm_x[k - 1] = 0.5;
<<<<<<< HEAD
<<<<<<< HEAD
            sm_fsum(k + 1);
        } else {
            sm_x[k - 1] = sm_xnu[sm_indices[k]][2 * sm_argind[k] + 1];
            sm_fsum(k + 1);
            sm_x[k - 1] = 1.0 - sm_x[k - 1];
            sm_fsum(k + 1);
=======
            dummy = sm_fsum(k + 1);
=======
            sm_fsum(k + 1);
>>>>>>> 62f2bf2 (Refactor Smolyak algorithms and add Genz test suite)
        } else {
            sm_x[k - 1] = sm_xnu[sm_indices[k]][2 * sm_argind[k] + 1];
            sm_fsum(k + 1);
            sm_x[k - 1] = 1.0 - sm_x[k - 1];
<<<<<<< HEAD
            dummy = sm_fsum(k + 1);
>>>>>>> f26fbd8 (Implement Smolyak quadrature subsystem with Clenshaw-Curtis integration)
=======
            sm_fsum(k + 1);
>>>>>>> 62f2bf2 (Refactor Smolyak algorithms and add Genz test suite)
        }
    }
    return sm_ftotal;
}
