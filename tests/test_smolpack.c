/**
 * @file test_smolpack.c
 * @brief Test driver for SMOLPACK Smolyak cubature library.
 * @date  30 April 2007
 *
 * Exercises both the delayed and standard Clenshaw-Curtis Smolyak
 * algorithms against the Genz test integrand family.  Reports
 * PASS / FAIL for each test based on a tolerance threshold.
 *
 * Usage:
 *   ./test_smolpack [seed]
 *
 * @author  Knut Petras
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "smolpack.h"
#include "genz.h"

/* ------------------------------------------------------------------ */
/*                          Helpers                                   */
/* ------------------------------------------------------------------ */

/**
 * @brief Array of integrand function pointers, indexed 1..7.
 */
static double (*funcs[8])(int, double[]) = {
    NULL, f1, f2, f3, f4, f5, f6, f7
};

/**
 * @brief Human-readable names for the test integrands.
 */
static const char *func_names[8] = {
    NULL,
    "Oscillatory",
    "Product peak",
    "Corner peak",
    "Gaussian",
    "Continuous",
    "Discontinuous",
    "exp(sum(x))"
};

/**
 * @brief Run a single integration test.
 *
 * @param fnum       Function index (1..7).
 * @param dim        Spatial dimension.
 * @param k          Number of Smolyak stages.
 * @param bs         Basic sequence (1 = delayed CC, 2 = standard CC).
 * @param tol        Absolute error tolerance.
 * @param pass_cnt   Pointer to pass counter (incremented on success).
 * @param fail_cnt   Pointer to fail counter (incremented on failure).
 */
static void run_test(int fnum, int dim, int k, int bs, double tol,
                     int *pass_cnt, int *fail_cnt) {
    int q = dim + k;
    double quad;
    double exact;
    double err;
    const char *alg = (bs == 1) ? "delayed-CC" : "standard-CC";

    if (bs == 1) {
        quad = int_smolyak(dim, q, funcs[fnum], 0);
    } else {
        quad = cc_int_smolyak(dim, q, funcs[fnum], 0);
    }

    exact = integral(fnum, dim);
    err = fabs(quad - exact);

    if (err <= tol) {
        printf("  PASS  %-14s dim=%2d k=%d %s  err=%.3e\n",
               func_names[fnum], dim, k, alg, err);
        (*pass_cnt)++;
    } else {
        printf("  FAIL  %-14s dim=%2d k=%d %s  err=%.3e (tol=%.1e)\n",
               func_names[fnum], dim, k, alg, err, tol);
        (*fail_cnt)++;
    }
}

/* ------------------------------------------------------------------ */
/*                            main                                    */
/* ------------------------------------------------------------------ */

int main(int argc, char *argv[]) {
    int seed;
    int dim;
    int fnum;
    int bs;
    int i;
    int pass_cnt = 0;
    int fail_cnt = 0;
    double sum;

    /* Parse optional seed argument. */
    if (argc >= 2) {
        seed = atoi(argv[1]);
    } else {
        seed = 12345;
    }

    printf("\n");
    timestamp();
    printf("\nSMOLPACK Test Suite\n");
    printf("  Seed = %d\n\n", seed);

    /*
     * Test configuration: small dimensions, moderate stages.
     * Tolerance is generous because some Genz functions are hard.
     */
    int dims[] = { 2, 3, 5 };
    int ndims = 3;
    int stages[] = { 3, 5, 7 };
    int nstages = 3;
    double tol;
    double tol_default = 1.0e-2;
    double tol_discont = 2.0e-1; /* f6 is discontinuous — sparse grids converge slowly */

    /* Test all 7 Genz functions x dimensions x stages x both algorithms. */
    for (fnum = 1; fnum <= 7; fnum++) {
        /* Randomise parameters for this integrand. */
        srand((unsigned int)(seed + fnum));
        sum = 0.0;

        for (i = 0; i < 40; i++) {
            c[i] = ((double)rand()) / (double)RAND_MAX;
            sum = sum + c[i];
        }
        for (i = 0; i < 40; i++) {
            c[i] = c[i] / sum * 9.0;
        }
        for (i = 0; i < 40; i++) {
            w[i] = ((double)rand()) / (double)RAND_MAX;
        }

        printf("--- %s (f%d) ---\n", func_names[fnum], fnum);
        tol = (fnum == 6) ? tol_discont : tol_default;

        for (int di = 0; di < ndims; di++) {
            dim = dims[di];
            for (int si = 0; si < nstages; si++) {
                int k = stages[si];
                for (bs = 1; bs <= 2; bs++) {
                    run_test(fnum, dim, k, bs, tol, &pass_cnt, &fail_cnt);
                }
            }
        }
        printf("\n");
    }

    /* Summary. */
    printf("========================================\n");
    printf("  Total: %d passed, %d failed\n", pass_cnt, fail_cnt);
    printf("========================================\n\n");

    timestamp();

    return (fail_cnt > 0) ? 1 : 0;
}
