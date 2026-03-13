/**
 * @file cc_weights.c
 * @brief Coefficient computation and weight functions for standard CC Smolyak.
 *
 * Contains the tree-based coefficient caching, memory management,
 * divide-and-conquer weight evaluation, and the qsort comparator.
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
 * @brief Allocate and return a new tree node.
 * @date  25 April 2007
 *
 * @return Pointer to the newly allocated tree node.
 */
struct cc_node *cc_talloc(void) {
    return (struct cc_node *)malloc(sizeof(struct cc_node));
}

/******************************************************************************/
/**
 * @brief Free the binary tree of cached coefficients.
 * @date  25 April 2007
 *
 * @param[in] p  Pointer to the tree node to free (recursively).
 */
void cc_frei(struct cc_node *p) {
    int i;

    if (!(p->left == NULL)) {
        for (i = 0; i < cc_maxind; i++) {
            cc_frei(p->left + i);
        }
        free(p->left);
    }

    if (!(p->right == NULL)) {
        for (i = 0; i <= cc_d; i++) {
            cc_frei(p->right + i);
        }
        free(p->right);
    }
}

/******************************************************************************/
/**
 * @brief Integer comparison callback for @c qsort.
 * @date  25 April 2007
 *
 * @param[in] a_ptr  Pointer to first integer (as @c const @c void*).
 * @param[in] b_ptr  Pointer to second integer.
 * @return           Negative, zero, or positive as @c *a < , == , > @c *b.
 */
int cc_less(const void *a_ptr, const void *b_ptr) {
    const int *a = (const int *)a_ptr;
    const int *b = (const int *)b_ptr;
    return (*a - *b);
}

/******************************************************************************/
/**
 * @brief Return the Smolyak coefficient, computing and caching as needed.
 * @date  25 April 2007
 *
 * Navigates the binary tree keyed by sorted node indices and their
 * frequencies.  If the coefficient at the leaf has not yet been computed,
 * it is evaluated via @c cc_calccoeff and cached.
 *
 * @return The coefficient value (retrieved from tree or freshly computed).
 */
double cc_coeff(void) {
    int i;
    int j;
    int num;
    struct cc_node *p;
    struct cc_node *pt;

    p = cc_root;
    /*
     * Initialize.
     */
    num = -1;
    /*
     * Frequency cc_anzw[] of 1-dim nodes.
     */
    for (i = 1; i <= cc_d; i++) {
        if (++cc_anzw[cc_lookind[cc_indices[i]][cc_argind[i]]] == 1) {
            cc_nodes_arr[++num] = cc_lookind[cc_indices[i]][cc_argind[i]];
        }
    }
    qsort(cc_nodes_arr, num + 1, sizeof(int), cc_less);
    /*
     * Search in the tree according to cc_anzw[..].
     */
    for (j = num; 0 <= j; j--) {
        /*
         * nodes_arr[j] to the LEFT.
         * Generate the node if it does not exist.
         */
        if (p->left == NULL) {
            p->left = (struct cc_node *)calloc(cc_maxind, sizeof(struct cc_node));
            pt = (p->left + cc_nodes_arr[j]);
            pt->left = NULL;
            pt->right = NULL;
            p = pt;
        } else {
            p = (p->left + cc_nodes_arr[j]);
        }
        /*
         * cc_anzw[nodes_arr[j]] to the RIGHT.
         */
        if (p->right == NULL) {
            p->right = (struct cc_node *)calloc(cc_d + 1, sizeof(struct cc_node));
            pt = (p->right + cc_anzw[cc_nodes_arr[j]]);
            pt->left = NULL;
            pt->right = NULL;
            p = pt;
        } else {
            p = (p->right + cc_anzw[cc_nodes_arr[j]]);
        }
        cc_anzw[cc_nodes_arr[j]] = 0;
    }
    /*
     * The coefficient must be calculated.
     */
    if (!(p->computed)) {
        (p->coeff) = cc_calccoeff(cc_q - cc_d);
        (p->computed) = 1;
    }

    return (p->coeff);
}

/******************************************************************************/
/**
 * @brief Sum formula indices across dimensions @p r through @p s.
 * @date  25 April 2007
 *
 * @param r  Start dimension.
 * @param s  End dimension.
 */
void cc_sumind(int r, int s) {
    int q;

    /*
     * R == S, one dimensional — do the calculation.
     */
    if (r == s) {
        cc_indsum[r][s] = cc_ninv[cc_indices[r]];
    } else {
        /*
         * R < S: compute average Q, split to [R,Q] + [Q+1,S],
         * and call cc_sumind recursively.
         */
        q = (r + s) / 2;
        cc_sumind(r, q);
        cc_sumind(q + 1, s);
        cc_indsum[r][s] = cc_indsum[r][q] + cc_indsum[q + 1][s];
    }
    return;
}

/******************************************************************************/
/**
 * @brief Calculate a Smolyak coefficient via divide-and-conquer.
 * @date  25 April 2007
 *
 * @param l  Index-sum budget.
 * @return   The computed coefficient.
 */
double cc_calccoeff(int l) {
    double value;

    cc_wcount++;
    /*
     * Calculate the subdivision parameters.
     */
    cc_sumind(1, cc_d);
    /*
     * Start the divide and conquer process.
     */
    value = cc_wl(1, cc_d, l);

    return value;
}

/******************************************************************************/
/**
 * @brief Weighted sum over dimensions @p r..@p s with index-sum budget @p l.
 * @date  25 April 2007
 *
 * Uses a divide-and-conquer split; left half evaluated with @c cc_we.
 *
 * @param r  Start dimension.
 * @param s  End dimension.
 * @param l  Remaining index budget.
 * @return   Accumulated weight-times-value product.
 */
double cc_wl(int r, int s, int l) {
    int i;
    int p;
    int q;
    double total;

    total = 0.0;
    /*
     * R = S, one dimensional — do the calculation.
     */
    if (r == s) {
        p = cc_lookind[cc_indices[r]][cc_argind[r]];

        for (i = cc_ninv[cc_indices[r]]; i <= l; i++) {
            if (cc_sw[i] < CC_FN) {
                if (i == 0) {
                    total = total + cc_dnu[0][0];
                } else {
                    if (cc_indices[r] == 0) {
                        total = total + cc_dnu[cc_sw[i]][0];
                    } else {
                        total = total + cc_dnu[cc_sw[i]][cc_invlook[cc_sw[i]][p]];
                    }
                }
            }
        }
    } else {
        /*
         * R < S: split to [R,Q] + [Q+1,S] and call cc_we and cc_wl recursively.
         */
        q = (r + s) / 2;
        for (i = cc_indsum[r][q]; i <= l - cc_indsum[q + 1][s]; i++) {
            total = total + cc_we(r, q, i) * cc_wl(q + 1, s, l - i);
        }
    }
    return total;
}

/******************************************************************************/
/**
 * @brief Exact-budget weighted sum over dimensions @p r..@p s.
 * @date  25 April 2007
 *
 * Like @c cc_wl but the index sum must equal @p l exactly.
 *
 * @param r  Start dimension.
 * @param s  End dimension.
 * @param l  Exact index budget.
 * @return   Product of one-dimensional weights.
 */
double cc_we(int r, int s, int l) {
    int i;
    int q;
    double total;

    total = 0.0;
    /*
     * R = S, one dimensional — do the calculation.
     */
    if (r == s) {
        if (cc_sw[l] < CC_FN) {
            if (cc_sw[l] == 0) {
                total = cc_dnu[0][0];
            } else {
                if (cc_indices[r] == 0) {
                    total = cc_dnu[cc_sw[l]][0];
                } else {
                    total = cc_dnu[cc_sw[l]][cc_invlook[cc_sw[l]][cc_lookind[cc_indices[r]][cc_argind[r]]]];
                }
            }
        }
    } else {
        /*
         * R < S: split to [R,Q] + [Q+1,S] and call cc_we recursively.
         */
        q = (r + s) / 2;
        for (i = cc_indsum[r][q]; i <= l - cc_indsum[q + 1][s]; i++) {
            total = total + cc_we(r, q, i) * cc_we(q + 1, s, l - i);
        }
    }
    return total;
}

/******************************************************************************/
/**
 * @brief Calculate coefficients via exhaustive enumeration (slow).
 * @date  25 April 2007
 *
 * @param k  Recursion depth.
 * @param l  Remaining index budget.
 * @return   Accumulated weight sum.
 */
double cc_calccoeff2(int k, int l) {
    int i;
    double wprod;

    if (k == 0) {
        cc_wcount++;
        cc_wsum = 0.0;
        cc_calccoeff2(1, l);
    } else if (k == cc_d + 1) {
        wprod = 1.0;
        for (i = 1; i <= cc_d; i++) {
            if (cc_indices[i] == 0) {
                wprod = wprod * cc_dnu[cc_wind[i]][0];
            } else {
                wprod = wprod *
                    cc_dnu[cc_wind[i]][cc_invlook[cc_wind[i]][cc_lookind[cc_indices[i]][cc_argind[i]]]];
            }
        }
        cc_wsum = cc_wsum + wprod;
    } else {
        i = cc_indices[k];
        while (cc_ninv[i] <= l) {
            cc_wind[k] = i;
            cc_calccoeff2(k + 1, l - cc_ninv[i]);
            i++;
        }
    }
    return cc_wsum;
}
