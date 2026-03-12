/**
 * @file smolyak_weights.c
 * @brief Coefficient computation and weight functions for Clenshaw-Curtis
 *        Smolyak quadrature.
 */

#include <stdio.h>
#include <stdlib.h>
#include "smolyak_internal.h"

/**
 * @brief Allocate a single sm_node.
 *
 * @return Pointer to the newly allocated sm_node.
 */
struct sm_node *sm_talloc(void) {
    return (struct sm_node *) malloc(sizeof(struct sm_node));
}

/**
 * @brief Recursively free a tree of sm_node structures.
 *
 * @param p Pointer to the root sm_node of the subtree to free.
 */
void sm_frei(struct sm_node *p) {
    int i;
    if (!(p->empty)) {
        free(p->coeff);
        free(p->computed);
    }
    if (!(p->left == NULL)) {
        for (i = 0; i < sm_maxind; i++) {
            if (!(p->left + i == NULL)) {
                sm_frei(p->left + i);
            }
        }
        free(p->left);
    }
    if (!(p->right == NULL)) {
        sm_frei(p->right);
        free(p->right);
    }
}

/**
 * @brief Return the Smolyak coefficient, caching results in a tree.
 *
 * @return The computed coefficient value.
 */
double sm_coeff(void) {
    int i;
    int j;
    int l;
    struct sm_node *p;
    struct sm_node *pt;
    int r;

    p = sm_root;

    for (i = 0; i < sm_maxind; i++) {
        sm_anzw[i] = 0;
    }

    for (i = 1; i <= sm_d; i++) {
        sm_anzw[sm_lookind[sm_indices[i]][sm_argind[i]]]++;
    }

    for (j = sm_maxind - 1; 1 <= j; j--) {
        if (p->left == NULL) {
            p->left = (struct sm_node *) calloc(maxdim, sizeof(struct sm_node));
            pt = (p->left + sm_anzw[j]);
            pt->left = NULL;
            pt->right = NULL;
            pt->empty = 1;
            p = pt;
        } else {
            p = (p->left + sm_anzw[j]);
        }

        if (p->right == NULL) {
            pt = sm_talloc();
            pt->empty = 1;
            pt->left = NULL;
            pt->right = NULL;
            if (j == 1) {
                pt->coeff = (double *) calloc(maxdim, sizeof(double));
                pt->computed = (int *) calloc(maxdim, sizeof(int));
                pt->empty = 0;
            }
            p->right = pt;
        }
        p = p->right;
    }

    if (!*(p->computed + sm_anzw[0])) {
        sm_wcount++;
        *(p->coeff + sm_anzw[0]) = sm_calccoeff(sm_q - sm_d);
        *(p->computed + sm_anzw[0]) = 1;
    }
    return *(p->coeff + sm_anzw[0]);
}

/**
 * @brief Divide-and-conquer computation of index sums.
 *
 * @param r Left index of the range.
 * @param s Right index of the range.
 */
void sm_sumind(int r, int s) {
    int q;
    if (r == s) {
        sm_indsum[r][s] = sm_ninv[sm_indices[r]];
    } else {
        q = (r + s) / 2;
        sm_sumind(r, q);
        sm_sumind(q + 1, s);
        sm_indsum[r][s] = sm_indsum[r][q] + sm_indsum[q + 1][s];
    }
}

/**
 * @brief Calculate the Smolyak coefficient via divide-and-conquer.
 *
 * @param l Level parameter.
 * @return The computed coefficient value.
 */
double sm_calccoeff(int l) {
    double value;
    sm_sumind(1, sm_d);
    value = sm_wl(1, sm_d, l);
    return value;
}

/**
 * @brief Divide-and-conquer weight computation ("left" branch).
 *
 * @param r Left index of the range.
 * @param s Right index of the range.
 * @param l Level parameter.
 * @return The accumulated weight total.
 */
double sm_wl(int r, int s, int l) {
    int i;
    int p;
    int q;
    double total;

    total = 0.0;

    if (r == s) {
        p = sm_lookind[sm_indices[r]][sm_argind[r]];
        for (i = sm_ninv[sm_indices[r]]; i <= l; i++) {
            if (sm_sw[i] < SM_FN) {
                if (i == 0) {
                    total = total + sm_dnu[0][0];
                } else {
                    if (sm_indices[r] == 0) {
                        total = total + sm_dnu[sm_sw[i]][0];
                    } else {
                        total = total + sm_dnu[sm_sw[i]][sm_invlook[sm_sw[i]][p]];
                    }
                }
            }
        }
    } else {
        q = (r + s) / 2;
        for (i = sm_indsum[r][q]; i <= l - sm_indsum[q + 1][s]; i++) {
            total = total + sm_we(r, q, i) * sm_wl(q + 1, s, l - i);
        }
    }
    return total;
}

/**
 * @brief Divide-and-conquer weight computation ("exact" branch).
 *
 * @param r Left index of the range.
 * @param s Right index of the range.
 * @param l Level parameter.
 * @return The accumulated weight total.
 */
double sm_we(int r, int s, int l) {
    int i;
    int q;
    double total;

    total = 0.0;

    if (r == s) {
        if (sm_sw[l] < SM_FN) {
            if (sm_sw[l] == 0) {
                total = sm_dnu[0][0];
            } else {
                if (sm_indices[r] == 0) {
                    total = sm_dnu[sm_sw[l]][0];
                } else {
                    total = sm_dnu[sm_sw[l]][sm_invlook[sm_sw[l]][sm_lookind[sm_indices[r]][sm_argind[r]]]];
                }
            }
        }
    } else {
        q = (r + s) / 2;
        for (i = sm_indsum[r][q]; i <= l - sm_indsum[q + 1][s]; i++) {
            total = total + sm_we(r, q, i) * sm_we(q + 1, s, l - i);
        }
    }
    return total;
}

/**
 * @brief Brute-force coefficient computation.
 *
 * @param k Current dimension index (0 to start, sm_d+1 to accumulate).
 * @param l Remaining level budget.
 * @return The accumulated weighted sum.
 */
double sm_calccoeff2(int k, int l) {
    double dummy;
    int i;
    double wprod;

    if (k == 0) {
        sm_wcount++;
        sm_wsum = 0.0;
        dummy = sm_calccoeff2(1, l);
    } else if (k == sm_d + 1) {
        wprod = 1.0;
        for (i = 1; i <= sm_d; i++) {
            if (sm_indices[i] == 0) {
                wprod = wprod * sm_dnu[sm_wind[i]][0];
            } else {
                wprod = wprod *
                    sm_dnu[sm_wind[i]][sm_invlook[sm_wind[i]][sm_lookind[sm_indices[i]][sm_argind[i]]]];
            }
        }
        sm_wsum = sm_wsum + wprod;
    } else {
        i = sm_indices[k];
        while (sm_ninv[i] <= l) {
            sm_wind[k] = i;
            dummy = sm_calccoeff2(k + 1, l - sm_ninv[i]);
            i++;
        }
    }
    return sm_wsum;
}
