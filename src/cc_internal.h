/**
 * @file cc_internal.h
 * @brief Internal shared state for the standard Clenshaw-Curtis Smolyak
 *        subsystem.
 *
 * This header is **private** to the three translation units that form
 * the standard-CC implementation:
 *
 *   - cc_init.c     — quadrature data tables and cc_init()
 *   - cc_weights.c  — coefficient and weight computation
 *   - cc_alg.c      — Smolyak combination algorithm and public API
 *
 * @warning Do NOT include this header from outside the cc_*
 *          translation units.
 */

#ifndef CC_INTERNAL_H
#define CC_INTERNAL_H

#include "smolpack.h"
#include <stdlib.h>  /* qsort */

/* ------------------------------------------------------------------ */
/*                       Subsystem constants                          */
/* ------------------------------------------------------------------ */

/** Maximum node slots per 1-D formula (standard CC: up to 1025). */
#define CC_UNIW   1025

/** Number of distinct 1-D base formulae in the standard CC sequence. */
#define CC_FN      12

/** Total formula slots (standard CC: freq[i] = 1 for all i). */
#define CC_GESFN   12

/* ------------------------------------------------------------------ */
/*                    Coefficient-tree node type                      */
/* ------------------------------------------------------------------ */

/**
 * @brief Binary-tree node for caching standard-CC Smolyak coefficients.
 *
 * Unlike the delayed-CC tree, each node stores a single scalar
 * coefficient and a single computed flag instead of arrays.
 */
struct cc_node {
    double          coeff;     /**< Cached coefficient value.         */
    int             computed;  /**< 1 when @c coeff has been set.     */
    struct cc_node *left;      /**< Left child pointer.               */
    struct cc_node *right;     /**< Right child pointer.              */
};

/* ------------------------------------------------------------------ */
/*             Shared state — defined in cc_init.c                    */
/* ------------------------------------------------------------------ */

extern double cc_quafo;                           /**< Accumulated cubature result.       */
extern double cc_x[maxdim];                       /**< Current function argument vector.  */
extern double cc_xnu[CC_FN][CC_UNIW];             /**< 1-D quadrature node positions.     */
extern double cc_dnu[CC_FN][CC_UNIW];             /**< 1-D delta difference weights.      */
extern double cc_ftotal;                          /**< Accumulated function-value sum.    */
extern double cc_wsum;                            /**< Coefficient accumulator.           */
extern double cc_summe;                           /**< Sub-formula value accumulator.     */
extern int    cc_d;                               /**< Current spatial dimension.         */
extern int    cc_q;                               /**< Cubature level parameter.          */
extern int    cc_n[CC_FN];                        /**< Half node-count per formula.       */
extern int    cc_ninv[CC_FN];                     /**< Inverse node-count mapping.        */
extern int    cc_sw[CC_GESFN];                    /**< Formula-index working array.       */
extern int    cc_wcount;                          /**< Weight-evaluation counter.         */
extern int    cc_indices[maxdim];                 /**< Formula index assigned per dim.    */
extern int    cc_argind[maxdim];                  /**< Node index assigned per dim.       */
extern int    cc_indsum[maxdim][maxdim];          /**< Partial index sums for D&C.        */
extern int    cc_anzw[CC_UNIW];                   /**< Node-frequency histogram.          */
extern int    cc_nodes_arr[maxdim];               /**< Sorted unique node indices.        */
extern int    cc_lookind[CC_FN][CC_UNIW];         /**< Formula/node → linear index table. */
extern int    cc_invlook[CC_FN][CC_UNIW];         /**< Linear index → node inverse table. */
extern int    cc_maxind;                          /**< Tree depth parameter.              */
extern int    cc_wind[maxdim];                    /**< Auxiliary array for slow coeff.    */

extern struct cc_node         *cc_root;           /**< Root of the coefficient tree.      */
extern double (*cc_f)(int, double x[]);           /**< User integrand function pointer.   */

/* ------------------------------------------------------------------ */
/*                       Internal prototypes                          */
/* ------------------------------------------------------------------ */

/* cc_init.c */
void             cc_setup(int dim);

/* cc_weights.c */
struct cc_node  *cc_talloc(void);
void             cc_frei(struct cc_node *p);
int              cc_less(const void *a, const void *b);
double           cc_coeff(void);
double           cc_calccoeff(int l);
double           cc_calccoeff2(int k, int l);
void             cc_sumind(int r, int s);
double           cc_wl(int r, int s, int l);
double           cc_we(int r, int s, int l);

/* cc_alg.c */
void             cc_formula(int k, int l);
double           cc_eval(int k);
double           cc_fsum(int k);

#endif /* CC_INTERNAL_H */
