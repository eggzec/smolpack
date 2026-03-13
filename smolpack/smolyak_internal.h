/**
 * @file smolyak_internal.h
 * @brief Internal shared state for the delayed Clenshaw-Curtis Smolyak
 *        subsystem.
 *
 * This header is **private** to the three translation units that form
 * the delayed-CC implementation:
 *
 *   - smolyak_init.c    — quadrature data tables and sm_init()
 *   - smolyak_weights.c — coefficient and weight computation
 *   - smolyak_alg.c     — Smolyak combination algorithm and public API
 *
 * It declares the subsystem constants, the coefficient-tree node type,
 * all shared state variables (defined once in smolyak_init.c), and
 * the internal function prototypes.
 *
 * @warning Do NOT include this header from outside the smolyak
 *          translation units.
 */

#ifndef SMOLYAK_INTERNAL_H
#define SMOLYAK_INTERNAL_H

#include "smolpack.h"

/* ------------------------------------------------------------------ */
/*                       Subsystem constants                          */
/* ------------------------------------------------------------------ */

/** Total node slots per 1-D formula (delayed CC uses at most 32). */
#define SM_UNIW   32

/** Number of distinct 1-D base formulae in the delayed CC sequence. */
#define SM_FN      6

/** Total formula slots including multiplicities (= sum of freq[]). */
#define SM_GESFN  50

/* ------------------------------------------------------------------ */
/*                    Coefficient-tree node type                      */
/* ------------------------------------------------------------------ */

/**
 * @brief Binary-tree node used to cache computed Smolyak coefficients.
 *
 * The tree is navigated by the node-frequency histogram (@c sm_anzw).
 * Each leaf stores an array of coefficients and corresponding validity
 * flags, one entry per possible value of @c sm_anzw[0].
 */
struct sm_node {
    int             empty;     /**< Non-zero when this node holds no data. */
    double         *coeff;     /**< Array of cached coefficient values.    */
    int            *computed;  /**< Per-entry validity flags (1 = valid).  */
    struct sm_node *left;      /**< Left child (indexed by node frequency).*/
    struct sm_node *right;     /**< Right child.                           */
};

/* ------------------------------------------------------------------ */
/*             Shared state — defined in smolyak_init.c               */
/* ------------------------------------------------------------------ */

extern double sm_quafo;                          /**< Accumulated cubature result.        */
extern double sm_x[maxdim];                      /**< Current function argument vector.   */
extern double sm_xnu[SM_FN][SM_UNIW];            /**< 1-D quadrature node positions.      */
extern double sm_dnu[SM_FN][SM_UNIW];            /**< 1-D delta difference weights.       */
extern double sm_ftotal;                         /**< Accumulated function-value sum.     */
extern double sm_wsum;                           /**< Coefficient accumulator.            */
extern double sm_summe;                          /**< Sub-formula value accumulator.      */
extern int    sm_d;                              /**< Current spatial dimension.          */
extern int    sm_q;                              /**< Cubature level parameter.           */
extern int    sm_n[SM_FN];                       /**< Half node-count per formula.        */
extern int    sm_ninv[SM_FN];                    /**< Inverse node-count mapping.         */
extern int    sm_sw[SM_GESFN];                   /**< Formula-index working array.        */
extern int    sm_wcount;                         /**< Weight-evaluation counter.          */
extern int    sm_indices[maxdim];                /**< Formula index assigned per dim.     */
extern int    sm_argind[maxdim];                 /**< Node index assigned per dim.        */
extern int    sm_indsum[maxdim][maxdim];         /**< Partial index sums for D&C.         */
extern int    sm_anzw[SM_UNIW];                  /**< Node-frequency histogram.           */
extern int    sm_lookind[SM_FN][SM_UNIW];        /**< Formula/node → linear index table.  */
extern int    sm_invlook[SM_FN][SM_UNIW];        /**< Linear index → node inverse table.  */
extern int    sm_maxind;                         /**< Tree depth parameter.               */
extern int    sm_wind[maxdim];                   /**< Auxiliary array for slow coeff.     */

extern struct sm_node         *sm_root;          /**< Root of the coefficient tree.       */
extern double (*sm_f)(int, double x[]);          /**< User integrand function pointer.    */

/* ------------------------------------------------------------------ */
/*                       Internal prototypes                          */
/* ------------------------------------------------------------------ */

/* smolyak_init.c */
void            sm_init(int dim);

/* smolyak_weights.c */
struct sm_node *sm_talloc(void);
void            sm_frei(struct sm_node *p);
double          sm_coeff(void);
double          sm_calccoeff(int l);
double          sm_calccoeff2(int k, int l);
void            sm_sumind(int r, int s);
double          sm_wl(int r, int s, int l);
double          sm_we(int r, int s, int l);

/* smolyak_alg.c */
void            sm_formula(int k, int l);
double          sm_eval(int k);
double          sm_fsum(int k);

#endif /* SMOLYAK_INTERNAL_H */
