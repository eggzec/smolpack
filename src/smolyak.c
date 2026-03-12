/**
 * @file smolyak.c
 * @brief Delayed Clenshaw-Curtis Smolyak sparse-grid quadrature.
 *
 * Implements the Smolyak combination technique using a "delayed"
 * Clenshaw-Curtis basic sequence.  The delayed construction yields
 * fewer function evaluations than the standard Clenshaw-Curtis
 * sequence for the same level of polynomial exactness.
 *
 * @author  Knut Petras
 * @date    2007
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "smolpack.h"

/** Total number of nodes across all 1-D quadrature formulae. */
#define uniw 32

/** Number of distinct basic formulae in the delayed CC sequence. */
#define fn 6

/** Total number of basic formulae including multiplicities. */
#define gesfn 50

/* ------------------------------------------------------------------ */
/*                        File-scope state                            */
/* ------------------------------------------------------------------ */

static double quafo;                            /**< Accumulated cubature result.       */
static double x[maxdim];                        /**< Current function argument vector.  */
static double xnu[fn][uniw];                    /**< Quadrature node positions.         */
static double dnu[fn][uniw];                    /**< Delta (difference) weights.        */
static double ftotal;                           /**< Accumulated function sum.          */
static double wsum;                             /**< Accumulated weight sum.            */
static double summe;                            /**< Sub-formula accumulator.           */
static int    d;                                /**< Current spatial dimension.         */
static int    q;                                /**< Cubature formula parameter.        */
static int    n[fn];                            /**< Node counts per formula.           */
static int    ninv[fn];                         /**< Inverse node-count mapping.        */
static int    sw[gesfn];                        /**< Working array for formula indices. */
static int    wcount;                           /**< Weight-evaluation counter.         */

static int    indices[maxdim];                  /**< Formula indices per dimension.     */
static int    argind[maxdim];                   /**< Nodal indices per dimension.       */
static int    indsum[maxdim][maxdim];           /**< Divide-and-conquer parameters.     */

static int    anzw[uniw];                       /**< Node-frequency histogram.          */
static int    lookind[fn][uniw];                /**< Node-number look-up table.         */
static int    invlook[fn][uniw];                /**< Inverse look-up table.             */
static int    maxind;                           /**< Tree-depth parameter.              */
static int    wind[maxdim];                     /**< Auxiliary for slow coeff calc.     */

/* ------------------------------------------------------------------ */
/*                    Tree node definition                            */
/* ------------------------------------------------------------------ */

/**
 * @brief Binary-tree node used to cache computed coefficients.
 *
 * The tree stores previously computed Smolyak coefficients so that
 * repeated queries with the same node-frequency pattern can be
 * answered without redundant arithmetic.
 */
struct tnode {
    int            empty;       /**< Nonzero if this node has no data. */
    double        *coeff;       /**< Array of cached coefficients.     */
    int           *computed;    /**< Flags: 1 if coeff[i] is valid.    */
    struct tnode  *left;        /**< Left child pointer.               */
    struct tnode  *right;       /**< Right child pointer.              */
};

static struct tnode *root;              /**< Root of the coefficient tree. */

/** Pointer to the user-supplied integrand function. */
static double (*f)(int, double x[]);

/* ------------------------------------------------------------------ */
/*                    Forward declarations                            */
/* ------------------------------------------------------------------ */

static void            formula(int k, int l);
static double          eval(int k);
static double          fsum(int k);
static void            init(int dim);
static double          calccoeff(int l);
static double          calccoeff2(int k, int l);
static double          wl(int r, int s, int l);
static double          we(int r, int s, int l);
static void            sumind(int r, int s);
static double          coeff(void);
static struct tnode   *talloc(void);
static void            frei(struct tnode *p);
 
/******************************************************************************/
/**
 * @brief Approximate a multidimensional integral (delayed Clenshaw-Curtis).
 *
 * Constructs a sparse grid from a "delayed" Clenshaw-Curtis basic
 * sequence and applies the Smolyak combination technique.  The
 * delayed construction typically requires fewer function evaluations
 * than the standard CC rule for the same polynomial exactness.
 *
 * @param[in] dim          Spatial dimension (1 <= dim < maxdim).
 * @param[in] qq           Level parameter; k = qq - dim stages.
 * @param[in] ff           Integrand callback.
 * @param[in] print_stats  Nonzero to print evaluation statistics.
 * @return Approximated integral value.
 */
double int_smolyak(int dim, int qq, double (*ff)(int, double xx[]),
                   int print_stats) {
/* 
  Make the parameters global.
*/
  d = dim;
  q = qq;
  f = ff;
/* 
  Initialize.
*/
  wcount = 0;
  count = 0;
  quafo = 0.0; 

  init ( dim );
/* 
  Call the Smolyak algorithm.  Q-DIM = K, the number of stages.
*/
  formula ( 1, q-dim );
/* 
  Free the allocated memory.
*/
  frei ( root );
/* 
  Print statistics, if desired.
*/
  if ( print_stats ) {
    printf ( "\n" );
    printf ( "  Number of function calls =  %d\n", count );
    printf ( "  Weight evaluations =        %d\n", wcount );
  }

  return quafo;
}
/******************************************************************************/
/**
 * @brief Carry out the Smolyak combination algorithm.
 *
 * Recursively distributes the level-sum budget across dimensions.
 * When all dimensions have been assigned a formula index the
 * product-formula evaluation @c eval() is invoked.
 */
static void formula(int k, int l) {
  int i;

  if ( k == d + 1 ) {
    quafo = quafo + eval ( 0 );
  } else {
    for ( i = 0; i <= l; i++ ) {
      if ( sw[i] < fn ) {
        indices[k] = sw[i];
        formula ( k+1, l-i );
      }
    }
  }
  return;
}
/******************************************************************************/
/**
 * @brief Calculate the value of a product formula.
 *
 * Iterates over all node combinations for the current formula
 * assignment, multiplying the weight coefficient by the symmetrised
 * function sum.
 */
static double eval(int k) {
  double dummy;
  int i;

  if ( k == 0 ) {
    summe = 0.0;
    dummy = eval ( 1 );
  }
/* 
  Summation corresponding to one coefficient.
*/
  else if ( k == d + 1 ) {
    summe = summe + coeff ( ) * fsum ( 0 ); 
  }
/*
  Choice of the nodes.
*/
  else {
    for ( i = 0; i <= n[indices[k]]; i++ ) {
      argind[k] = i;
      dummy = eval ( k + 1 );
    }
  }
  return summe;
}
/******************************************************************************/
/**
 * @brief Return a Smolyak coefficient, caching results in the tree.
 *
 * Looks up the coefficient for the current node-frequency pattern
 * in the binary tree.  Computes and stores it on a cache miss.
 */
static double coeff(void) {
  int i;
  int j;
  int l;
  struct tnode *p;
  struct tnode *pt;
  int r;

  p = root;
/*
  Initialize.
*/
  for ( i = 0; i < maxind; i++ ) {
    anzw[i] = 0;
  }
/* 
  Frequency anzw[  ] of 1-dim nodes
*/
  for ( i = 1; i <= d; i++ ) {
    anzw[lookind[indices[i]][argind[i]]]++;
  }
/* 
  Search in the tree according to anzw[..]
*/
  for ( j = maxind-1; 1 <= j; j-- ) {
     
/* 
  anzw[j] to the LEFT.
  Generate the node if it does not exist.
*/ 
    if ( p->left == NULL ) {
      p->left = (struct tnode *) calloc ( maxdim, sizeof(struct tnode) );
      pt = (p->left+anzw[j]);
      pt->left = NULL;
      pt->right = NULL;
      pt->empty = 1;
      p = pt;
    } else {
      p = ( p->left+anzw[j] );
    }
/* 
  one to the RIGHT
*/ 
    if ( p->right == NULL ) {
      pt = talloc ( );
      pt->empty = 1;
      pt->left = NULL;
      pt->right = NULL;
/* 
  leaf with coefficient.
*/
      if ( j == 1 ) {
        pt->coeff = ( double * ) calloc ( maxdim, sizeof(double) );   
        pt->computed = ( int * ) calloc ( maxdim, sizeof(int) );
        pt->empty = 0;
      }
      p->right = pt;
    }
    p = p->right;
  }
/*
  The coefficient must be computed.
*/
  if ( !*( p->computed+anzw[0] ) ) {
    wcount++;
    *( p->coeff +anzw[0] ) = calccoeff ( q-d );
    *( p->computed+anzw[0] ) = 1;
  }  
  return *( p->coeff +anzw[0] );
}
/******************************************************************************/
/**
 * @brief Recursively free the coefficient-tree structure.
 *
 * Releases all dynamically allocated memory in the subtree
 * rooted at @p p, including coefficient arrays.
 */
static void frei(struct tnode *p) {
  int i;

  if ( !( p->empty ) ) {
    free ( p->coeff );
    free ( p->computed );
  }

  if ( !( p->left == NULL ) ) {
    for ( i = 0; i < maxind; i++ ) {
      if ( !( p->left+i == NULL ) ) {
        frei ( p->left +i );
      }
    }
    free ( p->left );
  }

  if ( !( p->right == NULL ) ) {
     frei ( p->right );
     free ( p->right );
  }
  return;
}
/******************************************************************************/
/**
 * @brief Allocate a single tree node.
 *
 * Returns a pointer to a freshly allocated @c tnode or @c NULL
 * on allocation failure.
 */
static struct tnode *talloc(void) {
  return ( struct tnode * ) malloc ( sizeof(struct tnode) );
}
/******************************************************************************/
/**
 * @brief Sum formula indices using divide-and-conquer.
 *
 * Recursively partitions the dimension range [r, s] and stores
 * partial sums in @c indsum[][] for use by the coefficient
 * calculators.
 */
static void sumind(int r, int s) {
  int q;
/*
  R == S, one dimensional.
  Do the calculation.
*/
  if ( r == s ) {
    indsum[r][s] = ninv[indices[r]];
  }
/*
  R < S, compute average Q, split to [R,Q] + [Q+1,S], 
  and call SUMIND recursively.
*/
  else {
    q = ( r + s ) / 2;
    sumind ( r, q );
    sumind ( q+1, s );
    indsum[r][s] = indsum[r][q] + indsum[q+1][s];
  }
  return;
}
/******************************************************************************/
/**
 * @brief Calculate a Smolyak coefficient via divide-and-conquer.
 *
 * Computes the subdivision parameters and delegates to @c wl.
 */
static double calccoeff(int l) {
  double value;
/*
  Calculate the subdivision parameters.
*/
  sumind ( 1, d );
/* 
  Start the divide and conquer process.
*/
  value = wl ( 1, d, l ); 

  return value;
}
/******************************************************************************/
/**
 * @brief Divide-and-conquer weight sum ("left" branch).
 *
 * Sums contributions in dimensions r through s whose formula
 * indices do not exceed l.
 */
static double wl(int r, int s, int l) {
  int i;
  int p;
  int q;
  double total;

  total = 0.0;
/*
  R = S, one dimensional.
  Do the calculation.
*/
  if ( r == s ) {
    p = lookind[indices[r]][argind[r]];

    for ( i = ninv[indices[r]]; i <= l; i++ ) {
      if ( sw[i] < fn ) {
        if ( i == 0 ) {
          total = total + dnu[0][0];
        } else {
          if ( indices[r] == 0 ) {
            total = total + dnu[sw[i]][0];
          } else {
            total = total + dnu[sw[i]][invlook[sw[i]][p]]; 
          }
        }
      }
    }
  }
/*
  R < S, compute average Q, split to [R,Q] + [Q+1,S] 
  and call WE and WL recursively.
*/
  else {
    q = ( r + s ) / 2;

    for ( i = indsum[r][q]; i <= l-indsum[q+1][s]; i++ ) {
      total = total + we ( r, q, i ) * wl ( q+1, s, l-i );
    }
  }

  return total;
}
/******************************************************************************/
/**
 * @brief Divide-and-conquer weight sum ("exact" branch).
 *
 * Sums contributions in dimensions r through s whose formula
 * indices exactly equal l (one-dimensional) or recursively
 * partitions the range.
 */
static double we(int r, int s, int l) {
  int i;
  int q;
  double total;

  total = 0.0;
/*
  R = S, one dimensional.
  Do the calculation.
*/
  if ( r == s ) {
    if ( sw[l] < fn ) {
      if ( sw[l] == 0 ) {
        total = dnu[0][0];
      } else {
        if ( indices[r] == 0 ) {
          total = dnu[sw[l]][0];
        } else {
          total = dnu[sw[l]][invlook[sw[l]][lookind[indices[r]][argind[r]]]];
        }
      }
    }
  }
/*
  R < S, compute average Q, split to [R,Q] + [Q+1,S] and call WE recursively.
*/
  else {
    q = ( r + s ) / 2;
    for ( i = indsum[r][q]; i <= l-indsum[q+1][s]; i++ ) {
      total = total + we ( r, q, i ) * we ( q+1, s, l-i );
    }
  }
  return total;
}
/******************************************************************************/
/**
 * @brief Calculate coefficient by brute-force enumeration (slow).
 *
 * Provided for reference and debugging.  The tree-based method
 * (@c coeff / @c calccoeff) is preferred in production.
 */
static double calccoeff2(int k, int l) {
  double dummy;
  int i;
  double wprod;

  if ( k == 0 ) {
    wcount++;
    wsum = 0.0;
    dummy = calccoeff2 ( 1, l );
  } else if ( k == d + 1 ) {
    wprod = 1.0;

    for ( i = 1; i <= d; i++ ) {
      if ( indices[i] == 0 ) {
        wprod = wprod * dnu[wind[i]][0];
      } else {
        wprod = wprod *
          dnu[wind[i]][invlook[wind[i]][lookind[indices[i]][argind[i]]]];
      }
    }
    wsum = wsum + wprod;
  } else {
    i = indices[k];
    while ( ninv[i] <= l ) {
      wind[k] = i;
      dummy = calccoeff2 ( k+1, l-ninv[i] );
      i++;
    }
  }
  return wsum;
}
/******************************************************************************/
/**
 * @brief Compute symmetric sums of integrand values.
 *
 * Evaluates the integrand at nodes for which the same weight
 * applies, exploiting symmetry: @c f(x) and @c f(1-x) share
 * a weight.
 */
static double fsum(int k) {
  double dummy;

  if ( k == 0 ) {
    ftotal = 0.0;
    dummy = fsum ( 1 );
  } else if ( k == d + 1 ) {
    ftotal = ftotal + (*f) ( d, x );
  } else {
    if ( indices[k] == 0 ) {
      x[k-1] = 0.5;
      dummy = fsum ( k+1 );
    }
/*
  Use symmetry to get both X and -X.
*/
    else {
      x[k-1] = xnu[indices[k]][2*argind[k]+1];
      dummy = fsum ( k+1 );
      x[k-1] = 1.0 - x[k-1];
      dummy = fsum ( k+1 );
    }  
  }
  return ftotal;
}
/******************************************************************************/
/**
 * @brief Initialise the delayed Clenshaw-Curtis quadrature data.
 *
 * Sets up the node positions, delta weights, look-up tables,
 * and the coefficient-tree root for a new integration call.
 */
static void init(int dim) {
  int formfakt;
  int freq[fn] = { 1, 2, 3, 6, 12, 24 };
  int i;
  int j;
  int maxform;
  int nj[fn] = { 1, 3, 7, 15, 31, 63 };
/* 
  Some parameter calculations.
*/
  n[0] = 0;
  j = 0;
  ninv[0] = 0;

  for ( i = 1; i < fn; i++ ) {
    n[i] = ( nj[i] - nj[i-1] - 2 ) / 2;
    j = j + freq[i-1];
    ninv[i] = j;
  }

  for ( i = 0; i < gesfn; i++ ) {
    sw[i] = fn;
  }

  for ( i = 0; i < fn; i++ ) {
    sw[ninv[i]] = i;
  }
/* 
  Number of the 'largest' used formula.
*/
  i = q - dim;
  maxform = 0;
  while ( ninv[maxform+1] <= i ) {
    maxform = maxform + 1;
  }
/* 
  Total number of used 1-dim nodes.
*/
  maxind = ( nj[maxform] + 1 ) / 2;
/* 
  Table of 1-dim nodal numbers 0..MAXIND-1 corresponding 
  to a combination formula number/nodal number
  and inverse formula.
*/
  lookind[0][0] = 0;
  for ( i = 1; i <= maxform; i++ ) {
    formfakt = 1 << (maxform - i);
    for ( j = 0; j < (nj[i]+1)/4; j++ ) {
      lookind[i][j] = formfakt * ( 2 * j + 1 );
    }
/* 
  In a linear ordering of all used nodes, the (2J+1)-th 
  node of the I-th basic formula is LOOKIND[I][J]-th node.
*/
    for ( j = 0; j < (nj[i]+1)/2; j++ ) {
      invlook[i][formfakt*j] = j;
    }
/* 
  The LOOKIND[I][2^(MAXFORM-I)]-th node in a linear 
  ordering of all used nodes is the J-th node of the 
  I-th basic formula.  Note that MAXFORM is the number
  of used different basic formulas.
*/ 
  }

/* 
  Root of the coefficient TREE.
*/
  root = talloc ( );
  root->empty = 1;
  root->left = NULL;
  root->right = NULL;
/* 
  One dimensional formulas (Deltas) 
*/
  xnu[0][0] = 0.5;
  dnu[0][0] = 1.0;

  xnu[1][0] =  5.000000000000000E-001;
  xnu[1][1] =  8.8729833462074168851793E-001;
  dnu[1][0] = -5.5555555555555555555556E-001;
  dnu[1][1] =  2.7777777777777777777778E-001;

  xnu[2][0] =  5.0000000000000000000000E-001;
  xnu[2][1] =  7.1712187467340127900104E-001;
  xnu[2][2] =  8.8729833462074168851793E-001;
  xnu[2][3] =  9.8024563435401014171175E-001;
  dnu[2][0] = -2.1898617511520737327189E-001;
  dnu[2][1] =  2.0069870738798111145253E-001;
  dnu[2][2] = -1.4353373284361105741349E-001;
  dnu[2][3] =  5.2328113013233632596912E-002;

  xnu[3][0] =  5.0000000000000000000000E-001;
  xnu[3][1] =  6.1169334321448344081410E-001;
  xnu[3][2] =  7.1712187467340127900104E-001;
  xnu[3][3] =  8.1055147336861320147034E-001;
  xnu[3][4] =  8.8729833462074168851793E-001;
  xnu[3][5] =  9.4422961643612849944521E-001;
  xnu[3][6] =  9.8024563435401014171175E-001;
  xnu[3][7] =  9.9691598160637751110426E-001;
  dnu[3][0] = -1.1270301943013372747934E-001;
  dnu[3][1] =  1.0957842920079374820185E-001;
  dnu[3][2] = -1.0038444269948660093556E-001;
  dnu[3][3] =  8.5755954568195690393677E-002;
  dnu[3][4] = -6.7036417312274610184300E-002;
  dnu[3][5] =  4.6463597657562268842947E-002;
  dnu[3][6] = -2.6526471514693762748452E-002;
  dnu[3][7] =  8.5008598149701301695137E-003;

  xnu[4][0] =  5.0000000000000000000000E-001;
  xnu[4][1] =  5.5624447156659331287292E-001;
  xnu[4][2] =  6.1169334321448344081410E-001;
  xnu[4][3] =  6.6556769662898841654632E-001;
  xnu[4][4] =  7.1712187467340127900104E-001;
  xnu[4][5] =  7.6565987182218781198605E-001;
  xnu[4][6] =  8.1055147336861320147034E-001;
  xnu[4][7] =  8.5124810324576353930490E-001;
  xnu[4][8] =  8.8729833462074168851793E-001;
  xnu[4][9] =  9.1836296908443436775138E-001;
  xnu[4][10] =  9.4422961643612849944521E-001;
  xnu[4][11] =  9.6482742871487002833506E-001;
  xnu[4][12] =  9.8024563435401014171175E-001;
  xnu[4][13] =  9.9076557477687005343368E-001;
  xnu[4][14] =  9.9691598160637751110426E-001;
  xnu[4][15] =  9.9954906248383379883111E-001;
  dnu[4][0] = -5.6377621538718997889636E-002;
  dnu[4][1] =  5.5978436510476728440072E-002;
  dnu[4][2] = -5.4789218672831429083502E-002;
  dnu[4][3] =  5.2834946790117404871908E-002;
  dnu[4][4] = -5.0157125382596721131319E-002;
  dnu[4][5] =  4.6813554990632236808329E-002;
  dnu[4][6] = -4.2877994543200514816583E-002;
  dnu[4][7] =  3.8439810249501765521353E-002;
  dnu[4][8] = -3.3603750473896758409784E-002;
  dnu[4][9] =  2.8489754747061678706099E-002;
  dnu[4][10] = -2.3232151026683275572245E-002;
  dnu[4][11] =  1.7978551653564661048389E-002;
  dnu[4][12] = -1.2897842450451543066137E-002;
  dnu[4][13] =  8.2230249271939054668942E-003;
  dnu[4][14] = -4.2835769453095770463562E-003;
  dnu[4][15] =  1.2723903957809372077014E-003;

  xnu[5][0] =  5.0000000000000000000000E-001;
  xnu[5][1] =  5.2817215652329639498598E-001;
  xnu[5][2] =  5.5624447156659331287292E-001;
  xnu[5][3] =  5.8411762577610373249116E-001;
  xnu[5][4] =  6.1169334321448344081410E-001;
  xnu[5][5] =  6.3887491101091215753268E-001;
  xnu[5][6] =  6.6556769662898841654632E-001;
  xnu[5][7] =  6.9167966209936517345824E-001;
  xnu[5][8] =  7.1712187467340127900104E-001;
  xnu[5][9] =  7.4180901347292051378108E-001;
  xnu[5][10] =  7.6565987182218781198605E-001;
  xnu[5][11] =  7.8859785502602290742185E-001;
  xnu[5][12] =  8.1055147336861320147034E-001;
  xnu[5][13] =  8.3145483001239029773051E-001;
  xnu[5][14] =  8.5124810324576353930490E-001;
  xnu[5][15] =  8.6987802217634737933861E-001;
  xnu[5][16] =  8.8729833462074168851793E-001;
  xnu[5][17] =  9.0347026597510880592815E-001;
  xnu[5][18] =  9.1836296908443436775138E-001;
  xnu[5][19] =  9.3195396909684523857321E-001;
  xnu[5][20] =  9.4422961643612849944521E-001;
  xnu[5][21] =  9.5518557847850214624890E-001;
  xnu[5][22] =  9.6482742871487002833506E-001;
  xnu[5][23] =  9.7317142918670145257425E-001;
  xnu[5][24] =  9.8024563435401014171175E-001;
  xnu[5][25] =  9.8609143737429089828903E-001;
  xnu[5][26] =  9.9076557477687005343368E-001;
  xnu[5][27] =  9.9434237877371473996926E-001;
  xnu[5][28] =  9.9691598160637751110426E-001;
  xnu[5][29] =  9.9860312968611097953823E-001;
  xnu[5][30] =  9.9954906248383379883111E-001;
  xnu[5][31] =  9.9993644406017880596898E-001;
  dnu[5][0] = -2.8188814180191987109744E-002;
  dnu[5][1] =  2.8138849915627150636298E-002;
  dnu[5][2] = -2.7989218255238568736295E-002;
  dnu[5][3] =  2.7740702178279681993919E-002;
  dnu[5][4] = -2.7394605263980886602235E-002;
  dnu[5][5] =  2.6952749667633031963438E-002;
  dnu[5][6] = -2.6417473395059144940870E-002;
  dnu[5][7] =  2.5791626976024229388405E-002;
  dnu[5][8] = -2.5078569652948020678807E-002;
  dnu[5][9] =  2.4282165203336599357974E-002;
  dnu[5][10] = -2.3406777495318230607005E-002;
  dnu[5][11] =  2.2457265826816098707127E-002;
  dnu[5][12] = -2.1438980012491308330637E-002;
  dnu[5][13] =  2.0357755058472159466947E-002;
  dnu[5][14] = -1.9219905124773999502032E-002;
  dnu[5][15] =  1.8032216390391286320054E-002;
  dnu[5][16] = -1.6801938573891486499334E-002;
  dnu[5][17] =  1.5536775555843982439942E-002;
  dnu[5][18] = -1.4244877374144904399846E-002;
  dnu[5][19] =  1.2934839663607373455379E-002;
  dnu[5][20] = -1.1615723310923858549074E-002;
  dnu[5][21] =  1.0297116957956355574594E-002;
  dnu[5][22] = -8.9892758695005258819409E-003;
  dnu[5][23] =  7.7033752332797489010654E-003;
  dnu[5][24] = -6.4518989979126939693347E-003;
  dnu[5][25] =  5.2491234548106609491364E-003;
  dnu[5][26] = -4.1115209485759406322653E-003;
  dnu[5][27] =  3.0577534110586231698391E-003;
  dnu[5][28] = -2.1084676488811257036154E-003;
  dnu[5][29] =  1.2895248973428441362139E-003;
  dnu[5][30] = -6.3981211766590320201509E-004;
  dnu[5][31] =  1.8161074092276532984679E-004;

  return;
}
