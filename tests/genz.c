/**
 * @file genz.c
 * @brief Genz test integrand family and support utilities.
 *
 * Contains the seven Genz test functions (f1–f7), exact-integral
 * computation, complementary error function, timestamp, and tuple
 * enumeration — all extracted from the original smolpack.c.
 *
 * @author  Knut Petras
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "genz.h"

/* ------------------------------------------------------------------ */
/*                   Global parameter definitions                     */
/* ------------------------------------------------------------------ */

double c[maxdim];
double w[maxdim];
double wp;
double derf[8][20];

/* ------------------------------------------------------------------ */
/*                        Genz integrands                             */
/* ------------------------------------------------------------------ */

/******************************************************************************/
/**
 * @brief Oscillatory: cos(2*pi*w[0] + sum(c[i]*x[i])).
 * @date  26 April 2007
 */
double f1(int dim, double x[]) {
    int i;
    double pi = 3.141592653589793;
    double value;

    count++;
    value = 2.0 * pi * w[0];
    for (i = 0; i < dim; i++) {
        value = value + c[i] * x[i];
    }
    value = cos(value);

    return value;
}

/******************************************************************************/
/**
 * @brief Product peak: 1 / prod((1/c[i])^2 + (x[i]-w[i])^2).
 * @date  26 April 2007
 */
double f2(int dim, double x[]) {
    int i;
    double value;

    count++;
    value = 1.0;

    for (i = 0; i < dim; i++) {
        value = value * (pow(1.0 / c[i], 2) + pow(x[i] - w[i], 2));
    }
    value = 1.0 / value;

    return value;
}

/******************************************************************************/
/**
 * @brief Corner peak: 1 / (1 + sum(c[i]*x[i]))^(dim+1).
 * @date  26 April 2007
 */
double f3(int dim, double x[]) {
    int i;
    double value;

    count++;
    value = 1.0;
    for (i = 0; i < dim; i++) {
        value = value + c[i] * x[i];
    }
    value = 1.0 / pow(value, dim + 1);

    return value;
}

/******************************************************************************/
/**
 * @brief Gaussian: exp(-sum((c[i]*(x[i]-w[i]))^2)).
 * @date  26 April 2007
 */
double f4(int dim, double x[]) {
    int i;
    double value;

    count++;
    value = 0.0;
    for (i = 0; i < dim; i++) {
        value = value - pow(c[i] * (x[i] - w[i]), 2);
    }
    value = exp(value);

    return value;
}

/******************************************************************************/
/**
 * @brief Continuous: exp(-sum(c[i]*|x[i]-w[i]|)).
 * @date  26 April 2007
 */
double f5(int dim, double x[]) {
    int i;
    double value;

    count++;
    value = 0.0;

    for (i = 0; i < dim; i++) {
        value = value - c[i] * fabs(x[i] - w[i]);
    }

    value = exp(value);

    return value;
}

/******************************************************************************/
/**
 * @brief Discontinuous: exp(sum(c[i]*x[i])) if inside box, else 0.
 * @date  30 April 2007
 */
double f6(int dim, double x[]) {
    int i;
    double value;

    count++;

    value = 0.0;

    if (dim == 1) {
        if (x[0] <= w[0]) {
            for (i = 0; i < dim; i++) {
                value = value + c[i] * x[i];
            }
            value = exp(value);
        }
    } else {
        if (x[0] <= w[0] && x[1] <= w[1]) {
            for (i = 0; i < dim; i++) {
                value = value + c[i] * x[i];
            }
            value = exp(value);
        }
    }

    return value;
}

/******************************************************************************/
/**
 * @brief Exponential sum: exp(sum(x[i])).
 * @date  26 April 2007
 */
double f7(int dim, double x[]) {
    int i;
    double value;

    count++;
    value = 0.0;
    for (i = 0; i < dim; i++) {
        value = value + x[i];
    }
    value = exp(value);

    return value;
}

/* ------------------------------------------------------------------ */
/*                       Support utilities                            */
/* ------------------------------------------------------------------ */

/******************************************************************************/
/**
 * @brief Initialise Taylor-coefficient table for @c my_erfc.
 * @date  26 April 2007
 */
void init_erf(void) {
    double hilf;
    int i;
    int j;
    double pi = 3.141592653589793;

    wp = sqrt(pi);
    hilf = 2.0 / wp;

    derf[0][0] = 1.0;
    derf[1][0] = 1.5729920705028513065877936491E-01;
    derf[2][0] = 4.6777349810472658379307436327E-03;
    derf[3][0] = 2.2090496998585441372776129582E-05;
    derf[4][0] = 1.5417257900280018852159673486E-08;
    derf[5][0] = 1.5374597944280348501883434853E-12;
    derf[6][0] = 2.1519736712498913116593350399E-17;
    derf[7][0] = 4.1838256077794143986140102238E-23;

    for (i = 0; i <= 7; i++) {
        derf[i][1] = -exp(-pow((double)i, 2)) * hilf;
        derf[i][2] = -(double)i * derf[i][1];

        for (j = 3; j <= 19; j++) {
            derf[i][j] = -2.0 * ((double)i * derf[i][j - 1]
                       + (double)(j - 2) * derf[i][j - 2]
                       / (double)(j - 1)) / (double)j;
        }
    }
    return;
}

/******************************************************************************/
/**
 * @brief Return the exact integral of Genz function @p fnum.
 * @date  30 April 2007
 *
 * @param fnum  Function index (1..7).
 * @param dim   Spatial dimension.
 * @return      Exact integral over [0,1]^dim.
 */
double integral(int fnum, int dim) {
    int a;
    double arg;
    double bot;
    double c_prod;
    int i;
    int *ivec;
    int ivec_sum;
    double pi = 3.141592653589793;
    double prod;
    int rank;
    double total;
    double value;

    /* #1: Oscillatory. */
    if (fnum == 1) {
        arg = 0.0;
        for (i = 0; i < dim; i++) {
            arg = arg + c[i];
        }

        prod = 1.0;
        for (i = 0; i < dim; i++) {
            prod = prod * sin(0.5 * c[i]) / c[i];
        }

        value = pow(2.0, dim) * cos(2.0 * pi * w[0] + 0.5 * arg) * prod;
    }
    /* #2: Product Peak. */
    else if (fnum == 2) {
        value = 1.0;
        for (i = 0; i < dim; i++) {
            value = value * c[i] * (atan(c[i] * (1.0 - w[i]))
                                  + atan(c[i] * w[i]));
        }
    }
    /* #3: Corner Peak. */
    else if (fnum == 3) {
        ivec = (int *)malloc(dim * sizeof(int));

        total = 0.0;
        rank = 0;

        for (;;) {
            tuple_next(0, 1, dim, &rank, ivec);

            if (rank == 0) {
                break;
            }

            ivec_sum = 0;
            for (i = 0; i < dim; i++) {
                ivec_sum = ivec_sum + ivec[i];
            }

            bot = 1.0;
            for (i = 0; i < dim; i++) {
                if (ivec[i] == 1) {
                    bot = bot + c[i];
                }
            }
            total = total + pow(-1.0, ivec_sum) / bot;
        }

        a = 1;
        for (i = 1; i <= dim; i++) {
            a = a * i;
        }

        c_prod = 1.0;
        for (i = 0; i < dim; i++) {
            c_prod = c_prod * c[i];
        }
        value = total / ((double)a * c_prod);

        free(ivec);
    }
    /* #4: Gaussian. */
    else if (fnum == 4) {
        init_erf();
        value = 1.0;
        for (i = 0; i < dim; i++) {
            value = value * sqrt(pi) / (2.0 * c[i])
                  * (my_erfc(-c[i] * w[i])
                   - my_erfc(c[i] * (1.0 - w[i])));
        }
    }
    /* #5: Continuous function. */
    else if (fnum == 5) {
        value = 1.0;
        for (i = 0; i < dim; i++) {
            value = value / c[i]
                * (2.0 - exp(-c[i] * w[i]) - exp(c[i] * (w[i] - 1.0)));
        }
    }
    /* #6: Discontinuous function. */
    else if (fnum == 6) {
        value = 1.0;

        if (dim < 2) {
            for (i = 0; i < dim; i++) {
                value = value * (exp(c[i] * w[i]) - 1.0) / c[i];
            }
        } else {
            for (i = 0; i <= 1; i++) {
                value = value * (exp(c[i] * w[i]) - 1.0) / c[i];
            }
            for (i = 2; i < dim; i++) {
                value = value * (exp(c[i]) - 1.0) / c[i];
            }
        }
    }
    /* #7: exp(sum(x(i))). */
    else if (fnum == 7) {
        value = pow(exp(1.0) - 1.0, dim);
    }
    /* Unexpected call! */
    else {
        value = 0.0;
        printf("\n");
        printf("INTEGRAL - Fatal error!\n");
        printf("  Input function index FNUM must be between 1 and 7.\n");
        printf("  This value was FNUM = %d\n", fnum);
        exit(1);
    }

    return value;
}

/******************************************************************************/
/**
 * @brief Evaluate the complementary error function erfc(x).
 * @date  26 April 2007
 *
 * @param x  Argument.
 * @return   erfc(x).
 */
double my_erfc(double x) {
    double hilf;
    int i;
    int j;
    double prod[60];
    double smme;
    double value;

    if (x < 0) {
        value = 2.0 - my_erfc(-x);
    } else if (x < 7.5) {
        i = (int)(x + 0.499999);
        double h = x - (double)i;
        value = h * derf[i][19] + derf[i][18];
        for (j = 17; 0 <= j; j--) {
            value = h * value + derf[i][j];
        }
    } else {
        hilf = 0.5 / pow(x, 2);
        prod[0] = 1.0;
        i = 2 + (int)(190.0 / x);
        for (j = 1; j <= i; j++) {
            prod[j] = -prod[j - 1] * (2 * j - 1) * hilf;
        }
        smme = 0.0;
        for (j = i; 0 <= j; j--) {
            smme = smme + prod[j];
        }
        value = (smme * exp(-pow(x, 2)) / (wp * x));
    }
    return value;
}

/******************************************************************************/
/**
 * @brief Print the current date and time as a timestamp.
 * @date  24 September 2003
 * @author John Burkardt
 */
void timestamp(void) {
#define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    size_t len;
    time_t now;

    now = time(NULL);
    tm = localtime(&now);

    len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);
    (void)len;

    printf("%s\n", time_buffer);

    return;
#undef TIME_SIZE
}

/******************************************************************************/
/**
 * @brief Compute the next element of a tuple space.
 * @date  29 April 2003
 * @author John Burkardt
 *
 * The elements are N-vectors with entries between M1 and M2,
 * produced in lexicographic order.
 *
 * @param m1    Minimum entry value.
 * @param m2    Maximum entry value.
 * @param n     Number of components.
 * @param rank  In/out rank counter (set to 0 on first call).
 * @param x     In/out tuple vector.
 */
void tuple_next(int m1, int m2, int n, int *rank, int x[]) {
    int i;
    int j;

    if (m2 < m1) {
        *rank = 0;
        return;
    }

    if (*rank <= 0) {
        for (i = 0; i < n; i++) {
            x[i] = m1;
        }
        *rank = 1;
    } else {
        *rank = *rank + 1;
        i = n - 1;

        for (;;) {
            if (x[i] < m2) {
                x[i] = x[i] + 1;
                break;
            }

            x[i] = m1;

            if (i == 0) {
                *rank = 0;
                for (j = 0; j < n; j++) {
                    x[j] = m1;
                }
                break;
            }
            i = i - 1;
        }
    }
    return;
}
