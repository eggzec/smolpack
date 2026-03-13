# Theory

This page covers the mathematical foundations behind `smolpack`: what
multidimensional integration is, why tensor-product rules fail at scale, and
how Smolyak's algorithm constructs sparse grids that overcome the curse of
dimensionality.

---

## 1) Background: Multidimensional integration

**Cubature** is the numerical approximation of a multidimensional integral:

$$
I[f] = \int_{[0,1]^d} f(\mathbf{x})\,d\mathbf{x},
$$

where $f: [0,1]^d \to \mathbb{R}$ is an integrand defined on the unit
hypercube. In one dimension this is simply **quadrature** — a well-understood
problem with efficient classical methods (Gauss, Clenshaw-Curtis, etc.).
Extending these methods to $d$ dimensions is the central challenge.

Cubature arises constantly in computational science and engineering:

- **Physics** — partition functions, path integrals, radiative transfer
- **Finance** — option pricing, risk quantification, portfolio optimization
- **Statistics** — Bayesian posterior expectations, marginal likelihoods
- **Engineering** — uncertainty quantification, reliability analysis

### The curse of dimensionality

The simplest approach to multidimensional integration is the **tensor-product
rule**: take a 1-D quadrature rule with $n$ nodes and form the Cartesian
product across all $d$ dimensions. This produces $n^d$ nodes — an exponential
explosion that quickly becomes infeasible:

| Dimension $d$ | Nodes ($n=10$) |
|---|---|
| 1 | 10 |
| 3 | 1,000 |
| 5 | 100,000 |
| 10 | 10,000,000,000 |

Even for moderate $d$, tensor-product rules are computationally intractable.
This exponential scaling is the **curse of dimensionality**, first identified
by Richard Bellman in the context of dynamic programming.

## 2) Smolyak's algorithm

In 1963, Sergey Smolyak proposed a construction that achieves polynomial
exactness comparable to tensor-product rules while requiring dramatically
fewer nodes. The key idea is to take a **carefully chosen linear combination**
of lower-dimensional tensor-product rules rather than a single high-dimensional
one.

### Construction

Given a family of 1-D quadrature rules $U^{(1)}, U^{(2)}, U^{(3)}, \ldots$
with increasing accuracy (called the **basic sequence**), the Smolyak formula
of level $q$ in $d$ dimensions is:

$$
A(q, d) = \sum_{q-d+1 \le |\mathbf{i}| \le q}
           (-1)^{q - |\mathbf{i}|}
           \binom{d-1}{q - |\mathbf{i}|}
           \left(U^{(i_1)} \otimes U^{(i_2)} \otimes \cdots \otimes U^{(i_d)}\right),
$$

where $|\mathbf{i}| = i_1 + i_2 + \cdots + i_d$ and $\otimes$ denotes the
tensor product.

This alternating-sign combination cancels out lower-order errors while
retaining high-order polynomial exactness. The resulting set of nodes forms
a **sparse grid** — a subset of the full tensor-product grid with far fewer
points.

### Node count comparison

For smooth integrands, the Smolyak formula with $q - d = k$ stages uses:

$$
N_{\text{sparse}} = O\!\left(\frac{(\log N)^{d-1}}{(d-1)!}\right)
$$

nodes, compared to $N_{\text{tensor}} = O(n^d)$ for the tensor-product rule
achieving the same polynomial degree. The savings grow exponentially with
dimension.

!!! info "Polynomial exactness"
    A Smolyak formula of level $q$ exactly integrates all polynomials of
    total degree at most $2(q - d) + 1$ when using Clenshaw-Curtis basic
    rules. This is precisely the same degree of exactness that the
    corresponding tensor-product rule achieves, but with far fewer nodes.

## 3) Clenshaw-Curtis quadrature

The Smolyak construction requires a nested family of 1-D quadrature rules as
its basic sequence. `smolpack` uses **Clenshaw-Curtis** rules, which are based
on Chebyshev nodes.

### Standard Clenshaw-Curtis

The nodes for a Clenshaw-Curtis rule of order $m$ on $[0,1]$ are:

$$
x_j = \frac{1}{2}\left(1 - \cos\frac{j\pi}{m}\right), \quad j = 0, 1, \ldots, m.
$$

The standard nested sequence uses orders $m = 0, 2, 4, 8, 16, 32, \ldots$
(i.e., $2^{k-1}$ for $k \ge 2$), giving node counts of 1, 3, 5, 9, 17, 33,
65, …. The nested property — every node in a lower-order rule also appears
in higher-order rules — is essential for the Smolyak construction to be
efficient, since shared nodes need not be re-evaluated.

### Delayed Clenshaw-Curtis

Petras (2003) introduced a **delayed** variant that reorders the basic sequence
to reduce the total number of distinct nodes for a given level $q$. The key
insight is that the standard progression doubles the number of nodes at each
level, which is wasteful when the Smolyak combination technique only needs
a modest accuracy increase per stage.

The delayed sequence distributes node growth more gradually, achieving the
same polynomial exactness with fewer total function evaluations.

!!! note "Delayed vs. standard"
    The delayed variant (`int_smolyak`) and standard variant
    (`cc_int_smolyak`) produce different sparse grids with different
    node counts but achieve the same polynomial exactness at the same
    level $q$. The delayed variant typically uses fewer nodes.

## 4) Parameters and constraints

### The level parameter $q$

The level parameter `qq` (denoted $q$) controls the accuracy of the
integration. The **number of stages** is:

$$
k = q - d,
$$

where $d$ is the dimension. Higher $k$ yields a denser sparse grid with
more nodes and higher polynomial exactness.

The constraints are:

- $1 \le d < 40$ (dimension)
- $q - d < 48$ (number of stages)

### Weights

The weights $w_i$ in the cubature formula

$$
I[f] \approx \sum_{i=1}^{N} w_i\,f(\mathbf{x}_i)
$$

are computed via the Smolyak combination coefficients and the underlying
Clenshaw-Curtis quadrature weights. The `print_stats` flag reports both
the number of function evaluations $N$ and the number of weight evaluations
performed during grid construction.

## 5) Error and convergence

For sufficiently smooth integrands (functions with bounded mixed partial
derivatives), the Smolyak cubature error satisfies:

$$
|I[f] - A(q,d)[f]| = O\!\left(N^{-r}\,(\log N)^{(d-1)(r+1)}\right),
$$

where $r$ depends on the smoothness of $f$ and $N$ is the number of
nodes. This is near-optimal — much better than Monte Carlo's $O(N^{-1/2})$
rate, and only slightly worse than the optimal 1-D rate $O(N^{-r})$ due to
the logarithmic factors.

!!! tip "When Smolyak excels"
    Smolyak cubature is most effective for **smooth** integrands in
    **moderate** dimensions (roughly $d \le 20$). For very high dimensions
    or non-smooth integrands, Monte Carlo or quasi-Monte Carlo methods
    may be more appropriate.
