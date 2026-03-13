# smolpack

**Sparse-grid Smolyak cubature over $[0,1]^d$ for Python.**

---

## What is smolpack?

`smolpack` is a Python library for efficient numerical integration (cubature)
over the unit hypercube $[0,1]^d$ using
[Smolyak's algorithm](https://en.wikipedia.org/wiki/Sparse_grid) with
Clenshaw-Curtis quadrature rules. The numerical core is written in C and
compiled via `f2py`, giving near-native performance while providing a clean,
NumPy-based Python API.

**Numerical integration** (quadrature in 1-D, cubature in higher dimensions) is
a fundamental problem in computational science. Standard tensor-product rules
scale exponentially with dimension — a 1-D rule with $n$ nodes becomes $n^d$
nodes in $d$ dimensions, quickly making the computation intractable. Smolyak's
algorithm circumvents this *curse of dimensionality* by constructing **sparse
grids** that achieve comparable polynomial exactness with dramatically fewer
nodes.

`smolpack` approximates integrals of the form:

$$
I[f] = \int_{[0,1]^d} f(\mathbf{x})\,d\mathbf{x},
$$

using the Smolyak combination technique with Clenshaw-Curtis basic rules.

## Available solvers

| Solver | Underlying 1-D rule | Characteristics |
|---|---|---|
| `int_smolyak` | Delayed Clenshaw-Curtis | Fewer function evaluations for a given accuracy level |
| `cc_int_smolyak` | Standard Clenshaw-Curtis | Classical nested rule (1, 3, 5, 9, 17, 33, 65, … nodes) |

!!! tip "Choosing a solver"
    For most applications, start with `int_smolyak` — the delayed variant
    achieves the same polynomial exactness with fewer integrand evaluations.
    Use `cc_int_smolyak` when compatibility with the classical Clenshaw-Curtis
    node hierarchy is required.

## Quick example

```python
import numpy as np
import smolpack

# Integrate exp(x1 + x2 + x3) over [0,1]^3
# Exact value: (e - 1)^3 ≈ 5.073214

def my_func(dim, x):
    return np.exp(np.sum(x))

result = smolpack.int_smolyak(my_func, dim=3, qq=5)
print(f"Integral ≈ {result:.6f}")
```

## License

smolpack is distributed under the
[LGPL-2.1](https://github.com/eggzec/smolpack/blob/main/LICENSE) license.
