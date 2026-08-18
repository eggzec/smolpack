![SmolPack](assets/smolpack-banner.png)

# SmolPack

**Multidimensional Quadrature Using Sparse Grids for Python**

---

## Overview

`smolpack` is a Python library for efficient numerical integration (cubature)
over the unit hypercube $[0,1]^d$ using
[Smolyak's algorithm](https://en.wikipedia.org/wiki/Sparse_grid) with
Clenshaw-Curtis quadrature rules. The package allows you to **easily** and
**efficiently** integrate high-dimensional functions without the exponential
cost of classical tensor-product rules.

Numerical integration in multiple dimensions is a fundamental problem in
computational science. Standard tensor-product rules scale exponentially with
dimension — a 1-D rule with $n$ nodes becomes $n^d$ nodes in $d$ dimensions,
quickly making the computation intractable. Smolyak's algorithm constructs
**sparse grids** that achieve comparable polynomial exactness with dramatically
fewer nodes. `smolpack` approximates integrals of the form:

$$
I[f] = \int_{[0,1]^d} f(\mathbf{x})\,d\mathbf{x},
$$

using the Smolyak combination technique with Clenshaw-Curtis basic rules.

## Requirements

- [NumPy](http://www.numpy.org/)

## Example Usage

```python
import numpy as np
import smolpack


# Integrate exp(x1 + x2 + x3) over [0,1]^3
# Exact value: (e - 1)^3 ≈ 5.073214
def exp_sum(dim, x):
    return np.exp(np.sum(x))


result = smolpack.int_smolyak(exp_sum, dim=3, qq=5)
print(f"Result = {result:.6f}")

# Convergence with increasing level
exact = (np.e - 1.0) ** 3
for qq in range(4, 9):
    result = smolpack.int_smolyak(exp_sum, dim=3, qq=qq)
    print(
        f"qq={qq}  k={qq - 3}  result={result:.10f}  error={abs(result - exact):.2e}"
    )

# Compare both solvers
r1 = smolpack.int_smolyak(exp_sum, dim=3, qq=7)
r2 = smolpack.cc_int_smolyak(exp_sum, dim=3, qq=7)
print(f"Delayed CC:  {r1:.10f}")
print(f"Standard CC: {r2:.10f}")

# Query function-evaluation count
n = smolpack.get_count()
print(f"Evaluations: {n}")
```

## Main Features

1. **Sparse-grid cubature** over $[0,1]^d$ for dimensions $1 \le d < 40$.
2. Two solver variants: delayed and standard Clenshaw-Curtis Smolyak rules.
3. High-performance C core compiled via `f2py`, with a clean NumPy-based Python API.
4. Simple integrand signature — any callable `f(dim, x) -> float` works.
5. Built-in function-evaluation counter via `get_count()`.
6. `print_stats` flag for runtime diagnostics.

## See Also

- [NumPy](http://www.numpy.org/): Array library used for the integrand interface
- [SciPy](http://scipy.org): General scientific computing — complements `smolpack` for pre/post-processing
- [SMOLPACK C library](https://people.math.sc.edu/Burkardt/c_src/smolpack/smolpack.html): The original C implementation by Knut Petras, distributed by John Burkardt

## Acknowledgements

The author thanks [Knut Petras](https://people.math.sc.edu/Burkardt/c_src/smolpack/smolpack.html)
for the original SMOLPACK C library, which provides the numerical core of this
package and is distributed via [John Burkardt's page](https://people.math.sc.edu/Burkardt/c_src/smolpack/smolpack.html).

## References

- Smolyak, S., "Quadrature and Interpolation Formulas for Tensor Products of Certain Classes of Functions," *Doklady Akademii Nauk SSSR*, 4, 240–243, 1963.
- Petras, K., "Fast Calculation of Coefficients in the Smolyak Algorithm," *Numerical Algorithms*, 26(2), 93–109, 2001.
- Petras, K., "Smolyak Cubature of Given Polynomial Degree with Few Nodes for Increasing Dimension," *Numerische Mathematik*, 93(4), 729–753, 2003.
