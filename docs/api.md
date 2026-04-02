# API Reference

`smolpack` provides a Python interface to the SMOLPACK C library (Petras, 2001; 2003)
for sparse-grid Smolyak cubature. The numerical core is implemented in C and compiled
via `f2py`, providing near-native performance with a NumPy-based Python interface.
See the [Theory](theory.md) and [Quickstart](quickstart.md) for mathematical background
and usage examples.

---

## Main Features

- Sparse-grid cubature over $[0,1]^d$ for any dimension $1 \le d < 40$
- Two solver variants: delayed and standard Clenshaw-Curtis Smolyak rules
- Simple integrand signature — any callable `f(dim, x) -> float` works
- Built-in function-evaluation counter via `get_count()`
- `print_stats` flag for runtime diagnostics

## Solvers

| Function | Underlying 1-D rule | Integrates |
|---|---|---|
| `int_smolyak` | Delayed Clenshaw-Curtis | $\int_{[0,1]^d} f(\mathbf{x})\,d\mathbf{x}$ |
| `cc_int_smolyak` | Standard Clenshaw-Curtis | $\int_{[0,1]^d} f(\mathbf{x})\,d\mathbf{x}$ |
| `get_count` | — | Returns the function-evaluation counter |

---

## Common interface

Both integration routines share the same calling convention and return the
approximate integral as a `float`.

### Arguments

| Argument | Type | Description |
|---|---|---|
| `f` | callable | Integrand function. Signature: `f(dim, x) -> float` |
| `dim` | `int` | Spatial dimension. Must satisfy $1 \le dim < 40$ |
| `qq` | `int` | Level parameter. Number of stages: $k = qq - dim$. Must satisfy $qq - dim < 48$ |
| `print_stats` | `bool` | If `True`, print function-call and weight-evaluation counts. Default: `False` |

!!! warning "Dimension limit"
    The maximum supported dimension is 39. Passing `dim >= 40` results in
    undefined behavior at the C level.

!!! warning "Level parameter"
    The constraint $qq - dim < 48$ must be satisfied. Violating this
    will produce incorrect results or crash.

### Callback signature

```python
def f(dim: int, x: numpy.ndarray) -> float: ...
```

| Parameter | Type | Description |
|---|---|---|
| `dim` | `int` | Spatial dimension |
| `x` | `numpy.ndarray` | 1-D array of shape `(dim,)` with each component in $[0,1]$ |
| **return** | `float` | Value of the integrand at point `x` |

!!! note "Lambda functions"
    Lambda functions work as integrands:
    `smolpack.int_smolyak(lambda d, x: np.exp(np.sum(x)), dim=3, qq=5)`

---

## `int_smolyak(f, dim, qq, print_stats=False)`

Approximate an integral over $[0,1]^{dim}$ using the **delayed
Clenshaw-Curtis** Smolyak algorithm.

The delayed variant uses a specially constructed basic sequence that typically
requires fewer function evaluations than the standard Clenshaw-Curtis sequence
for the same level of polynomial exactness.

$$
\int_{[0,1]^d} f(\mathbf{x})\,d\mathbf{x} \approx \sum_{i=1}^{N} w_i\,f(\mathbf{x}_i),
$$

where $\{(\mathbf{x}_i, w_i)\}$ are the sparse-grid nodes and weights
constructed from the delayed Clenshaw-Curtis basic sequence via the Smolyak
combination technique.

**Returns:** `float` — the approximate integral value.

```python
import numpy as np
import smolpack


def f(dim, x):
    return np.exp(np.sum(x))


result = smolpack.int_smolyak(f, dim=3, qq=5)
```

---

## `cc_int_smolyak(f, dim, qq, print_stats=False)`

Approximate an integral over $[0,1]^{dim}$ using the **standard
Clenshaw-Curtis** Smolyak algorithm.

This variant uses the classical Clenshaw-Curtis nested sequence
(1, 3, 5, 9, 17, 33, 65, … nodes) as the underlying 1-D rule. It may
require more function evaluations than `int_smolyak` at the same level but
can be more accurate for extremely smooth integrands.

**Returns:** `float` — the approximate integral value.

```python
import numpy as np
import smolpack


def f(dim, x):
    return np.exp(np.sum(x))


result = smolpack.cc_int_smolyak(f, dim=3, qq=5)
```

---

## `get_count()`

Return the value of the global function-evaluation counter from the
underlying C library.

The counter is reset to zero at the beginning of each integration call
(`int_smolyak` or `cc_int_smolyak`).

**Returns:** `int` — number of function evaluations recorded.

```python
import smolpack

result = smolpack.int_smolyak(f, dim=3, qq=5)
n = smolpack.get_count()
print(f"Evaluations: {n}")
```

---

## Example Workflows

`smolpack` supports any callable integrand. You can:

- Integrate smooth analytic functions in low to moderate dimensions
- Compare accuracy between `int_smolyak` and `cc_int_smolyak` at the same level
- Monitor the number of function evaluations with `get_count()` or `print_stats=True`
- Use lambda functions for quick interactive experiments

See the [Quickstart](quickstart.md) for full code examples, including:

- Exponential sum and product integrands
- High-dimensional integration (10-D)
- Genz oscillatory test functions
- Convergence study with increasing `qq`
- Solver comparison

All examples demonstrate both `int_smolyak` and `cc_int_smolyak`, along with
the `get_count()` and `print_stats` diagnostics.
