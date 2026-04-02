# Quickstart

This page walks through progressively richer examples of using `smolpack`.
Each example is self-contained and can be pasted directly into a Python script
or interactive session. For mathematical background, see the [Theory](theory.md) section.

---

## Exponential sum (3-D)

The simplest non-trivial test — integrate $e^{x_1+x_2+x_3}$ over the unit cube:

$$
I = \int_{[0,1]^3} e^{x_1+x_2+x_3}\,dx = (e-1)^3 \approx 5.073214.
$$

```python
import numpy as np
import smolpack


def exp_sum(dim, x):
    return np.exp(np.sum(x))


result = smolpack.int_smolyak(exp_sum, dim=3, qq=5)
print(f"Result = {result:.6f}")
```

## Constant function (sanity check)

A constant integrand should return exactly 1 for any dimension:

$$
\int_{[0,1]^d} 1\,dx = 1.
$$

```python
import smolpack

result = smolpack.int_smolyak(lambda d, x: 1.0, dim=10, qq=11)
print(f"Result = {result:.12f}")  # 1.000000000000
```

## Comparing the two solvers

Both solvers applied to the same integrand in 3-D:

```python
import numpy as np
import smolpack


def exp_sum(dim, x):
    return np.exp(np.sum(x))


exact = (np.e - 1.0) ** 3

r1 = smolpack.int_smolyak(exp_sum, dim=3, qq=7)
r2 = smolpack.cc_int_smolyak(exp_sum, dim=3, qq=7)

print(f"Delayed CC:  {r1:.10f}  (error {abs(r1 - exact):.2e})")
print(f"Standard CC: {r2:.10f}  (error {abs(r2 - exact):.2e})")
print(f"Exact:       {exact:.10f}")
```

## Convergence with increasing `qq`

The level parameter `qq` controls accuracy. The number of Smolyak stages is
$k = qq - dim$; higher $k$ yields more accurate results at the cost of
additional function evaluations.

```python
import numpy as np
import smolpack


def exp_sum(dim, x):
    return np.exp(np.sum(x))


exact = (np.e - 1.0) ** 3

for qq in range(4, 10):
    result = smolpack.int_smolyak(exp_sum, dim=3, qq=qq)
    print(
        f"qq={qq}  k={qq - 3}  result={result:.10f}  error={abs(result - exact):.2e}"
    )
```

## Product integrand (separable)

For the product $f(\mathbf{x}) = x_1 x_2 \cdots x_d$, the exact value is $(1/2)^d$:

```python
import numpy as np
import smolpack


def product_x(dim, x):
    return np.prod(x)


dim = 5
result = smolpack.int_smolyak(product_x, dim=dim, qq=dim + 3)
exact = 0.5**dim
print(f"Result = {result:.10f}  (exact = {exact})")
```

## High-dimensional integration

Smolyak's strength is moderate- to high-dimensional problems. Here we integrate
the product of cosines in 10 dimensions:

$$
\int_{[0,1]^{10}} \prod_{i=1}^{10} \cos(x_i)\,dx = \sin(1)^{10} \approx 1.5867 \times 10^{-4}.
$$

```python
import numpy as np
import smolpack


def product_cosine(dim, x):
    return np.prod(np.cos(x))


result = smolpack.int_smolyak(product_cosine, dim=10, qq=12)
exact = np.sin(1.0) ** 10
print(f"Result = {result:.10e}")
print(f"Exact  = {exact:.10e}")
```

## Genz oscillatory test function

The [Genz test functions](https://www.math.wsu.edu/faculty/genz/software/software.html)
are standard benchmarks for cubature methods. The oscillatory family is:

$$
f(\mathbf{x}) = \cos\!\left(2\pi u + \sum_{i=1}^d a_i x_i\right).
$$

```python
import numpy as np
import smolpack


def genz_oscillatory(dim, x):
    a = np.ones(dim)
    return np.cos(2.0 * np.pi * a[0] + np.sum(a * x))


result = smolpack.int_smolyak(genz_oscillatory, dim=5, qq=8)
print(f"Oscillatory integral = {result:.8f}")
```

## Printing solver statistics

Set `print_stats=True` to see function-call and weight-evaluation counts:

```python
import numpy as np
import smolpack


def exp_sum(dim, x):
    return np.exp(np.sum(x))


result = smolpack.int_smolyak(exp_sum, dim=3, qq=5, print_stats=True)
```

The counter can also be queried programmatically after any integration call:

```python
n = smolpack.get_count()
print(f"Function evaluations: {n}")
```

## Tips

!!! tip "Choosing the level parameter"
    Start with $qq = dim + 2$ and increase until the result stabilizes.
    Each additional level adds more sparse-grid nodes, improving accuracy
    at the cost of more integrand evaluations.

!!! tip "When to use `cc_int_smolyak`"
    Prefer the standard Clenshaw-Curtis variant when reproducibility against
    published sparse-grid node tables is important, or when the integrand is
    extremely smooth and benefits from the standard node hierarchy.

!!! tip "Lambda functions"
    Lambda functions work as integrands for quick experiments:

    ```python
    smolpack.int_smolyak(lambda d, x: sum(x), dim=3, qq=5)
    ```
