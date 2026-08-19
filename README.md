![SmolPack](https://raw.githubusercontent.com/eggzec/smolpack/main/docs/assets/smolpack-banner.png)

# SmolPack

**Multidimensional Quadrature Using Sparse Grids for Python**

[![Tests](https://github.com/eggzec/smolpack/actions/workflows/test.yml/badge.svg)](https://github.com/eggzec/smolpack/actions/workflows/test.yml)
[![Documentation](https://github.com/eggzec/smolpack/actions/workflows/docs.yml/badge.svg)](https://github.com/eggzec/smolpack/actions/workflows/docs.yml)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

[![codecov](https://codecov.io/github/eggzec/smolpack/graph/badge.svg)](https://codecov.io/github/eggzec/smolpack)
[![License: LGPL-2.1](https://img.shields.io/badge/License-LGPL--2.1-blue.svg)](LICENSE)

[![PyPI Downloads](https://img.shields.io/pypi/dm/smolpack.svg?label=PyPI%20downloads)](https://pypi.org/project/smolpack/)
[![Python versions](https://img.shields.io/pypi/pyversions/smolpack.svg)](https://pypi.org/project/smolpack/)

`smolpack` is a high-performance Python library for numerical integration
(cubature) over the unit hypercube [0,1]^d using
[Smolyak's algorithm](https://en.wikipedia.org/wiki/Sparse_grid) with
Clenshaw-Curtis quadrature rules. It approximates integrals of the form:

$$I[f] = \int_{[0,1]^d} f(\mathbf{x})\,d\mathbf{x}$$

Two solvers are provided: a **delayed Clenshaw-Curtis** variant (fewer
evaluations) and a **standard Clenshaw-Curtis** variant (classical rule).

## Quick example

```python
import numpy as np
import smolpack


def my_func(dim, x):
    return np.exp(np.sum(x))


result = smolpack.int_smolyak(my_func, dim=3, qq=5)
```

## Installation

```bash
pip install smolpack
```

Requires Python 3.10+ and NumPy. See the
[full installation guide](https://eggzec.github.io/smolpack/installation/) for
uv, poetry, and source builds.

## Documentation

- [Theory](https://eggzec.github.io/smolpack/theory/) — sparse grids, Smolyak construction, Clenshaw-Curtis rules
- [Quickstart](https://eggzec.github.io/smolpack/quickstart/) — runnable examples
- [API Reference](https://eggzec.github.io/smolpack/api/) — solver signatures and parameters
- [References](https://eggzec.github.io/smolpack/references/) — literature citations

## Attribution

The original C library SMOLPACK was written by Knut Petras and is distributed
via [John Burkardt's page](https://people.math.sc.edu/Burkardt/c_src/smolpack/smolpack.html).

## License

LGPL-2.1 — see [LICENSE](LICENSE).
