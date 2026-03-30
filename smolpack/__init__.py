"""SMOLPACK - Multidimensional Quadrature Using Smolyak Sparse Grids.

Efficient numerical integration (cubature) over the unit hypercube
[0,1]^d using Smolyak's algorithm with Clenshaw-Curtis quadrature rules.

Two algorithms are provided:

* :func:`int_smolyak` — delayed Clenshaw-Curtis (fewer evaluations).
* :func:`cc_int_smolyak` — standard Clenshaw-Curtis (classical rule).

Example
-------
>>> import numpy as np
>>> import smolpack
>>>
>>> def my_func(dim, x):
...     return np.exp(np.sum(x))
>>> result = smolpack.int_smolyak(my_func, dim=3, qq=5)
"""

from ._smolpack import delayed_cc as _delayed_cc
from ._smolpack import get_count
from ._smolpack import standard_cc as _standard_cc


def int_smolyak(
    f: callable, dim: int, qq: int, *, print_stats: bool = False
) -> float:
    """Approximate an integral over [0,1]^dim (delayed Clenshaw-Curtis).

    Parameters
    ----------
    f : callable
        Integrand ``f(dim, x) -> float``, where *x* is a 1-D NumPy array.
    dim : int
        Spatial dimension (1 ≤ dim < 40).
    qq : int
        Level parameter.  The number of stages is ``k = qq - dim``.
    print_stats : bool, optional
        If *True*, print function-call and weight-evaluation counts.

    Returns
    -------
    float
        Approximated integral value.

    See Also
    --------
    cc_int_smolyak : Standard Clenshaw-Curtis variant.
    """
    return _delayed_cc(dim, qq, int(print_stats), f)


def cc_int_smolyak(
    f: callable, dim: int, qq: int, *, print_stats: bool = False
) -> float:
    """Approximate an integral over [0,1]^dim (standard Clenshaw-Curtis).

    Parameters
    ----------
    f : callable
        Integrand ``f(dim, x) -> float``, where *x* is a 1-D NumPy array.
    dim : int
        Spatial dimension (1 ≤ dim < 40).
    qq : int
        Level parameter.  The number of stages is ``k = qq - dim``.
    print_stats : bool, optional
        If *True*, print function-call and weight-evaluation counts.

    Returns
    -------
    float
        Approximated integral value.

    See Also
    --------
    int_smolyak : Delayed Clenshaw-Curtis variant (fewer evaluations).
    """
    return _standard_cc(dim, qq, int(print_stats), f)


__all__ = ["cc_int_smolyak", "get_count", "int_smolyak"]
