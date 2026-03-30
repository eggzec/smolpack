import numpy as np
import pytest

import smolpack


def _exp_sum(dim, x):
    """exp(x1 + x2 + ... + xd).  Exact = (e - 1)^d."""
    return np.exp(np.sum(x))


def _constant(dim, x) -> float:
    """Constant function 1.  Exact = 1 for any dimension."""
    return 1.0


def _sum_linear(dim, x):
    """x1 + x2 + ... + xd.  Exact = d / 2."""
    return np.sum(x)


def _product_x(dim, x):
    """x1 * x2 * ... * xd.  Exact = (1/2)^d."""
    return np.prod(x)


def _sum_squares(dim, x):
    """x1^2 + x2^2 + ... + xd^2.  Exact = d / 3."""
    return np.sum(x**2)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestIntSmolyakConstant:
    @pytest.mark.parametrize("dim", [1, 2, 5, 10])
    def test_constant(self, dim) -> None:
        result = smolpack.int_smolyak(_constant, dim=dim, qq=dim + 1)
        assert result == pytest.approx(1.0, abs=1e-12)


class TestIntSmolyakLinear:
    @pytest.mark.parametrize("dim", [1, 2, 3, 5])
    def test_linear(self, dim) -> None:
        exact = dim / 2.0
        result = smolpack.int_smolyak(_sum_linear, dim=dim, qq=dim + 2)
        assert result == pytest.approx(exact, abs=1e-12)


class TestIntSmolyakQuadratic:
    @pytest.mark.parametrize("dim", [1, 2, 3, 5])
    def test_sum_squares(self, dim) -> None:
        exact = dim / 3.0
        result = smolpack.int_smolyak(_sum_squares, dim=dim, qq=dim + 3)
        assert result == pytest.approx(exact, abs=1e-10)


class TestIntSmolyakExpSum:
    @pytest.mark.parametrize(
        "dim, qq, tol",
        [(1, 5, 1e-10), (2, 5, 1e-6), (3, 5, 1e-3), (3, 7, 1e-6), (5, 8, 1e-3)],
    )
    def test_exp_sum(self, dim, qq, tol) -> None:
        exact = (np.e - 1.0) ** dim
        result = smolpack.int_smolyak(_exp_sum, dim=dim, qq=qq)
        assert result == pytest.approx(exact, abs=tol)


class TestIntSmolyakProduct:
    @pytest.mark.parametrize("dim", [1, 2, 3, 4])
    def test_product(self, dim) -> None:
        exact = 0.5**dim
        result = smolpack.int_smolyak(_product_x, dim=dim, qq=dim + 3)
        assert result == pytest.approx(exact, abs=1e-10)


class TestIntSmolyakPrintStats:
    def test_print_stats_false(self) -> None:
        r1 = smolpack.int_smolyak(_exp_sum, dim=2, qq=4, print_stats=False)
        r2 = smolpack.int_smolyak(_exp_sum, dim=2, qq=4, print_stats=True)
        assert r1 == pytest.approx(r2, abs=1e-15)
