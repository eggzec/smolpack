import numpy as np
import pytest

import smolpack


def _exp_sum(dim, x):
    """exp(x1 + x2 + ... + xd)."""
    return np.exp(np.sum(x))


def _constant(dim, x) -> float:
    """Constant function 1."""
    return 1.0


class TestGetCount:
    def test_count_is_int(self) -> None:
        smolpack.int_smolyak(_exp_sum, dim=2, qq=4)
        assert isinstance(smolpack.get_count(), int)

    def test_count_nonnegative(self) -> None:
        smolpack.int_smolyak(_exp_sum, dim=2, qq=4)
        assert smolpack.get_count() >= 0

    def test_count_after_cc(self) -> None:
        smolpack.cc_int_smolyak(_exp_sum, dim=2, qq=4)
        assert smolpack.get_count() >= 0


class TestCrossAlgorithm:
    @pytest.mark.parametrize("dim", [1, 2, 3])
    def test_convergence_agreement(self, dim) -> None:
        qq = dim + 5
        r1 = smolpack.int_smolyak(_exp_sum, dim=dim, qq=qq)
        r2 = smolpack.cc_int_smolyak(_exp_sum, dim=dim, qq=qq)
        exact = (np.e - 1.0) ** dim
        # Both should be close to exact
        assert r1 == pytest.approx(exact, abs=1e-4)
        assert r2 == pytest.approx(exact, abs=1e-4)

    def test_1d_match(self) -> None:
        """In 1-D both algorithms should give very close results."""
        r1 = smolpack.int_smolyak(_exp_sum, dim=1, qq=6)
        r2 = smolpack.cc_int_smolyak(_exp_sum, dim=1, qq=6)
        assert r1 == pytest.approx(r2, rel=1e-8)


class TestEdgeCases:
    def test_dim_1(self) -> None:
        result = smolpack.int_smolyak(_constant, dim=1, qq=2)
        assert result == pytest.approx(1.0, abs=1e-14)

    def test_high_dim_constant(self) -> None:
        result = smolpack.int_smolyak(_constant, dim=20, qq=21)
        assert result == pytest.approx(1.0, abs=1e-10)
