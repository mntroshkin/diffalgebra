from hypothesis import given, strategies as st, settings
from .strategies import polynomial_small, polynomial_medium
import pytest

import diffalgebra as da


@given(polynomial_medium)
def test_adding_one_inequality(f: da.ConstantPolynomial):
    assert f + 1 != f


@given(polynomial_medium)
def test_multiply_polynomial_by_zero(f: da.ConstantPolynomial):
    assert (f * 0) == 0


@given(polynomial_medium)
def test_zeroth_power_is_one(f: da.ConstantPolynomial):
    assert (f ** 0) == 1


@settings(deadline=None)
@pytest.mark.slow
@given(polynomial_small,
       st.integers(min_value=0, max_value=5),
       st.integers(min_value=0, max_value=5))
def test_polynomial_exponents_add_under_product(f: da.ConstantPolynomial, a: int, b: int):
    assert (f ** a) * (f ** b) == f ** (a + b)


@settings(deadline=None)
@pytest.mark.slow
@given(polynomial_small,
       st.integers(min_value=1, max_value=4),
       st.integers(min_value=1, max_value=4))
def test_polynomial_exponents_multiply_under_powers(f: da.ConstantPolynomial, a: int, b: int):
    assert (f ** a) ** b == f ** (a * b)