from hypothesis import given, strategies as st, settings
from .strategies import diff_polynomial
import diffalgebra as da


A = da.DifferentialRing(functions=["u", "v"])
diff_polynomial_small = diff_polynomial(ring=A, max_terms=5, max_nonlinearity=2)
diff_polynomial_medium = diff_polynomial(ring=A, max_terms=10, max_nonlinearity=3)


@given(diff_polynomial_medium)
def test_multiply_diff_polynomial_by_zero(f: da.DifferentialPolynomial):
    assert (f * 0) == 0


@given(diff_polynomial_medium)
def test_zeroth_power_is_one(f: da.DifferentialPolynomial):
    assert (f ** 0) == 1


@settings(deadline=None)
@given(diff_polynomial_small,
       st.integers(min_value=0, max_value=5),
       st.integers(min_value=0, max_value=5))
def test_diff_polynomial_exponents_add_under_product(f: da.DifferentialPolynomial, a: int, b: int):
    assert (f ** a) * (f ** b) == f ** (a + b)


@given(diff_polynomial_medium)
def test_zeroth_derivative_is_original(f: da.DifferentialPolynomial):
    assert f.diff(order=0) == f

@settings(deadline=None)
@given(diff_polynomial_medium,
       st.integers(min_value=1, max_value=3),
       st.integers(min_value=1, max_value=3))
def test_derivative_composition(f: da.DifferentialPolynomial, a: int, b: int):
    assert f.diff(order=a).diff(order=b) == f.diff(order=a+b)


@settings(deadline=None)
@given(diff_polynomial_small, diff_polynomial_small)
def test_product_rule(f: da.DifferentialPolynomial, g: da.DifferentialPolynomial):
    assert (f * g).diff() == f * g.diff() + f.diff() * g


@given(diff_polynomial_medium)
def test_trivial_equality(f: da.DifferentialPolynomial):
    assert f == f

@given(diff_polynomial_medium)
def test_variational_derivative_of_differential(f: da.DifferentialPolynomial):
    u, v = A.gens()
    assert f.diff().delta(u) == 0
    assert f.diff().delta(v) == 0

@given(diff_polynomial_medium)
def test_antiderivative_of_differential(f: da.DifferentialPolynomial):
    f0 = f - f.coefficient(1)
    assert f0.diff().integral() == f0