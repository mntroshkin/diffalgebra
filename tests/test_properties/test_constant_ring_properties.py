from hypothesis import given, strategies as st, settings
import diffalgebra as da


R = da.ConstantPolyRing(constants=["t"], ring_name="R")
t = R.gen("t")


@st.composite
def polynomial(draw, max_terms: int) -> da.ConstantPolynomial:
    coefs = draw(st.lists(st.integers(min_value=-10, max_value=10), max_size=max_terms))
    f = R.promote(sum(coef * (t ** i) for i, coef in enumerate(coefs)))
    return f


@given(polynomial(max_terms=20))
def test_adding_one_inequality(f: da.ConstantPolynomial):
    assert f + 1 != f


@given(polynomial(max_terms=20))
def test_multiply_polynomial_by_zero(f: da.ConstantPolynomial):
    assert (f * 0) == 0


@given(polynomial(max_terms=20))
def test_zeroth_power_is_one(f: da.ConstantPolynomial):
    assert (f ** 0) == 1


@given(polynomial(max_terms=10),
       st.integers(min_value=0, max_value=5),
       st.integers(min_value=0, max_value=5))
def test_polynomial_exponents_add_under_product(f: da.ConstantPolynomial, a: int, b: int):
    assert (f ** a) * (f ** b) == f ** (a + b)


@given(polynomial(max_terms=10),
       st.integers(min_value=0, max_value=5),
       st.integers(min_value=0, max_value=5))
def test_polynomial_exponents_multiply_under_powers(f: da.ConstantPolynomial, a: int, b: int):
    assert (f ** a) ** b == f ** (a * b)