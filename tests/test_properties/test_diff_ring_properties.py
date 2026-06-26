from hypothesis import given, strategies as st, settings
import diffalgebra as da

A = da.DifferentialRing(functions=["u"], ring_name="A")
u = A.gen("u")

@st.composite
def diff_polynomial(draw, max_terms: int) -> da.DifferentialPolynomial:
    coefs_and_exps = draw(st.lists(st.tuples(st.integers(min_value=-10, max_value=10), 
                                             st.integers(min_value=0, max_value=10)),
                                   max_size=max_terms))
    f = A.promote(sum(coef * (u[i]) ** exp for i, (coef, exp) in enumerate(coefs_and_exps)))
    return f


@given(diff_polynomial(max_terms=20))
def test_multiply_diff_polynomial_by_zero(f: da.DifferentialPolynomial):
    assert (f * 0) == 0


@given(diff_polynomial(max_terms=20))
def test_zeroth_power_is_one(f: da.DifferentialPolynomial):
    assert (f ** 0) == 1

@settings(deadline=None)
@given(diff_polynomial(max_terms=5),
       st.integers(min_value=0, max_value=5),
       st.integers(min_value=0, max_value=5))
def test_diff_polynomial_exponents_add_under_product(f: da.DifferentialPolynomial, a: int, b: int):
    assert (f ** a) * (f ** b) == f ** (a + b)


@given(diff_polynomial(max_terms=20))
def test_zeroth_derivative_is_original(f: da.DifferentialPolynomial):
    assert f.diff(order=0) == f


@given(diff_polynomial(max_terms=5),
       st.integers(min_value=0, max_value=5),
       st.integers(min_value=0, max_value=5))
def test_derivative_composition(f: da.DifferentialPolynomial, a: int, b: int):
    assert f.diff(order=a).diff(order=b) == f.diff(order=a+b)


@settings(deadline=None)
@given(diff_polynomial(max_terms=10), diff_polynomial(max_terms=10))
def test_product_rule(f: da.DifferentialPolynomial, g: da.DifferentialPolynomial):
    assert (f * g).diff() == f * g.diff() + f.diff() * g