from hypothesis import given, strategies as st, settings
from .strategies import diff_polynomial, polynomial
import diffalgebra as da

R = da.ConstantPolyRing(constants=["t", "s"])

polynomial_small = polynomial(ring=R, max_terms=3)

B = da.DifferentialRing(functions=["w"], base_ring=R)
w = B.gen("w")
diff_polynomial_tiny = diff_polynomial(ring=B, max_terms=2, max_nonlinearity=2)


@given(polynomial_small, diff_polynomial_tiny)
def test_evolution_of_constant_is_zero(f: da.ConstantPolynomial, image: da.DifferentialPolynomial):
    F = B.promote(f)
    D = da.EvolutionOperator(ring=B, mapping={w: image})
    assert D(F) == 0
    

@settings(deadline=None)
@given(diff_polynomial_tiny, diff_polynomial_tiny, diff_polynomial_tiny)
def test_product_rule_for_evolution(F: da.DifferentialPolynomial, G: da.DifferentialPolynomial,
                                     image: da.DifferentialPolynomial):
    D = da.EvolutionOperator(ring=B, mapping={w: image})
    assert D(F * G) == D(F) * G + F * D(G)


@settings(deadline=None)
@given(diff_polynomial_tiny, diff_polynomial_tiny,
       st.integers(min_value=1, max_value=5))
def test_evolution_commutes_with_differential(F: da.DifferentialPolynomial,
                                              image: da.DifferentialPolynomial,
                                              order: int):
    D = da.EvolutionOperator(ring=B, mapping={w: image})
    assert D(F.diff(order)) == D(F).diff(order)