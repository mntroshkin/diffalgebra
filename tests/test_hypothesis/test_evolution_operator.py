from hypothesis import given, strategies as st, settings
from .strategies import diff_polynomial, polynomial
import diffalgebra as da

R = da.ConstantPolyRing(constants=["a", "b", "c", "d", "e"])

polynomial_small = polynomial(ring=R, max_terms=3)

A = da.DifferentialRing(functions=["u", "v"], base_ring=R)
diff_polynomial_small = diff_polynomial(ring=A, max_terms=5, max_nonlinearity=2)
diff_polynomial_medium = diff_polynomial(ring=A, max_terms=10, max_nonlinearity=3)

@given(polynomial_small, st.tuples(diff_polynomial_small, diff_polynomial_small))
def test_derivation_of_constant_is_zero(f: da.ConstantPolynomial,
                                        images: tuple[da.DifferentialPolynomial, ...]):
    F = A.promote(f)
    D = da.EvolutionOperator(ring=A, mapping={var: image for var, image in zip(A.gens(), images)})
    assert D(F) == 0