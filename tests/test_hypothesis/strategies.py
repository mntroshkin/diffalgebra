from hypothesis import strategies as st
import diffalgebra as da


@st.composite
def polynomial(draw, ring: da.ConstantPolyRing, max_terms: int = 5) -> da.ConstantPolynomial:
    gens = ring.gens()
    gen_count = len(gens)

    monomial_list = draw(st.lists(st.lists(st.integers(min_value=0, max_value=10), 
                                    min_size=gen_count, max_size=gen_count),
                        min_size=0, max_size=max_terms))
    term_count = len(monomial_list)
    monomial_list = map(tuple, monomial_list)
    coef_list = draw(st.lists(st.integers(min_value=-5, max_value=5),
                        min_size=term_count, max_size=term_count))
    terms = list(zip(monomial_list, coef_list))
    return da.ConstantPolynomial(ring, terms)


R = da.ConstantPolyRing(constants=["a", "b", "c", "d", "e"])

polynomial_small = polynomial(ring=R, max_terms=3)
polynomial_medium = polynomial(ring=R, max_terms=10)


@st.composite
def diff_polynomial(draw, ring: da.DifferentialRing, 
                    max_terms: int = 5, max_nonlinearity: int = 3) -> da.DifferentialPolynomial:
    gens = ring.gens()
    gen_count = len(gens)

    monomial_list = draw(st.lists(st.lists(st.lists(st.tuples(st.integers(min_value=0, max_value=3),
                                                                 st.integers(min_value=0, max_value=3)),
                                                       min_size=0, max_size=max_nonlinearity),
                                              min_size=gen_count, max_size=gen_count), 
                                     min_size=0, max_size=max_terms))
    term_count = len(monomial_list)
    monomial_list = map(tuple, monomial_list)
    coef_list = draw(st.lists(st.integers(min_value=-5, max_value=5),
                        min_size=term_count, max_size=term_count))
    terms = list(zip(monomial_list, coef_list))
    return da.DifferentialPolynomial(ring, terms)


A = da.DifferentialRing(functions=["u", "v"])
diff_polynomial_small = diff_polynomial(ring=A, max_terms=5, max_nonlinearity=2)
diff_polynomial_medium = diff_polynomial(ring=A, max_terms=10, max_nonlinearity=3)
