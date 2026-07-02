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


@st.composite
def diff_polynomial(draw, ring: da.DifferentialRing, 
                    max_terms: int = 5, max_nonlinearity: int = 2) -> da.DifferentialPolynomial:
    gens = ring.gens()
    gen_count = len(gens)

    terms = draw(st.lists(st.tuples(st.lists(st.lists(st.integers(min_value=1, max_value=3),
                                                    min_size=0, 
                                                    max_size=max_nonlinearity).map(enumerate).map(list),
                                            min_size=gen_count, max_size=gen_count).map(tuple),
                                    st.integers(min_value=-5, max_value=5)),
                                    min_size=0, max_size=max_terms))
    return da.DifferentialPolynomial(ring, terms)
