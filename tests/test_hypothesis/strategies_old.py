from hypothesis import strategies as st
import diffalgebra as da


@st.composite
def polynomial(draw, ring: da.ConstantPolyRingOld, max_terms: int = 5) -> da.ConstantPolynomialOld:
    gens = ring.gens()
    gen_count = len(gens)

    exponents_list = draw(st.lists(st.lists(st.integers(min_value=0, max_value=10), 
                                    min_size=gen_count, max_size=gen_count),
                        min_size=0, max_size=max_terms))
    terms = len(exponents_list)
    coef_list = draw(st.lists(st.integers(min_value=-5, max_value=5),
                        min_size=terms, max_size=terms))
    monomials = []
    for coef, exponents in zip(coef_list, exponents_list):
        monomial = coef
        for exponent, var in zip(exponents, gens):
            monomial *= var ** exponent
        monomials.append(monomial)
    return ring.promote(sum(monomials))


R = da.ConstantPolyRingOld(constants=["a", "b", "c", "d", "e"])

polynomial_old_small = polynomial(ring=R, max_terms=3)
polynomial_old_medium = polynomial(ring=R, max_terms=10)

@st.composite
def diff_polynomial(draw, ring: da.DifferentialRing, 
                    max_terms: int = 5, max_nonlinearity: int = 3) -> da.DifferentialPolynomial:
    gens = ring.gens()
    gen_count = len(gens)

    exp_degrees_list = draw(st.lists(st.lists(st.lists(st.tuples(st.integers(min_value=0, max_value=3),
                                                                 st.integers(min_value=0, max_value=3)),
                                                       min_size=0, max_size=max_nonlinearity),
                                              min_size=gen_count, max_size=gen_count), 
                                     min_size=0, max_size=max_terms))
    terms = len(exp_degrees_list)
    coef_list = draw(st.lists(st.integers(min_value=-10, max_value=10),
                              min_size=terms, max_size=terms))
    monomials = []
    for coef, exp_degrees_term in zip(coef_list, exp_degrees_list):
        monomial = coef
        for exp_degrees_var, var in zip(exp_degrees_term, gens):
            for exp, degree in exp_degrees_var:
                monomial *= var[degree] ** exp
        monomials.append(monomial)
    return ring.promote(sum(monomials))


A = da.DifferentialRing(functions=["u", "v"])
diff_polynomial_small = diff_polynomial(ring=A, max_terms=5, max_nonlinearity=2)
diff_polynomial_medium = diff_polynomial(ring=A, max_terms=10, max_nonlinearity=3)
