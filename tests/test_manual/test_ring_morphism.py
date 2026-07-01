import diffalgebra as da

def test_product_under_ring_morphism():
    R = da.ConstantPolyRing(constants=["r", "s"])
    S = da.ConstantPolyRing(constants=["t"])
    r = R.gen("r")
    s = R.gen("s")
    t = S.gen("t")
    phi = da.RingMorphism(source=R, target=S, mapping={r: 1 + t, s: 2 + t})
    assert phi(r * s) == (1 + t) * (2 + t)

def test_powers_under_ring_morphism():
    R = da.ConstantPolyRing(constants=["t"])
    S = da.ConstantPolyRing(constants=["s"])
    t = R.gen("t")
    s = S.gen("s")
    phi = da.RingMorphism(source=R, target=S, mapping={t: 1 + s})
    assert phi(t ** 20) == (1 + s) ** 20

def test_derivative_under_diff_ring_morphism():
    A = da.DifferentialRing(functions=["u"])
    B = da.DifferentialRing(functions=["v"])
    u = A.gen("u")
    v = B.gen("v")
    F = da.DiffRingMorphism(source=A, target=B, mapping={u: v ** 2})
    assert F(u[1]) == 2 * v * v[1]

def test_diff_morphism_implicit_base_identity():
    R = da.ConstantPolyRing(constants=["t"])
    t = R.gen("t")
    A = da.DifferentialRing(functions=["u"], base_ring=R)
    B = da.DifferentialRing(functions=["v"], base_ring=R)
    u = A.gen("u")
    v = B.gen("v")
    F = da.DiffRingMorphism(source=A, target=B, mapping={u: v})
    assert F(t) == t

def test_diff_morphism_implicit_base_from_QQ():
    R = da.ConstantPolyRing(constants=["t"])
    A = da.DifferentialRing(functions=["u"])
    B = da.DifferentialRing(functions=["v"], base_ring=R)
    u = A.gen("u")
    v = B.gen("v")
    F = da.DiffRingMorphism(source=A, target=B, mapping={u: v})
    assert F(1) == 1