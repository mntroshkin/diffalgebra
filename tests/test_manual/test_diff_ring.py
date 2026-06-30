import diffalgebra as da

def test_zero_representation():
    A = da.DifferentialRing(functions=["u"], ring_name="A")
    assert str(A.promote(0)) == "0"

def test_derivative_representation():
    A = da.DifferentialRing(functions=["u"], ring_name="A")
    u = A.gen("u")
    assert str(u[0]) == "u"
    assert str(u[1]) == "u_x"
    assert str(u[2]) == "u_xx"
    assert str(u[3]) == "u_xxx"
    assert str(u[5]) == "u_5"

def test_product_rule():
    A = da.DifferentialRing(functions=["u", "v"], ring_name="A")
    u, v = A.gens()
    assert (u * v).diff() == u[1] * v + u * v[1]

def test_partial_derivative_zero():
    A = da.DifferentialRing(functions=["u", "v"], ring_name="A")
    u, v = A.gens()
    f = u * u[2] + u[1] ** 2
    assert f.d(v) == 0

def test_partial_derivative():
    A = da.DifferentialRing(functions=["u"], ring_name="A")
    u = A.gen("u")
    assert (u * u[1] ** 2 * u[2]).d(u[1]) == 2 * u * u[1] * u[2]

def test_coefficient():
    R = da.ConstantPolyRing(constants=["t"], ring_name="R")
    t = R.gen("t")
    A = da.DifferentialRing(functions=["u"], base_ring=R, ring_name="A")
    u = A.gen("u")
    f = (1 + t) + 2 * u + t ** 2 * u * u[1]
    assert f.coefficient(1) == 1 + t
    assert f.coefficient(u) == 2
    assert f.coefficient(u[1]) == 0
    assert f.coefficient(u * u[1]) == t ** 2  