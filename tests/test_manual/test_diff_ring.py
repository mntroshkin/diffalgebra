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
    u = A.gen("u")
    v = A.gen("v")
    assert (u * v).diff() == u[1] * v + u * v[1]