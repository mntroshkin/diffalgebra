from fractions import Fraction
import diffalgebra as da

def test_integral_partial():
    R = da.ConstantPolyRing(constants=["t"])
    t = R.gen("t")
    f = t ** 2
    assert f.integral(t) == Fraction(1, 3) * t ** 3

def test_integral_total():
    A = da.DifferentialRing(functions=["u"])
    u = A.gen("u")
    assert u[1].integral() == u
    assert (2 * u * u[1]).integral() == u ** 2

def test_integral_does_not_exist():
    A = da.DifferentialRing(functions=["u", "v"])
    u, v = A.gens()
    assert u.integral() is None
    assert (u * u[1] ** 2).integral() is None
    assert (u * u[1] ** 2 + u[4]).integral() is None
    assert (u[1] + v).integral() is None