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
