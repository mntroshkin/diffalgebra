from fractions import Fraction
import math
import diffalgebra as da

def test_antiderivative_constant():
    R = da.ConstantPolyRing(constants=["t"])
    t = R.gen("t")
    f = t ** 2
    assert f.int(t) == Fraction(1, 3) * t ** 3