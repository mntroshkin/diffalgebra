import diffalgebra as da

def test_square_of_binomial():
    R = da.ConstantPolyRing(constants=["t"], ring_name="R")
    t = R.gen("t")
    assert (1 + t) ** 2 == t ** 2 + 2 * t + 1

def test_square_of_binomial_representation():
    R = da.ConstantPolyRing(constants=["t"], ring_name="R")
    t = R.gen("t")
    assert str((1 + t) ** 2) == "t^2+2*t+1"

def test_negative_coefs_representation():
    R = da.ConstantPolyRing(constants=["t"], ring_name="R")
    t = R.gen("t")
    assert(str(-t ** 2 - 2 * t - 3) == "-t^2-2*t-3")

def test_zero_representation():
    R = da.ConstantPolyRing(constants=["t"], ring_name="R")
    assert str(R.promote(0)) == "0"