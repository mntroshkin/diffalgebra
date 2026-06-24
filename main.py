from src.diffalgebra import ConstantPolyRing

A = ConstantPolyRing(constants=["t", "s"], ring_name="A")
t = A.get_constant("t")
s = A.get_constant("s")

print((t + s) * (t - s))