from .diff_ring import FuncGenerator, DifferentialPolynomial

def KdV_flow(u: FuncGenerator, order: int = 1) -> DifferentialPolynomial:
    if order < 0:
        raise ValueError
    K = u[1]
    for i in range(order):
        L = K.integral()
        K = L.diff(3) + 4 * u * L.diff() + 2 * u[1] * L
    return K