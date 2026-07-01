from typing import Optional

from .diff_ring import DifferentialRing, DifferentialPolynomial, FuncGenerator, Expression
from .exceptions import DefinitionError

class EvolutionOperator:
    _ring: DifferentialRing
    _mapping: tuple[DifferentialPolynomial]

    def __init__(self, ring: DifferentialRing,
                 mapping: dict[FuncGenerator, Expression],
                 name: Optional[str] = None):
        self._ring = ring
        _mapping = []
        for gen_name in ring._func_names:
            generator = ring.gen(gen_name)
            if generator not in mapping:
                raise DefinitionError(f"Generator {gen_name} has no defined image")
            image = mapping[generator]
            if not ring.is_element(image):
                raise TypeError(f"{image} is not an element of ring {ring}")
            _mapping.append(ring.promote(image))
        for generator in mapping.keys():
            if not ring.is_generator(generator):
                raise TypeError(f"{generator} is not a generator of {ring}")
        
        self._mapping = tuple(_mapping)
        self._name = name

    def apply(self, expression: Expression):
        if not self._ring.is_element(expression):
            raise TypeError(f"{expression} is not an element of {self._ring}")
        result = self._ring.promote(0)
        expression = self._ring.promote(expression)
        for i, var in enumerate(self._ring.gens()):
            factor = self._mapping[i]
            for j in range(expression._highest_derivative(var) + 1):
                result += expression.d(var[j]) * factor
                factor = factor.diff()
        return result
    
    def __call__(self, argument: Expression) -> DifferentialPolynomial:
        return self.apply(argument)