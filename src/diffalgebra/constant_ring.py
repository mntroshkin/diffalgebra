from typing import Optional
from fractions import Fraction


type Rational = int | Fraction


class ConstantRing[ElementType]:
    ring_name: str

    def __init__(self, ring_name: Optional[str] = None):
        self.ring_name = ring_name if ring_name else "anonymous constant ring"

    def is_element(self, expression) -> bool:
        raise NotImplementedError
    
    def promote(self, expression) -> ElementType:
        raise NotImplementedError
    
    def get_constant(self, const_name) -> ElementType:
        raise NotImplementedError


class Rationals(ConstantRing[Rational]):
    def __init__(self):
        ring_name = "QQ"
        super().__init__(ring_name)

    def is_element(self, expression) -> bool:
        return isinstance(expression, (int, Fraction))

    def promote(self, expression):
        if not self.is_element(expression):
            raise TypeError(f"Expression {expression} is not element of {self.ring_name} and can't be promoted")
        return Fraction(expression)
    
    def get_constant(self, const_name):
        raise ValueError
    
QQ = Rationals()
