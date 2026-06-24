from fractions import Fraction
from typing import Sequence, Optional

from .constant_ring import ConstantRing

type Rational = int | Fraction


class ConstantMonomial:
    ring: ConstantPolyRing
    exponents: tuple[int, ...]
    coefficient: Rational

    def __init__(self, ring: ConstantPolyRing,
                exponents: Optional[Sequence[int]], coefficient: Rational = 1):
        self.ring = ring
        self.coefficient = coefficient
        if exponents is None:
            self.exponents = tuple(0 for const_name in self.ring.constants)
        else:
            self.exponents = tuple(exponents)

    def __mul__(self, other: ConstantMonomial) -> ConstantMonomial:
        l = len(self.ring.constants)
        product_exponents = tuple(self.exponents[i] + other.exponents[i] for i in range(l))
        product_coefficient = self.coefficient * other.coefficient
        return ConstantMonomial(ring=self.ring,
                                exponents=product_exponents,
                                coefficient=product_coefficient)
    
    def __eq__(self, other) -> bool:
        if isinstance(other, ConstantMonomial):
            return (self.ring, self.exponents, self.coefficient) == (other.ring, other.exponents, other.coefficient)
        raise TypeError
    
    def __str__(self) -> str:
        if self.exponents == tuple(0 for const_name in self.ring.constants):
            return str(self.coefficient)
        
        factors = []
        for i, const_name in enumerate(self.ring.constants):
            if self.exponents[i] == 1:
                factors.append(const_name)
            if self.exponents[i] > 1:
                factors.append(f"{const_name}^{self.exponents[i]}")
        
        if self.coefficient == 1:
            return "*".join(factors)
        elif self.coefficient == -1:
            return "-" + "*".join(factors)
        else:
            return str(self.coefficient) + "*" + "*".join(factors)


class ConstantPolynomial:
    ring: ConstantPolyRing
    terms: list[ConstantMonomial]

    def _normalize_terms(self, terms: Sequence[ConstantMonomial]) -> list[ConstantMonomial]:
        terms = list(terms)
        terms.sort(key=lambda term: term.exponents, reverse=True)
        normalized_terms = []
        i = 0
        while i < len(terms):
            exponents = terms[i].exponents
            coefficient_sum = 0
            while i < len(terms) and terms[i].exponents == exponents:
                coefficient_sum += terms[i].coefficient
                i += 1
            if coefficient_sum != 0:
                normalized_terms.append(ConstantMonomial(self.ring, exponents, coefficient_sum))
        return normalized_terms

    def __init__(self, ring: ConstantPolyRing, terms: Sequence[ConstantMonomial]):
        self.ring = ring
        self.terms = self._normalize_terms(terms)

    def __add__(self, other):
        if self.ring.is_element(other):
            if not isinstance(other, ConstantPolynomial):
                other = self.ring.promote(other)
            return ConstantPolynomial(self.ring, terms=self.terms + other.terms)
        else:
            return NotImplemented
        
    def __radd__(self, other):
        return self + other
    
    def __mul__(self, other):
        if self.ring.is_element(other):
            if not isinstance(other, ConstantPolynomial):
                other = self.ring.promote(other)
            product_terms = [term1 * term2 for term1 in self.terms for term2 in other.terms]
            return ConstantPolynomial(self.ring, terms=product_terms)
        else:
            return NotImplemented
        
    def __rmul__(self, other):
        return self * other
    
    def __neg__(self):
        return self * (-1)
    
    def __sub__(self, other):
        return self + other * (-1)
    
    def __rsub__(self, other):
        return self * (-1) + other
    
    def __pow__(self, other):
        if not isinstance(other, int):
            raise TypeError
        if other < 0:
            raise ValueError
        if other == 0:
            return 1
        if other % 2 == 0:
            sqrt = self ** (other // 2)
            return sqrt * sqrt
        else:
            return self * self ** (other - 1)
            
    def __eq__(self, other):
        if self.ring.is_element(other):
            if not isinstance(other, ConstantPolynomial):
                other = self.ring.promote(other)
            return self.terms == other.terms
        raise TypeError
    
    def __str__(self) -> str:
        if len(self.terms) == 0:
            return "0"
        terms_str = [str(term) for term in self.terms]
        for i in range(1, len(self.terms)):
            if terms_str[i][0] != "-":
                terms_str[i] = "+" + terms_str[i]
        return ''.join(terms_str)
    
    def parenthesis_str(self) -> str:
        return f"({str(self)})" if len(self.terms) > 1 else str(self)


class ConstantPolyRing(ConstantRing[ConstantPolynomial]):
    ring_name: str
    constants: list[str]

    def __init__(self, constants: Sequence[str], ring_name: Optional[str] = None):
        super().__init__(ring_name)

        self.constants = []
        for const_name in constants:
            if not const_name:
                raise ValueError(f"Error defining {self.ring_name}: Empty constant names are not allowed")
            if const_name in self.constants:
                raise ValueError(f"Error defining {self.ring_name}: Repeating constant names are not allowed")
            self.constants.append(const_name)

    def is_element(self, expression) -> bool:
        if isinstance(expression, (int, Fraction)):
            return True
        elif isinstance(expression, ConstantPolynomial) and expression.ring == self:
            return True
        else:
            return False

    def promote(self, expression) -> ConstantPolynomial:
        if not self.is_element(expression):
            raise TypeError(f"Expression {expression} is not element of {self.ring_name} and can't be promoted")
        if isinstance(expression, ConstantPolynomial):
            return expression
        elif isinstance(expression, (int, Fraction)):
            monomial = ConstantMonomial(ring=self, exponents=None, coefficient=expression)
            return ConstantPolynomial(ring=self, terms=[monomial])
        raise TypeError
    
    def get_constant(self, const_name: str):
        if const_name not in self.constants:
            raise ValueError
        exponents = tuple(1 if const_name == const else 0 for const in self.constants)
        monomial = ConstantMonomial(ring=self, exponents=exponents)
        return ConstantPolynomial(ring=self, terms=[monomial])