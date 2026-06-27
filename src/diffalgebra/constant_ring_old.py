from fractions import Fraction
from typing import Sequence, Optional

from .exceptions import SymbolNameError, RingMismatchError

type Rational = int | Fraction
type Constant = int | Fraction | ConstantPolynomial


class ConstantMonomial:
    _ring: ConstantRing
    _exponents: tuple[int, ...]
    _coefficient: Rational

    def __init__(self, ring: ConstantRing,
                exponents: Optional[Sequence[int]], coefficient: Rational = 1):
        self._ring = ring
        self._coefficient = coefficient
        if exponents is None:
            self._exponents = tuple(0 for gen in self._ring._generators)
        else:
            self._exponents = tuple(exponents)

    def __mul__(self, other: ConstantMonomial) -> ConstantMonomial:
        l = len(self._ring._gen_names)
        product_exponents = tuple(self._exponents[i] + other._exponents[i] for i in range(l))
        product_coefficient = self._coefficient * other._coefficient
        return ConstantMonomial(ring=self._ring,
                                exponents=product_exponents,
                                coefficient=product_coefficient)

    def _d(self, var: ConstantGenerator) -> ConstantMonomial:
        index = self._ring._gen_names.index(var._gen_name)
        exponent = self._exponents[index]
        if exponent == 0:
            return ConstantMonomial(self._ring, exponents=None, coefficient=0)
        else:
            new_exponents = tuple(exp - 1 if i == index else exp for i, exp in enumerate(self._exponents))
            return ConstantMonomial(self._ring, exponents=new_exponents, coefficient=self._coefficient * exponent)

    def _int(self, var: ConstantGenerator) -> ConstantMonomial:
        index = self._ring._gen_names.index(var._gen_name)
        factor = Fraction(1, self._exponents[index] + 1)
        new_exponents = tuple(exp + 1 if i == index else exp for i, exp in enumerate(self._exponents))
        return ConstantMonomial(self._ring, exponents=new_exponents, coefficient=factor * self._coefficient)

    
    def __eq__(self, other) -> bool:
        if isinstance(other, ConstantMonomial):
            return (self._ring, self._exponents, self._coefficient) == (other._ring, other._exponents, other._coefficient)
        return False
    
    def __str__(self) -> str:
        if self._exponents == tuple(0 for gen in self._ring._generators):
            return str(self._coefficient)
        
        factors = []
        for i, gen_name in enumerate(self._ring._gen_names):
            if self._exponents[i] == 1:
                factors.append(gen_name)
            if self._exponents[i] > 1:
                factors.append(f"{gen_name}^{self._exponents[i]}")
        
        if self._coefficient == 1:
            return "*".join(factors)
        elif self._coefficient == -1:
            return "-" + "*".join(factors)
        else:
            return str(self._coefficient) + "*" + "*".join(factors)


class ConstantPolynomial:
    _ring: ConstantRing
    _terms: list[ConstantMonomial]

    def _normalize_terms(self, terms: Sequence[ConstantMonomial]) -> list[ConstantMonomial]:
        terms = list(terms)
        terms.sort(key=lambda term: term._exponents, reverse=True)
        normalized_terms = []
        i = 0
        while i < len(terms):
            exponents = terms[i]._exponents
            coefficient_sum = 0
            while i < len(terms) and terms[i]._exponents == exponents:
                coefficient_sum += terms[i]._coefficient
                i += 1
            if coefficient_sum != 0:
                normalized_terms.append(ConstantMonomial(self._ring, exponents, coefficient_sum))
        return normalized_terms

    def __init__(self, ring: ConstantRing, terms: Sequence[ConstantMonomial]):
        self._ring = ring
        self._terms = self._normalize_terms(terms)

    def __add__(self, other):
        if isinstance(other, ConstantPolynomial) and self._ring != other._ring:
            raise RingMismatchError
        if self._ring.is_element(other):
            if not isinstance(other, ConstantPolynomial):
                other = self._ring.promote(other)
            return ConstantPolynomial(self._ring, terms=self._terms + other._terms)
        else:
            return NotImplemented
        
    def __radd__(self, other):
        return self + other
    
    def __mul__(self, other):
        if isinstance(other, ConstantPolynomial) and self._ring != other._ring:
            raise RingMismatchError
        if self._ring.is_element(other):
            if not isinstance(other, ConstantPolynomial):
                other = self._ring.promote(other)
            product_terms = [term1 * term2 for term1 in self._terms for term2 in other._terms]
            return ConstantPolynomial(self._ring, terms=product_terms)
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
    
    def __pow__(self, other) -> ConstantPolynomial:
        if not isinstance(other, int):
            raise TypeError
        if other < 0:
            raise ValueError
        if other == 0:
            return self._ring.promote(1)
        if other % 2 == 0:
            sqrt = self ** (other // 2)
            return sqrt * sqrt
        else:
            return self * self ** (other - 1)
        
    def d(self, var: ConstantGenerator, order: int = 1) -> ConstantPolynomial:
        if order < 0:
            raise ValueError
        if order == 0:
            return self
        if order == 1:
            if self._ring != var._ring:
                raise RingMismatchError
            return ConstantPolynomial(ring=self._ring,
                                    terms = [term._d(var) for term in self._terms])
        return self.d(var, order - 1).d(var)

    def int(self, var: ConstantGenerator) -> ConstantPolynomial:
        return ConstantPolynomial(ring=self._ring, terms=[term._int(var) for term in self._terms])
            
    def __eq__(self, other):
        if self._ring.is_element(other):
            if not isinstance(other, ConstantPolynomial):
                other = self._ring.promote(other)
            return self._terms == other._terms
        return False
    
    def __str__(self) -> str:
        if len(self._terms) == 0:
            return "0"
        terms_str = [str(term) for term in self._terms]
        for i in range(1, len(self._terms)):
            if terms_str[i][0] != "-":
                terms_str[i] = "+" + terms_str[i]
        return ''.join(terms_str)
    
    def parenthesis_str(self) -> str:
        return f"({str(self)})" if len(self._terms) > 1 else str(self)


class ConstantGenerator(ConstantPolynomial):
    _gen_name: str

    def __init__(self, ring: ConstantRing, gen_name: str):
        if gen_name not in ring._gen_names:
            raise SymbolNameError
        exponents = tuple(1 if gen_name == const else 0 for const in ring._gen_names)
        monomial = ConstantMonomial(ring, exponents)
        terms = [monomial]

        super().__init__(ring, terms)
        self._gen_name = gen_name

    def __hash__(self) -> int:
        return hash((id(self._ring), self._gen_name))


class ConstantRing:
    _name: str
    _description: str
    _gen_names: list[str]
    _generators: dict[str, ConstantGenerator]

    def __init__(self, generators: Sequence[str], ring_name: Optional[str] = None):
        self._gen_names = []
        for const_name in generators:
            if not const_name:
                raise SymbolNameError(f"Empty generator names are not allowed")
            if const_name in self._gen_names:
                raise SymbolNameError(f"Repeatin generator names are not allowed")
            self._gen_names.append(const_name)

        self._generators = {gen_name: ConstantGenerator(ring=self, gen_name=gen_name) for gen_name in self._gen_names}

        ring_description = f"QQ[{', '.join(self._gen_names)}]"

        self._name = str(ring_name)
        self._description_brief = ring_description
        self._description = f"{ring_name} = {ring_description}" if ring_name else ring_description

    def __str__(self) -> str:
        return self._description

    def is_element(self, expression) -> bool:
        if isinstance(expression, (int, Fraction)):
            return True
        elif isinstance(expression, ConstantPolynomial) and expression._ring == self:
            return True
        else:
            return False
        
    def is_generator(self, variable) -> bool:
        if isinstance(variable, ConstantGenerator) and variable._ring == self:
            return True
        else:
            return False

    def promote(self, expression) -> ConstantPolynomial:
        if not self.is_element(expression):
            raise TypeError(f"Expression {expression} is not an element of {self._name} and can't be promoted")
        if isinstance(expression, ConstantPolynomial):
            return expression
        elif isinstance(expression, (int, Fraction)):
            monomial = ConstantMonomial(ring=self, exponents=None, coefficient=expression)
            return ConstantPolynomial(ring=self, terms=[monomial])
        raise TypeError
    
    def gen(self, gen_name: str) -> ConstantGenerator:
        if not gen_name in self._gen_names:
            raise SymbolNameError(f"{gen_name} is a not a generator of {self}")
        else:
            return self._generators[gen_name]
    
    def gens(self, *gen_names: str) -> tuple[ConstantGenerator, ...]:
        if len(gen_names) == 0:
            return tuple(self._generators.values())
        else:
            return tuple(self.gen(gen_name) for gen_name in gen_names)


class ConstantPolyRing(ConstantRing):
    def __init__(self, constants: Sequence[str], ring_name: str | None = None):
        if not constants:
            raise SymbolNameError(f"The list of ring generators must be non-empty")
        super().__init__(constants, ring_name)

class Rationals(ConstantRing):
    def __init__(self):
        super().__init__(generators=[], ring_name="QQ")
        self._name = "QQ"
        self._description = "QQ"
        self._description_brief = "QQ"

QQ = Rationals()