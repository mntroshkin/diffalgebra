from fractions import Fraction
from typing import Sequence, Optional

from .exceptions import SymbolNameError, RingMismatchError


type Constant = int | Fraction | ConstantPolynomial
type Rational = int | Fraction

type Monomial = tuple[int, ...]
type Term = tuple[Monomial, Rational]

def multiply_terms(left: Term, right: Term) -> Term:
    monomial_left, coefficient_left = left
    monomial_right, coefficient_right = right
    monomial = tuple(exp_left + exp_right for exp_left, exp_right in zip(monomial_left, monomial_right))
    coefficient = coefficient_left * coefficient_right
    return monomial, coefficient

def partial_term(term: Term, var_index: int) -> Term:
    monomial, coefficient = term
    exponent = monomial[var_index] 
    if exponent == 0:
        return (tuple(0 for var in monomial), 0)
    new_monomial = tuple(exp - 1 if i == var_index else exp for i, exp in enumerate(monomial))
    return (new_monomial, coefficient * exponent)

def integrate_term(term: Term, var_index: int) -> Term:
    monomial, coefficient = term
    factor = Fraction(1, monomial[var_index] + 1)
    new_monomial = tuple(exp + 1 if i == var_index else exp for i, exp in enumerate(monomial))
    return (new_monomial, coefficient * factor)

def str_term(ring: ConstantRing, term: Term) -> str:
    monomial, coefficient = term
    if monomial == tuple(0 for gen in ring._generators):
        return str(coefficient)
    
    factors = []
    for i, gen_name in enumerate(ring._gen_names):
        if monomial[i] == 1:
            factors.append(gen_name)
        if monomial[i] > 1:
            factors.append(f"{gen_name}^{monomial[i]}")
    
    if coefficient == 1:
        return "*".join(factors)
    elif coefficient == -1:
        return "-" + "*".join(factors)
    else:
        return str(coefficient) + "*" + "*".join(factors)

def normalize_terms(terms: list[Term]) -> list[Term]:
    terms.sort(reverse=True)
    normalized_terms: list[Term] = []
    i = 0
    while i < len(terms):
        monomial = terms[i][0]
        coefficient_sum = 0
        while i < len(terms) and terms[i][0] == monomial:
            coefficient_sum += terms[i][1]
            i += 1
        if coefficient_sum != 0:
            normalized_terms.append((monomial, coefficient_sum))
    return normalized_terms


class ConstantPolynomial:
    _ring: ConstantRing
    _terms: list[Term]

    def __init__(self, ring: ConstantRing, terms: list[Term]):
        self._ring = ring
        self._terms = normalize_terms(terms)

    def __add__(self, other):
        if isinstance(other, ConstantPolynomial) and self._ring != other._ring:
            raise RingMismatchError
        if self._ring.is_element(other):
            if not isinstance(other, ConstantPolynomial):
                other = self._ring.promote(other)
            return ConstantPolynomial(ring=self._ring, terms = self._terms + other._terms)
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
            product_terms = [multiply_terms(left, right) for left in self._terms for right in other._terms]
            return ConstantPolynomial(ring=self._ring, terms=product_terms)
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
            var_index = self._ring._gen_names.index(var._gen_name)
            return ConstantPolynomial(ring=self._ring,
                                    terms = [partial_term(term, var_index) for term in self._terms])
        return self.d(var, order - 1).d(var)
    
    def integral(self, var: ConstantGenerator) -> ConstantPolynomial:
        var_index = self._ring._gen_names.index(var._gen_name)
        return ConstantPolynomial(ring=self._ring,
                                  terms = [integrate_term(term, var_index) for term in self._terms])

    def __eq__(self, other):
        if self._ring.is_element(other):
            if not isinstance(other, ConstantPolynomial):
                other = self._ring.promote(other)
            return self._terms == other._terms
        return False

    def __str__(self) -> str:
        if len(self._terms) == 0:
            return "0"
        terms_str = [str_term(self._ring, term) for term in self._terms]
        for i in range(1, len(self._terms)):
            if terms_str[i][0] != "-":
                terms_str[i] = "+" + terms_str[i]
        return ''.join(terms_str)
    
    def parenthesis_str(self) -> str:
        return f"({str(self)})" if len(self._terms) > 1 else str(self)
    
    def __lt__(self, other: Constant) -> bool:
        if not self._ring.is_element(other):
            raise RingMismatchError
        other = self._ring.promote(other)
        return self._terms < other._terms

def partial(expression: Constant, var: ConstantGenerator, order: int = 1) -> Constant:
    if order < 0:
        raise ValueError
    if order == 0:
        return expression
    if isinstance(expression, (int, Fraction)):
        return 0
    if isinstance(expression, ConstantPolynomial):
        return expression.d(var, order)

class ConstantGenerator(ConstantPolynomial):
    _gen_name: str

    def __init__(self, ring: ConstantRing, gen_name: str):
        if gen_name not in ring._gen_names:
            raise SymbolNameError
        monomial = tuple(1 if gen_name == const else 0 for const in ring._gen_names)
        super().__init__(ring, terms=[(monomial, 1)])
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
            monomial = tuple(0 for gen in self._gen_names)
            return ConstantPolynomial(ring=self, terms=[(monomial, expression)])
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