from typing import Optional, Sequence
from fractions import Fraction

from .constant_ring import ConstantRing, QQ
from .constant_poly_ring import ConstantPolynomial

type ConstantCoefficient = int | Fraction | ConstantPolynomial

type DiffFactor = tuple[int, int]
type DiffFactors = list[DiffFactor]

class DifferentialMonomial:
    ring: DifferentialRing
    factors: tuple[DiffFactors, ...]

    @staticmethod
    def _normalize_factors(factors: DiffFactors) -> DiffFactors:
        factors.sort()
        normalized_factors: DiffFactors = []
        i = 0
        while i < len(factors):
            derivative_order = factors[i][0]
            total_power = 0
            while i < len(factors) and factors[i][0] == derivative_order:
                total_power += factors[i][1]
                i += 1
            if total_power != 0:
                normalized_factors.append((derivative_order, total_power))
        return normalized_factors

    def __init__(self, ring: DifferentialRing,
                 factors: Optional[Sequence[DiffFactors]],
                 coefficient: ConstantCoefficient = 1):
        self.ring = ring
        self.coefficient = self.ring.base_ring.promote(coefficient)

        if factors is None:
            self.factors = tuple([] for func in self.ring.functions)
        else:
            self.factors = tuple(self._normalize_factors(func_factors) for func_factors in factors)

    def __mul__(self, other: DifferentialMonomial) -> DifferentialMonomial:
        l = len(self.ring.functions)
        product_factors = tuple(self.factors[i] + other.factors[i] for i in range(l))
        product_coefficient = self.coefficient * other.coefficient
        return DifferentialMonomial(ring=self.ring,
                                factors=product_factors,
                                coefficient=product_coefficient)

    def __eq__(self, other) -> bool:
        if isinstance(other, DifferentialMonomial):
            return (self.ring, self.factors, self.coefficient) == (other.ring, other.factors, other.coefficient)
        raise TypeError
    
    def __str__(self) -> str:
        if self.factors == tuple([] for func_name in self.ring.functions):
            return str(self.coefficient)
        
        factors = []
        for i, func_name in enumerate(self.ring.functions):
            for derivative_order, power in self.factors[i]:
                if power == 1:
                    if derivative_order == 0:
                        factors.append(func_name)
                    elif derivative_order <= 3:
                        factors.append(f"{func_name}_{"x" * derivative_order}")
                    else:
                        factors.append(f"{func_name}_{derivative_order}")
                else:
                    if derivative_order == 0:
                        factors.append(f"{func_name}^{power}")
                    elif derivative_order <= 3:
                        factors.append(f"({func_name}_{"x" * derivative_order})^{power}")
                    else:
                        factors.append(f"({func_name}_{derivative_order})^{power}")
        
        if self.coefficient == 1:
            return "*".join(factors)
        elif self.coefficient == -1:
            return "-" + "*".join(factors)
        else:
            if isinstance(self.coefficient, ConstantPolynomial):
                coefficient_str = self.coefficient.parenthesis_str()
            else:
                coefficient_str = str(self.coefficient)
            return coefficient_str + "*" + "*".join(factors)
        
    def diff(self) -> list[DifferentialMonomial]:
        diff_terms: list[DifferentialMonomial] = []

        l = len(self.ring.functions)
        for i in range(l):
            for j, factor in enumerate(self.factors[i]):
                derivative_order, power = factor
                new_factors_ith = self.factors[i].copy()
                new_factors_ith[j] = (derivative_order, power - 1)
                new_factors_ith.append((derivative_order + 1, 1))
                new_factors = tuple(new_factors_ith if k == i else self.factors[k] for k in range(l))
                diff_terms.append(DifferentialMonomial(ring=self.ring,
                                                       coefficient=self.coefficient * power,
                                                       factors=new_factors))
        return diff_terms


class DifferentialPolynomial:
    ring: DifferentialRing
    terms: list[DifferentialMonomial]

    def _normalize_terms(self, terms: Sequence[DifferentialMonomial]) -> list[DifferentialMonomial]:
        terms = list(terms)
        terms.sort(key=lambda term: term.factors)
        normalized_terms = []
        i = 0
        while i < len(terms):
            factors = terms[i].factors
            coefficient_sum = 0
            while i < len(terms) and terms[i].factors == factors:
                coefficient_sum += terms[i].coefficient
                i += 1
            if coefficient_sum != 0:
                normalized_terms.append(DifferentialMonomial(self.ring, factors, coefficient_sum))
        return normalized_terms

    def __init__(self, ring: DifferentialRing, terms: Sequence[DifferentialMonomial]):
        self.ring = ring
        self.terms = self._normalize_terms(terms)

    def __add__(self, other):
        if self.ring.is_element(other):
            if not isinstance(other, DifferentialPolynomial):
                other = self.ring.promote(other)
            return DifferentialPolynomial(self.ring, terms=self.terms + other.terms)
        else:
            return NotImplemented
        
    def __radd__(self, other):
        return self + other
    
    def __mul__(self, other):
        if self.ring.is_element(other):
            if not isinstance(other, DifferentialPolynomial):
                other = self.ring.promote(other)
            product_terms = [term1 * term2 for term1 in self.terms for term2 in other.terms]
            return DifferentialPolynomial(self.ring, terms=product_terms)
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
        else:
            return self * (self ** (other - 1))
            
    def __eq__(self, other):
        if self.ring.is_element(other):
            if not isinstance(other, DifferentialPolynomial):
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

    def diff(self, order: int = 1) -> DifferentialPolynomial:
        if order < 0:
            raise ValueError
        if order == 0:
            return self
        if order == 1:
            return DifferentialPolynomial(ring=self.ring, 
                                          terms=[diff_term for term in self.terms for diff_term in term.diff()])
        else:
            return self.diff(order - 1).diff()
        
def diff(expression: ConstantCoefficient | DifferentialPolynomial, order: int = 1):
    if order < 0:
        raise ValueError
    if order == 0:
        return expression
    if order >= 1:
        if isinstance(expression, (int, Fraction, ConstantPolynomial)):
            return 0
        else:
            return expression.diff(order)

class DifferentialRing:
    ring_name: str
    base_ring: ConstantRing[ConstantCoefficient]
    functions: list[str]

    def __init__(self, functions: list[str], base_ring: ConstantRing = QQ, ring_name: Optional[str] = None):
        self.base_ring = base_ring
        self.ring_name = ring_name if ring_name else "anonymous differential ring ring"
        
        self.functions = []
        for func_name in functions:
            if not func_name:
                raise ValueError(f"Error defining {self.ring_name}: Empty function names are not allowed")
            if func_name in self.functions:
                raise ValueError(f"Error defining {self.ring_name}: Repeating function names are not allowed")
            self.functions.append(func_name)

    def is_element(self, expression) -> bool:
        if self.base_ring.is_element(expression):
            return True
        if isinstance(expression, DifferentialPolynomial) and expression.ring == self:
            return True
        return False
    
    def promote(self, expression) -> DifferentialPolynomial:
        if not self.is_element(expression):
            raise TypeError(f"Expression {expression} is not element of {self.ring_name} and can't be promoted")
        if isinstance(expression, DifferentialPolynomial):
            return expression
        if isinstance(expression, (int, Fraction, ConstantPolynomial)):
            expression = self.base_ring.promote(expression)
            monomial = DifferentialMonomial(ring=self, factors=None, coefficient=expression)
            return DifferentialPolynomial(ring=self, terms=[monomial])
        raise TypeError
    
    def get_constant(self, const_name: str) -> ConstantCoefficient:
        return self.base_ring.get_constant(const_name)
    
    def get_function(self, func_name: str) -> DifferentialPolynomial:
        if func_name not in self.functions:
            raise ValueError
        factors = tuple([(0, 1)] if func_name == func else [] for func in self.functions)
        monomial = DifferentialMonomial(ring=self, factors=factors)
        return DifferentialPolynomial(ring=self, terms=[monomial])