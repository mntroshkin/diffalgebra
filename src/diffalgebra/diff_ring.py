from typing import Optional, Sequence
from fractions import Fraction

from .exceptions import RingMismatchError, SymbolNameError
from .constant_ring import ConstantRing, QQ
from .constant_poly_ring import ConstantPolynomial, ConstantVariable

type ConstantCoefficient = int | Fraction | ConstantPolynomial

type DiffFactor = tuple[int, int]
type DiffFactors = list[DiffFactor]

class DifferentialMonomial:
    _ring: DifferentialRing
    _factors: tuple[DiffFactors, ...]
    _coefficient: ConstantCoefficient

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
        self._ring = ring
        self._coefficient = self._ring._base_ring.promote(coefficient)

        if factors is None:
            self._factors = tuple([] for func in self._ring._functions)
        else:
            self._factors = tuple(self._normalize_factors(func_factors) for func_factors in factors)

    def __mul__(self, other: DifferentialMonomial) -> DifferentialMonomial:
        l = len(self._ring._functions)
        product_factors = tuple(self._factors[i] + other._factors[i] for i in range(l))
        product_coefficient = self._coefficient * other._coefficient
        return DifferentialMonomial(ring=self._ring,
                                factors=product_factors,
                                coefficient=product_coefficient)

    def __eq__(self, other) -> bool:
        if isinstance(other, DifferentialMonomial):
            return (self._ring, self._factors, self._coefficient) == (other._ring, other._factors, other._coefficient)
        raise TypeError
    
    def __str__(self) -> str:
        if self._factors == tuple([] for func_name in self._ring._functions):
            return str(self._coefficient)
        
        factors = []
        for i, func_name in enumerate(self._ring._functions):
            for derivative_order, power in self._factors[i]:
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
        
        if self._coefficient == 1:
            return "*".join(factors)
        elif self._coefficient == -1:
            return "-" + "*".join(factors)
        else:
            if isinstance(self._coefficient, ConstantPolynomial):
                coefficient_str = self._coefficient.parenthesis_str()
            else:
                coefficient_str = str(self._coefficient)
            return coefficient_str + "*" + "*".join(factors)
        
    def _diff(self) -> list[DifferentialMonomial]:
        diff_terms: list[DifferentialMonomial] = []

        l = len(self._ring._functions)
        for i in range(l):
            for j, factor in enumerate(self._factors[i]):
                derivative_order, power = factor
                new_factors_ith = self._factors[i].copy()
                new_factors_ith[j] = (derivative_order, power - 1)
                new_factors_ith.append((derivative_order + 1, 1))
                new_factors = tuple(new_factors_ith if k == i else self._factors[k] for k in range(l))
                diff_terms.append(DifferentialMonomial(ring=self._ring,
                                                       coefficient=self._coefficient * power,
                                                       factors=new_factors))
        return diff_terms

    def _d(self, var: FuncVariable) -> DifferentialMonomial:
        i = self._ring._functions.index(var._func_name)
        l = len(self._ring._functions)

        new_coefficient = 0
        new_factors_ith = self._factors[i].copy()
        for j, (derivative, exp) in enumerate(new_factors_ith):
            if derivative == var._derivative:
                new_factors_ith[j] = (derivative, exp - 1)
                new_coefficient = self._coefficient * exp
        new_factors = tuple(new_factors_ith if k == i else self._factors[k] for k in range(l))
        return DifferentialMonomial(ring=self._ring, factors=new_factors, coefficient=new_coefficient)

class DifferentialPolynomial:
    _ring: DifferentialRing
    _terms: list[DifferentialMonomial]

    def _normalize_terms(self, terms: Sequence[DifferentialMonomial]) -> list[DifferentialMonomial]:
        terms = list(terms)
        terms.sort(key=lambda term: term._factors)
        normalized_terms = []
        i = 0
        while i < len(terms):
            factors = terms[i]._factors
            coefficient_sum = 0
            while i < len(terms) and terms[i]._factors == factors:
                coefficient_sum += terms[i]._coefficient
                i += 1
            if coefficient_sum != 0:
                normalized_terms.append(DifferentialMonomial(self._ring, factors, coefficient_sum))
        return normalized_terms

    def __init__(self, ring: DifferentialRing, terms: Sequence[DifferentialMonomial]):
        self._ring = ring
        self._terms = self._normalize_terms(terms)

    def __add__(self, other):
        if isinstance(other, DifferentialPolynomial) and other._ring != self._ring:
            raise RingMismatchError
        if isinstance(other, ConstantPolynomial) and other._ring != self._ring._base_ring:
            raise RingMismatchError
        if self._ring.is_element(other):
            if not isinstance(other, DifferentialPolynomial):
                other = self._ring.promote(other)
            return DifferentialPolynomial(self._ring, terms=self._terms + other._terms)
        else:
            return NotImplemented
        
    def __radd__(self, other):
        return self + other
    
    def __mul__(self, other):
        if isinstance(other, DifferentialPolynomial) and other._ring != self._ring:
            raise RingMismatchError
        if isinstance(other, ConstantPolynomial) and other._ring != self._ring._base_ring:
            raise RingMismatchError
        if self._ring.is_element(other):
            if not isinstance(other, DifferentialPolynomial):
                other = self._ring.promote(other)
            product_terms = [term1 * term2 for term1 in self._terms for term2 in other._terms]
            return DifferentialPolynomial(self._ring, terms=product_terms)
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
            return self._ring.promote(1)
        if other % 2 == 0:
            sqrt = self ** (other // 2)
            return sqrt * sqrt
        else:
            return self * self ** (other - 1)
            
    def __eq__(self, other):
        if self._ring.is_element(other):
            if not isinstance(other, DifferentialPolynomial):
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

    def diff(self, n: int = 1) -> DifferentialPolynomial:
        if n < 0:
            raise ValueError
        if n == 0:
            return self
        if n == 1:
            return DifferentialPolynomial(ring=self._ring, 
                                          terms=[diff_term for term in self._terms for diff_term in term._diff()])
        else:
            return self.diff(n - 1).diff()
        
    def d(self, var: ConstantVariable | FuncVariable) -> DifferentialPolynomial:
        if isinstance(var, ConstantVariable):
            if self._ring._base_ring != var._ring:
                raise RingMismatchError
            new_terms = [DifferentialMonomial(ring=self._ring, 
                                              factors=term._factors, 
                                              coefficient=term._coefficient.d(var)) for term in self._terms]
            return DifferentialPolynomial(ring=self._ring, terms=new_terms)
        if isinstance(var, FuncVariable):
            if self._ring != var._ring:
                raise RingMismatchError
            return DifferentialPolynomial(ring=self._ring,
                                  terms = [term._d(var) for term in self._terms])
            
        

def total_derivative(expression: ConstantCoefficient | DifferentialPolynomial, order: int = 1):
    if order < 0:
        raise ValueError
    if order == 0:
        return expression
    if order >= 1:
        if isinstance(expression, (int, Fraction, ConstantPolynomial)):
            return 0
        else:
            return expression.diff(order)


class FuncVariable(DifferentialPolynomial):
    _func_name: str
    _derivative: int

    def __init__(self, ring: DifferentialRing, func_name: str, derivative: int = 0):
        if func_name not in ring._functions:
            raise SymbolNameError
        if derivative < 0:
            raise ValueError
        factors = tuple([(derivative, 1)] if func_name == func else [] for func in ring._functions)
        monomial = DifferentialMonomial(ring=ring, factors=factors)
        terms = [monomial]

        super().__init__(ring, terms)
        self._func_name = func_name
        self._derivative = derivative
    
    def diff(self, n: int = 1) -> FuncVariable:
        if n < 0:
            raise ValueError
        else:
            return FuncVariable(ring=self._ring, func_name=self._func_name, derivative=self._derivative + n)

    def __getitem__(self, key) -> FuncVariable:
        if isinstance(key, int) and key >= 0:
            return self.diff(n=key)
        else:
            raise ValueError


class DifferentialRing:
    _ring_name: str
    _base_ring: ConstantRing[ConstantCoefficient]
    _functions: list[str]

    def __init__(self, functions: list[str], base_ring: ConstantRing = QQ, ring_name: Optional[str] = None):
        self._base_ring = base_ring
        self._ring_name = ring_name if ring_name else "anonymous differential ring"
        
        self._functions = []
        for func_name in functions:
            if not func_name:
                raise ValueError(f"Error defining {self._ring_name}: Empty function names are not allowed")
            if func_name in self._functions:
                raise ValueError(f"Error defining {self._ring_name}: Repeating function names are not allowed")
            self._functions.append(func_name)

    def is_element(self, expression) -> bool:
        if self._base_ring.is_element(expression):
            return True
        if isinstance(expression, DifferentialPolynomial) and expression._ring == self:
            return True
        return False
    
    def promote(self, expression) -> DifferentialPolynomial:
        if not self.is_element(expression):
            raise TypeError(f"Expression {expression} is not element of {self._ring_name} and can't be promoted")
        if isinstance(expression, DifferentialPolynomial):
            return expression
        if isinstance(expression, (int, Fraction, ConstantPolynomial)):
            expression = self._base_ring.promote(expression)
            monomial = DifferentialMonomial(ring=self, factors=None, coefficient=expression)
            return DifferentialPolynomial(ring=self, terms=[monomial])
        raise TypeError
    
    def get_constant(self, const_name: str) -> ConstantCoefficient:
        return self._base_ring.get_constant(const_name)
    
    def get_function(self, func_name: str) -> FuncVariable:
        return FuncVariable(ring=self, func_name=func_name)
