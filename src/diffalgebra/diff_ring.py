from typing import Optional, Sequence
from fractions import Fraction

from .exceptions import RingMismatchError, SymbolNameError
from .constant_ring import ConstantRing, ConstantPolynomial, ConstantGenerator, QQ, Constant, partial as partial_for_const


type Expression = int | Fraction | ConstantPolynomial | DifferentialPolynomial
type Generator = ConstantGenerator | FuncGenerator

type DiffFactor = tuple[int, int]
type DiffFactors = list[DiffFactor]

type DiffMonomial = tuple[DiffFactors, ...]
type DiffTerm = tuple[DiffMonomial, Constant]

def normalize_factors(factors: DiffFactors) -> DiffFactors:
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

def normalize_terms(terms: list[DiffTerm]) -> list[DiffTerm]:
    for i, term in enumerate(terms):
        monomial, coefficient = term
        monomial = tuple(factors for factors in monomial)
        terms[i] = (monomial, coefficient)
    terms.sort()
    normalized_terms: list[DiffTerm] = []
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

def multiply_terms(left: DiffTerm, right: DiffTerm) -> DiffTerm:
    monomial_left, coefficient_left = left
    monomial_right, coefficient_right = right
    monomial = tuple(normalize_factors(factors_left + factors_right)
                     for factors_left, factors_right
                     in zip(monomial_left, monomial_right))
    coefficient = coefficient_left * coefficient_right
    return monomial, coefficient

def diff_termwise(term: DiffTerm) -> list[DiffTerm]:
    monomial, coefficient = term
    diff_terms: list[DiffTerm] = []
    for i in range(len(monomial)):
        for j, factor in enumerate(monomial[i]):
            derivative_order, power = factor
            new_factors_ith = monomial[i].copy()
            new_factors_ith[j] = (derivative_order, power - 1)
            new_factors_ith.append((derivative_order + 1, 1))
            new_factors_ith = normalize_factors(new_factors_ith)
            new_factors = tuple(new_factors_ith if k == i else factors 
                                for k, factors in enumerate(monomial))
            diff_terms.append((new_factors, coefficient * power))
    return diff_terms


def exponent_in_term(term: DiffTerm, var_index: int, derivative_order: int) -> int:
    monomial, coefficient = term
    for derivative, power in monomial[var_index]:
        if derivative == derivative_order:
            return power
    return 0


def partial_termwise(term: DiffTerm, var_index: int, derivative_order: int) -> DiffTerm:
    monomial, coefficient = term

    new_coefficient = 0
    new_factors_ith = monomial[var_index].copy()
    for j, (derivative, power) in enumerate(new_factors_ith):
        if derivative == derivative_order:
            new_factors_ith[j] = (derivative, power - 1)
            new_coefficient = coefficient * power
    new_factors_ith = normalize_factors(new_factors_ith)
    new_monomial = tuple(new_factors_ith if k == var_index else factors 
                        for k, factors in enumerate(monomial))
    if new_coefficient == 0:
        return (tuple([] for var in monomial), 0)
    return (new_monomial, new_coefficient)
    

def str_term(ring: DifferentialRing, term: DiffTerm) -> str:
    monomial, coefficient = term
    if monomial == tuple([] for func_name in ring._func_names):
        return str(coefficient)
    
    factors = []
    for i, func_name in enumerate(ring._func_names):
        for derivative_order, power in monomial[i]:
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
    
    if coefficient == 1:
        return "*".join(factors)
    elif coefficient == -1:
        return "-" + "*".join(factors)
    else:
        if isinstance(coefficient, ConstantPolynomial):
            coefficient_str = coefficient.parenthesis_str()
        else:
            coefficient_str = str(coefficient)
        return coefficient_str + "*" + "*".join(factors)
    

class DifferentialPolynomial:
    _ring: DifferentialRing
    _terms: list[DiffTerm]

    def __init__(self, ring: DifferentialRing, terms: list[DiffTerm]):
        self._ring = ring
        self._terms = normalize_terms(terms)

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
            product_terms = [multiply_terms(left, right) for left in self._terms for right in other._terms]
            return DifferentialPolynomial(ring=self._ring, terms=product_terms)     
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
    
    def __pow__(self, other: int) -> DifferentialPolynomial:
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
        terms_str = [str_term(self._ring, term) for term in self._terms]
        for i in range(1, len(self._terms)):
            if terms_str[i][0] != "-":
                terms_str[i] = "+" + terms_str[i]
        return ''.join(terms_str)

    def __repr__(self) -> str:
        return str(self)
    
    def diff(self, order: int = 1) -> DifferentialPolynomial:
        if order < 0:
            raise ValueError
        if order == 0:
            return self
        if order == 1:
            return DifferentialPolynomial(ring=self._ring, 
                                          terms=[diff_term 
                                                 for term in self._terms 
                                                 for diff_term in diff_termwise(term)])
        else:
            return self.diff(order - 1).diff()
        
    def d(self, var: Generator, order: int = 1) -> DifferentialPolynomial:
        if order < 0:
            raise ValueError
        if order == 0:
            return self
        if order == 1:
            if isinstance(var, ConstantGenerator):
                if self._ring._base_ring != var._ring:
                    raise RingMismatchError
                new_terms = [(monomial, partial_for_const(coefficient, var)) 
                             for (monomial, coefficient) in self._terms]
                return DifferentialPolynomial(ring=self._ring, terms=new_terms)
            if isinstance(var, FuncGenerator):
                if self._ring != var._ring:
                    raise RingMismatchError
                var_index = self._ring._func_names.index(var._func_name)
                derivative_order = var._derivative
                return DifferentialPolynomial(ring=self._ring,
                                    terms = [partial_termwise(term, var_index, derivative_order) 
                                             for term in self._terms])
        return self.d(var, order - 1).d(var)

    def _highest_derivative(self, var: FuncGenerator) -> int:
        if self._ring != var._ring:
            raise RingMismatchError
        var_index = self._ring._func_names.index(var._func_name)
        highest_derivative = -1
        for term in self._terms:
            highest_for_term = max((derivative for derivative, power in term[0][var_index]), default=-1)
            highest_derivative = max(highest_derivative, highest_for_term)
        return highest_derivative
    
    def delta(self, var: FuncGenerator) -> DifferentialPolynomial:
        if not self._ring.is_generator(var):
            raise ValueError("{var} is not a generator of {self._ring}")
        d = self._highest_derivative(var)
        result_terms: list[DiffTerm] = []
        for i in range(d + 1):
            summand = (-1) ** i * self.d(var[i]).diff(order=i)
            result_terms.extend(self._ring.promote(summand)._terms)
        return DifferentialPolynomial(ring=self._ring, terms=result_terms)

    def coefficient(self, monomial: Expression) -> Constant:
        monomial = self._ring.promote(monomial)
        if len(monomial._terms) != 1:
            raise ValueError
        exponents, coef = monomial._terms[0]
        for term in self._terms:
            if term[0] == exponents:
                return term[1]
        return 0
    
    def integral(self) -> DifferentialPolynomial | None:
        if self == 0:
            return self._ring.promote(0)
        if self.coefficient(1) != 0:
            return
        for var_index, var in enumerate(self._ring.gens()):
            if self._highest_derivative(var) != -1:
                new_integrand = self._ring.promote(0)
                integral_found = self._ring.promote(0)
                k = self._highest_derivative(var)
                if k == 0:
                    return
                for term in self._terms:
                    monomial, coefficient = term
                    j = exponent_in_term(term, var_index, k)
                    if j == 0:
                        new_integrand += DifferentialPolynomial(self._ring, terms=[term])
                    elif j > 1:
                        return
                    else:
                        j = exponent_in_term(term, var_index, k - 1)
                        factor = Fraction(1, j + 1) * var[k - 1] ** (j + 1)
                        h_factors_ith = [(derivative, power)
                                           for derivative, power in monomial[var_index]
                                           if derivative < k - 1]
                        h_monomial = tuple(h_factors_ith
                                           if k == var_index else factors
                                           for k, factors in enumerate(monomial))
                        h = DifferentialPolynomial(self._ring, terms=[(h_monomial, coefficient)])
                        integral_found += h * factor
                        new_integrand -= h.diff() * factor
                integral_remaining = new_integrand.integral()
                if integral_remaining is None:
                    return
                return integral_found + integral_remaining
                        

def total_derivative(expression: Expression, order: int = 1) -> Expression:
    if order < 0:
        raise ValueError
    if order == 0:
        return expression

    if isinstance(expression, (int, Fraction, ConstantPolynomial)):
        return 0
    else:
        return expression.diff(order)

def partial_derivative(expression: Expression,
                       var: Generator,
                       order: int = 1) -> Expression:
    if isinstance(expression, (int, Fraction)):
        return 0
    if isinstance(expression, ConstantPolynomial) and isinstance(var, ConstantGenerator):
        return expression.d(var, order)
    if isinstance(expression, DifferentialPolynomial) and isinstance(var, (ConstantGenerator, FuncGenerator)):
        return expression.d(var, order)
    raise TypeError


class FuncGenerator(DifferentialPolynomial):
    _func_name: str
    _derivative: int

    def __init__(self, ring: DifferentialRing, func_name: str, derivative: int = 0):
        if func_name not in ring._func_names:
            raise SymbolNameError
        if derivative < 0:
            raise ValueError
        monomial = tuple([(derivative, 1)] if func_name == func else [] for func in ring._func_names)
        super().__init__(ring, terms=[(monomial, 1)])

        self._func_name = func_name
        self._derivative = derivative
    
    def diff(self, order: int = 1) -> FuncGenerator:
        if order < 0:
            raise ValueError
        else:
            return FuncGenerator(ring=self._ring, func_name=self._func_name, derivative=self._derivative + order)

    def __getitem__(self, key) -> FuncGenerator:
        if isinstance(key, int) and key >= 0:
            return self.diff(order=key)
        else:
            raise ValueError
        
    def __hash__(self) -> int:
        return hash((id(self._ring), self._func_name, self._derivative))


class DifferentialRing:
    _name: str
    _base_ring: ConstantRing
    _func_names: list[str]
    _generators: dict[str, FuncGenerator]

    def __init__(self, functions: Sequence[str], base_ring: ConstantRing = QQ, ring_name: Optional[str] = None):
        self._base_ring = base_ring
        
        self._func_names = []
        for func_name in functions:
            if not func_name:
                raise SymbolNameError(f"Empty function names are not allowed")
            if func_name in self._func_names:
                raise SymbolNameError(f"Repeating function names are not allowed")
            if func_name in self._base_ring._gen_names:
                raise SymbolNameError(f"Function name coincides with a constant name")
            self._func_names.append(func_name)

        self._generators = {func_name: FuncGenerator(ring=self, func_name=func_name) for func_name in self._func_names}
        
        ring_definition = f"{self._base_ring._description_brief}{{{ ', '.join(self._func_names) }}}"

        self._name = str(ring_name)
        self._description_brief = ring_definition
        self._description = f"{ring_name} = {ring_definition}" if ring_name else ring_definition

    def __str__(self) -> str:
        return self._description

    def __repr__(self) -> str:
        return f"{str(self)}: differential ring with id={id(self)}" 

    def is_element(self, expression) -> bool:
        if self._base_ring.is_element(expression):
            return True
        if isinstance(expression, DifferentialPolynomial) and expression._ring == self:
            return True
        return False
    
    def is_generator(self, variable) -> bool:
        if isinstance(variable, FuncGenerator) and variable._ring == self and variable._derivative == 0:
            return True
        else:
            return False
    
    def promote(self, expression) -> DifferentialPolynomial:
        if not self.is_element(expression):
            raise TypeError(f"Expression {expression} is not an element of {self._name} and can't be promoted")
        if isinstance(expression, DifferentialPolynomial):
            return expression
        if isinstance(expression, (int, Fraction, ConstantPolynomial)):
            monomial = tuple([] for func in self._func_names)
            return DifferentialPolynomial(ring=self, terms=[(monomial, expression)])
        raise TypeError
    
    def constants(self, *const_names: str) -> tuple[ConstantGenerator, ...]:
        return self._base_ring.gens(*const_names)
    
    def gen(self, func_name: str) -> FuncGenerator:
        if not func_name in self._func_names:
            raise SymbolNameError(f"{func_name} is a not a generator of {self}")
        else:
            return self._generators[func_name]
    
    def gens(self, *func_names: str) -> tuple[FuncGenerator, ...]:
        if len(func_names) == 0:
            return tuple(self._generators.values())
        else:
            return tuple(self.gen(func_name) for func_name in func_names)
