from fractions import Fraction
from typing import Type

class GeneralAlgebra:
    _empty_exponent = 0
    _default_exponent = 1
    CoefficientType = Fraction
    MonomialType: Type['GeneralMonomial']
    PolynomialType: Type['GeneralPolynomial']
    
    def __init__(self, variable_identifiers):
        self.variables = dict()
        for var_id in variable_identifiers:
            if var_id in self.variables.keys():
                raise KeyError("Variable identifiers must be unique within each algebra")
            self.variables[var_id] = Variable(self, var_id)

    def get_all_ids(self):
        return tuple([var_id for var_id in self.variables.keys()])

    def is_coefficient(self, expression):
        return isinstance(expression, int) or isinstance(expression, Fraction)

    def is_element(self, expression):
        if self.is_coefficient(expression):
            return True
        if isinstance(expression, self.PolynomialType) and expression.algebra == self:
            return True
        return False

    def get_variable(self, var_id):
        return self.get_monomial({var_id : self._default_exponent})

    def get_monomial(self, exponents):
        return self.MonomialType(self, exponents).to_polynomial()

    def normalize_coefficient(self, coefficient):
        return coefficient

    def normalize_exponent(self, exponent):
        pass

    def normalize_exponent_vector(self, exponent_vector):
        normalized_exponent_vector = self.empty_exponent_vector()
        if not isinstance(exponent_vector, dict):
            raise TypeError
        for var_id in exponent_vector.keys():
            if not var_id in self.variables.keys():
                raise KeyError("Variable identifier not found in the algebra")
            normalized_exponent_vector[var_id] = self.normalize_exponent(exponent_vector[var_id])
        return normalized_exponent_vector

    def empty_exponent_vector(self):
        return {var_id: self._empty_exponent for var_id in self.variables.keys()}


class Variable:
    def __init__(self, algebra, identifier, representation=None):
        representation = representation or identifier
        if not isinstance(algebra, GeneralAlgebra):
            raise TypeError
        if not isinstance(identifier, str) or not isinstance(representation, str):
            raise TypeError
        if identifier == "" or representation == "":
            raise ValueError
        self.algebra = algebra
        self.identifier = identifier
        self.representation = representation or identifier

    def set_representation(self, new_rep):
        if not isinstance(new_rep, str):
            raise TypeError
        if not new_rep:
            raise ValueError("Representation string must be non-empty.")
        self.representation = new_rep

    def monomial(self, exponent):
        return self.algebra.get_monomial({self.identifier : exponent})

    def to_polynomial(self):
        return self.algebra.get_variable(self.identifier)


class GeneralMonomial:
    AlgebraType = GeneralAlgebra

    def __init__(self, algebra, exponents, coefficient=1):
        if not isinstance(algebra, self.AlgebraType):
            raise TypeError
        self.algebra = algebra
        if not self.algebra.is_coefficient(coefficient):
            raise TypeError
        if isinstance(exponents, dict):
            if coefficient != 0:
                self.exponents = self.algebra.normalize_exponent_vector(exponents)
            else:
                self.exponents = self.algebra.empty_exponent_vector()
            self.coefficient = self.algebra.normalize_coefficient(coefficient)
        elif self.algebra.is_coefficient(exponents):
            self.exponents = self.algebra.empty_exponent_vector()
            self.coefficient = self.algebra.normalize_coefficient(exponents * coefficient)
        else:
            raise TypeError

    def exponent_tuple(self):
        return tuple(zip(self.exponents.keys(), self.exponents.values()))

    def to_polynomial(self):
        return self.algebra.PolynomialType(self.algebra, [self])


class GeneralPolynomial:
    AlgebraType = GeneralAlgebra

    def __init__(self, algebra, argument):
        if not isinstance(algebra, self.AlgebraType):
            raise TypeError
        self.algebra = algebra
        if isinstance(argument, list):
            self.monomials = []
            for monomial in argument:
                if not isinstance(monomial, self.algebra.MonomialType):
                    raise TypeError
                if monomial.algebra != algebra:
                    raise TypeError
            argument.sort(key=lambda monomial: monomial.exponent_tuple())
            i = 0
            while i < len(argument):
                exponent_tuple = argument[i].exponent_tuple()
                exponent_dict = argument[i].exponents
                coefficient_sum = 0
                while i < len(argument) and argument[i].exponent_tuple() == exponent_tuple:
                    coefficient_sum += argument[i].coefficient
                    i += 1
                if coefficient_sum != 0:
                    self.monomials.append(self.algebra.MonomialType(algebra, exponent_dict, coefficient_sum))
        elif isinstance(argument, type(self)):
            if self.algebra != argument.algebra:
                raise TypeError
            self.monomials = argument.monomials[:]
        elif self.algebra.is_coefficient(argument):
            if argument != 0:
                self.monomials = [self.algebra.MonomialType(algebra, argument)]
            else:
                self.monomials = []
        else:
            raise TypeError

    def __str__(self):
        if len(self.monomials) == 0:
            return '0'
        else:
            return ' + '.join(map(str, self.monomials))

    def str_parenthesis(self):
        if len(self.monomials) <= 1:
            return str(self)
        else:
            return f"({str(self)})"

    def __eq__(self, other):
        if self.algebra.is_coefficient(other):
            return self == type(self)(self.algebra, other)
        if isinstance(other, type(self)):
            if not self.algebra == other.algebra:
                raise TypeError
            exponents_eq = (tuple([m1.exponent_tuple() for m1 in self.monomials]) == tuple([m2.exponent_tuple() for m2 in other.monomials]))
            coefficients_eq = (tuple([m1.coefficient for m1 in self.monomials]) == tuple([m2.coefficient for m2 in other.monomials]))
            return exponents_eq and coefficients_eq
        return NotImplemented

    def __add__(self, other):
        if not self.algebra.is_element(other):
            return NotImplemented
        other = type(self)(self.algebra, other)
        return type(self)(self.algebra, self.monomials + other.monomials)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + other * (-1)

    def __rsub__(self, other):
        return self * (-1) + other

    def __mul__(self, other):
        if not self.algebra.is_element(other):
            return NotImplemented
        other = type(self)(self.algebra, other)
        product_monomials = []
        for m1 in self.monomials:
            for m2 in other.monomials:
                product_exponents = {var_id: m1.exponents[var_id] + m2.exponents[var_id] for var_id in self.algebra.get_all_ids()}
                product_coefficient = m1.coefficient * m2.coefficient
                product_monomial = self.algebra.MonomialType(self.algebra, product_exponents, product_coefficient)
                product_monomials.append(product_monomial)
        return type(self)(self.algebra, product_monomials)

    def __rmul__(self, other):
        return self * other

    def coefficient(self, monomial):
        if isinstance(monomial, self.algebra.MonomialType):
            if not monomial.algebra == self.algebra:
                raise TypeError
            for m in self.monomials:
                if m.exponents == monomial.exponents:
                    return m.coefficient
            return self.algebra.normalize_coefficient(0)
        elif isinstance(monomial, self.algebra.PolynomialType):
            if not monomial.algebra == self.algebra:
                raise TypeError
            if len(monomial.monomials) != 1:
                raise ValueError
            return self.coefficient(monomial.monomials[0])
        else:
            raise TypeError


GeneralAlgebra.CoefficientType = Fraction
GeneralAlgebra.MonomialType = GeneralMonomial
GeneralAlgebra.PolynomialType = GeneralPolynomial


class Algebra(GeneralAlgebra):
    def __init__(self, variable_identifiers):
        super().__init__(variable_identifiers)
    
    def normalize_exponent(self, exponent):
        if not isinstance(exponent, int):
            raise TypeError
        if exponent < 0:
            raise ValueError
        return exponent
    
class Monomial(GeneralMonomial):
    def __init__(self, algebra, exponents, coefficient=1):
        super().__init__(algebra, exponents, coefficient)

    def total_degree(self):
        return sum([self.exponents[var_id] for var_id in self.algebra.get_all_ids()])

    def __str__(self):
        factors = []
        if self.coefficient != 1 or self.total_degree() == 0:
            factors.append(str(self.coefficient))
        for var_id in self.algebra.get_all_ids():
            var_name = self.algebra.variables[var_id].representation
            degree = self.exponents[var_id]
            if degree > 0:
                if degree == 1:
                    factors.append(var_name)
                else:
                    factors.append(f"{var_name}^{degree}")
        return '*'.join(factors)
        

class Polynomial(GeneralPolynomial):
    def __init__(self, algebra, argument):
        super().__init__(algebra, argument)

    def max_degree(self):
        return max([m.total_degree() for m in self.monomials])


Algebra.CoefficientType = Fraction
Algebra.MonomialType = Monomial
Algebra.PolynomialType = Polynomial
