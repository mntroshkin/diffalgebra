from .PolynomialAlgebra import GeneralAlgebra, GeneralMonomial, GeneralPolynomial
from .DifferentialAlgebra import DiffAlgebra, DiffPolynomial


class ExteriorDiffAlgebra(GeneralAlgebra):
    _empty_exponent = tuple()
    _default_exponent = (0, )

    def __init__(self, variable_identifiers, coefficient_algebra=None):
        super().__init__(variable_identifiers)
        if not coefficient_algebra:
            coefficient_algebra = self.CoefficientAlgebraType([])
        if not isinstance(coefficient_algebra, self.CoefficientAlgebraType):
            raise TypeError
        self.coefficient_algebra = coefficient_algebra

    def is_coefficient(self, expression):
        return self.coefficient_algebra.is_element(expression)

    def normalize_coefficient(self, coefficient):
        return self.CoefficientType(self.coefficient_algebra, coefficient)

    def normalize_exponent_with_sign(self, exponent):
        if not (isinstance(exponent, list) or isinstance(exponent, tuple)):
            raise TypeError
        for diff_degree in exponent:
            if not isinstance(diff_degree, int):
                raise TypeError
            if diff_degree < 0:
                raise ValueError
        sorted_exponent, sign = signed_sort(exponent)
        sorted_exponent = tuple(sorted_exponent)
        return sorted_exponent, sign

    def normalize_exponent_vector_with_sign(self, exponent_vector):
        normalized_exponent_vector = self.empty_exponent_vector()
        total_sign = 1
        if not isinstance(exponent_vector, dict):
            raise TypeError
        for var_id in exponent_vector.keys():
            if not var_id in self.variables.keys():
                raise KeyError("Variable identifier not found in the algebra")
            normalized_exponent, sign = self.normalize_exponent_with_sign(exponent_vector[var_id])
            total_sign *= sign
            normalized_exponent_vector[var_id] = normalized_exponent
        if total_sign == 0:
            normalized_exponent_vector = self.empty_exponent_vector()
        return normalized_exponent_vector, total_sign

    def is_element(self, expression):
        if self.coefficient_algebra.is_element(expression):
            return True
        if isinstance(expression, ExteriorDiffPolynomial) and expression.algebra == self:
            return True
        return False


class ExteriorDiffMonomial(GeneralMonomial):
    AlgebraType = ExteriorDiffAlgebra

    def __init__(self, algebra, exponents, coefficient=1):
        super().__init__(algebra, exponents, coefficient)
        if not isinstance(algebra, ExteriorDiffAlgebra):
            raise TypeError
        if not algebra.coefficient_algebra.is_element(coefficient):
            raise TypeError
        self.algebra = algebra
        if isinstance(exponents, dict):
            if coefficient != 0:
                self.exponents, sign = self.algebra.normalize_exponent_vector_with_sign(exponents)
                coefficient *= sign
            else:
                self.exponents = self.algebra.empty_exponent_vector()
            self.coefficient = self.algebra.normalize_coefficient(coefficient)
        elif self.algebra.is_coefficient(exponents):
            self.exponents = self.algebra.empty_exponent_vector()
            self.coefficient = self.algebra.normalize_coefficient(exponents * coefficient)
        else:
            raise TypeError

    def __str__(self):
        factors = []
        if self.coefficient != 1 or self.max_diff_degree() == -1:
            factors.append(self.coefficient.str_parenthesis())
        for var_id in self.algebra.get_all_ids():
            var_name = self.algebra.variables[var_id].representation
            for p in range(0, self.max_diff_degree() + 1):
                p_count = self.exponents[var_id].count(p)
                if p_count > 0:
                    if p == 0:
                        factors.append(f"{var_name}")
                    else:
                        factors.append(f"{var_name}_{p}")
        return '*'.join(factors)

    def max_diff_degree(self):
        return max([max(self.exponents[var_id] + (-1,)) for var_id in self.algebra.get_all_ids()])


class ExteriorDiffPolynomial(GeneralPolynomial):
    def __init__(self, algebra, argument):
        super().__init__(algebra, argument)

    def __mul__(self, other):
        if not self.algebra.is_element(other):
            return NotImplemented
        other = type(self)(self.algebra, other)
        product_monomials = []
        for m1 in self.monomials:
            for m2 in other.monomials:
                product_exponents = dict()
                product_sign = 1
                swaps = 0
                var_number = len(self.algebra.get_all_ids())
                for i in range(var_number):
                    var_i = self.algebra.get_all_ids()[i]
                    for j in range(i + 1, var_number):
                        var_j = self.algebra.get_all_ids()[j]
                        swaps += len(m2.exponents[var_i]) * len(m1.exponents[var_j])
                    for j in range(i):
                        var_j = self.algebra.get_all_ids()[j]
                        swaps += len(m2.exponents[var_i]) * len(m2.exponents[var_j])
                sign = (-1) ** swaps
                product_exponents = {var_id: m1.exponents[var_id] + m2.exponents[var_id] for var_id in
                                     self.algebra.get_all_ids()}
                product_coefficient = m1.coefficient * m2.coefficient * sign
                product_monomial = self.algebra.MonomialType(self.algebra, product_exponents, product_coefficient)
                product_monomials.append(product_monomial)
        return type(self)(self.algebra, product_monomials)

    def max_diff_degree(self):
        return max([m.max_diff_degree() for m in self.monomials])

    def diff_x(self, n=1):
        if n == 0:
            return self
        if n == 1:
            summands = []
            for monomial in self.monomials:
                coefficient = monomial.coefficient.diff_x()
                exponents = monomial.exponents
                new_monomial = ExteriorDiffMonomial(self.algebra, exponents, coefficient)
                summands.append(ExteriorDiffPolynomial(self.algebra, [new_monomial]))
            for var_id in self.algebra.get_all_ids():
                for i in range(0, self.max_diff_degree() + 1):
                    dvar_i = self.algebra.get_monomial({var_id: (i + 1,)})
                    summands.append(dvar_i * self.diff_wrt_variable(var_id, i))
            return sum(summands)
        else:
            return self.diff_x(n - 1).diff_x()

    def diff_wrt_variable(self, var_id, degree):
        summands = []
        var_id_index = self.algebra.get_all_ids().index(var_id)
        for monomial in self.monomials:
            preceding_count = sum([len(monomial.exponents[var_id_other]) for var_id_other in self.algebra.get_all_ids()[:var_id_index]])
            var_count = monomial.exponents[var_id].count(degree)
            if var_count != 0:
                var_index = monomial.exponents[var_id].index(degree)
                new_coefficient = var_count * monomial.coefficient * (-1) ** (var_index + preceding_count)
                var_exponents = monomial.exponents[var_id][:var_index] + monomial.exponents[var_id][var_index + 1:]
                new_exponents = {var_id2: (monomial.exponents[var_id2] if var_id2 != var_id else var_exponents) for var_id2 in
                               self.algebra.get_all_ids()}
                summands.append(ExteriorDiffMonomial(self.algebra, new_exponents, new_coefficient))
        return ExteriorDiffPolynomial(self.algebra, summands)

    def variational_derivative(self, var_id):
        result = ExteriorDiffPolynomial(self.algebra, 0)
        for i in range(0, self.max_diff_degree() + 1):
            result += (-1) ** i + self.diff_wrt_variable(var_id, i).diff_x(i)
        return result


def signed_sort(lst):
    lst = list(lst)
    sign = 1
    for i in range(len(lst)):
        for j in range(i + 1, len(lst)):
            if lst[i] > lst[j]:
                sign *= -1
                lst[i], lst[j] = lst[j], lst[i]
    for i in range(len(lst) - 1):
        if lst[i] == lst[i + 1]:
            sign = 0
    return tuple(lst), sign


ExteriorDiffAlgebra.CoefficientAlgebraType = DiffAlgebra
ExteriorDiffAlgebra.CoefficientType = DiffPolynomial
ExteriorDiffAlgebra.MonomialType = ExteriorDiffMonomial
ExteriorDiffAlgebra.PolynomialType = ExteriorDiffPolynomial