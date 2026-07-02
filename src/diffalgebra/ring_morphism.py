from typing import Optional

from .constant_ring import ConstantRing, ConstantGenerator, ConstantPolynomial, Constant, Term, QQ
from .diff_ring import DifferentialRing, FuncGenerator, DifferentialPolynomial, Expression, DiffTerm
from .exceptions import DefinitionError

class RingMorphism:
    _source: ConstantRing
    _target: ConstantRing
    _mapping: tuple[ConstantPolynomial]

    def __init__(self, source: ConstantRing, target: ConstantRing, mapping: dict[ConstantGenerator, Constant],
                 name: Optional[str] = None):
        self._source = source
        self._target = target
        _mapping = []
        for gen_name in source._gen_names:
            generator = source.gen(gen_name)
            if generator not in mapping:
                raise DefinitionError(f"Generator {gen_name} has no defined image")
            image = mapping[generator]
            if not target.is_element(image):
                raise TypeError(f"{image} is not an element of target {target}")
            _mapping.append(target.promote(image))
        for generator in mapping.keys():
            if not source.is_generator(generator):
                raise TypeError(f"{generator} is not a generator of {source}")
        
        self._mapping = tuple(_mapping)
        self._name = name


    def apply(self, expression: Constant) -> ConstantPolynomial:
        if not self._source.is_element(expression):
            raise TypeError(f"{expression} is not an element of {self._source}")
        image_terms: list[Term] = []
        expression = self._source.promote(expression)
        for monomial, coefficient in expression._terms:
            image_term = coefficient
            for gen_image, exp in zip(self._mapping, monomial):
                image_term *= (gen_image ** exp)
            image_terms.extend(self._target.promote(image_term)._terms)
        return ConstantPolynomial(ring=self._target, terms=image_terms)

    def __call__(self, argument: Constant) -> ConstantPolynomial:
        return self.apply(argument)
        
    @classmethod
    def identity(cls, ring: ConstantRing) -> RingMorphism:
        return RingMorphism(source=ring, target=ring, mapping={gen: gen for gen in ring._generators.values()})


class DiffRingMorphism:
    _source: DifferentialRing
    _target: DifferentialRing
    _base: RingMorphism
    _mapping: tuple[DifferentialPolynomial]

    def __init__(self, source: DifferentialRing, 
                 target: DifferentialRing,
                 mapping: dict[FuncGenerator, Expression],
                 base: Optional[RingMorphism] = None,
                 name: Optional[str] = None):
        if base is None:
            if source._base_ring == target._base_ring:
                base = RingMorphism.identity(source._base_ring)
            elif source._base_ring == QQ:
                base = RingMorphism(source=QQ, target=target._base_ring, mapping=dict())
            else:
                raise DefinitionError(f"Base ring homomorphism between base rings {source._base_ring} and {target._base_ring} must be provided")
        self._base = base

        self._source = source
        self._target = target
        _mapping = []
        for gen_name in source._func_names:
            generator = source.gen(gen_name)
            if generator not in mapping:
                raise DefinitionError(f"Generator {gen_name} has no defined image")
            image = mapping[generator]
            if not target.is_element(image):
                raise TypeError(f"{image} is not an element of target {target}")
            _mapping.append(target.promote(image))
        for generator in mapping.keys():
            if not source.is_generator(generator):
                raise TypeError(f"{generator} is not a generator of {source}")
        
        self._mapping = tuple(_mapping)
        self._name = name

    def apply(self, expression: Expression) -> DifferentialPolynomial:
        if not self._source.is_element(expression):
            raise TypeError(f"{expression} is not an element of {self._source}")
        image_terms: list[DiffTerm] = []
        expression = self._source.promote(expression)
        for term in expression._terms:
            monomial, coefficient = term
            image_term = self._base(coefficient)
            for gen_image, factors in zip(self._mapping, monomial):
                for derivative, exp in factors:
                    image_term *= gen_image.diff(order=derivative) ** exp
            image_terms.extend(self._target.promote(image_term)._terms)
        return DifferentialPolynomial(ring=self._target, terms=image_terms)
    
    def __call__(self, argument: Expression) -> DifferentialPolynomial:
        return self.apply(argument)