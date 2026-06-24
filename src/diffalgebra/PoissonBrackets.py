from .DifferentialAlgebra import DiffAlgebra
from .ExteriorDifferentialAlgebra import ExteriorDiffAlgebra


class MultivectorAlgebra(ExteriorDiffAlgebra):
    def __init__(self, diff_algebra: DiffAlgebra):
        variable_identifiers = ["theta{" + var_id + "}" for var_id in diff_algebra.get_all_ids()]
        super().__init__(variable_identifiers, diff_algebra)
