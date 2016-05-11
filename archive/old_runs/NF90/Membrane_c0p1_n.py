from proteus import *
from proteus.default_n import *
from Membrane_p import *
from Membrane_Model import *
from parameters import *
from diffusion import *

nDTout = 1

timeIntegration = NoIntegration

rtol_u[0] = 1.0e-3
atol_u[0] = 1.0e-3

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

massLumping = True

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton
maxNonlinearIts = 25
#maxLineSearches = 2
fullNewtonFlag = True
nl_atol_res = 1.0e-7
l_atol_res = 1.0e-10
tolFac = 1.0e-7


conservativeFlux = {0:'pwl-bdm'}

auxiliaryVariables = [BoundaryDiffusiveFlux()]