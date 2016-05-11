from proteus import *
from proteus.default_n import *
from Membrane_p import *
from Membrane_Model import *
from parameters import *
from diffusion import *

nDTout = 1


timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

matrix = SparseMatrix

massLumping = True


#class MyDummy_exterior_Strong(ConstantAdvection_Diffusion_SIPG_exterior):
#	#useStrongDirichletConstraints=True
#	penalty_constant = 100000.0


class MyDummy_exterior_Strong(ConstantAdvection_Diffusion_IIPG_exterio):
    def __init__(self,
                 vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        Advection_DiagonalUpwind_Diffusion_IIPG_exterior.__init__(self,
                                                                  vt,
                                                                  getPointwiseBoundaryConditions,
                                                                  getAdvectiveFluxBoundaryConditions,
                                                                  getDiffusiveFluxBoundaryConditions)
        self.penalty_constant = 100000.0
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=-1.0


if parallel:
    multilevelLinearSolver = KSP_petsc4py
    #for petsc do things lie
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    levelLinearSolver = KSP_petsc4py#
    #for petsc do things like
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    #levelLinearSolver = PETSc#
    #pick number of layers to use in overlap
    nLayersOfOverlapForParallel = 0
    #type of partition
    parallelPartitioningType = MeshParallelPartitioningTypes.node
    #parallelPartitioningType = MeshParallelPartitioningTypes.element
    #have to have a numerical flux in parallel
    numericalFluxType = MyDummy_exterior_Strong#ConstantAdvection_Diffusion_SIPG_exterior#Advection_DiagonalUpwind_Diffusion_IIPG_exterior
    #for true residual test
    linearSolverConvergenceTest = 'r-true'
    #to allow multiple models to set different ksp options
    #linear_solver_options_prefix = 'poisson_'
    linearSmoother = None
else:
	pass    
multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton
maxNonlinearIts = 1
fullNewtonFlag = False


nl_atol_res = 1.0e-7
l_atol_res = 1.0e-20
tolFac = 0.0


conservativeFlux = {0:'pwl-bdm'}

auxiliaryVariables = [BoundaryDiffusiveFlux()]