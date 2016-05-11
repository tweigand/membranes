from proteus import *
from proteus.default_p import *
from Membrane_Model import *
import numpy as np
from diffusion import *

name="Membrane_Model_Pores_Outflow"

L=(100.0,1.0,1.0)  
nd = 2
nc = 1

Length = 5000.

parallel = True

polyfile = "Domain"
#back out approximate length scale for mesh
#refinementLevel = 32
#he = L[0]/float(refinementLevel)
#if want to specify area constraint from command line
#use_uniform_refinement = False
#if use_uniform_refinement:
#    triangleOptions = "VApq30Dena%8.8f" % ((he**2)/2.0,)
#else:
#    triangleOptions = "VApq20Dena"
triangleOptions = "Apq30DenP" 
genMesh = True

T = 1.

Diff = {}

coefficients = Membrane_Model()

initialConditions = None

eps  = 1.e-6
def getDBC_Water(x,flag):
	if (flag == 2):
		return lambda x,t: 1.0-den_dry_mem/den_wet_mem
	if (flag == 3):
		return lambda x,t: (1.-den_dry_mem/den_wet_mem)*np.exp(-(Vol*dP)/(R_gas*Temp))


def noflux(x,flag):
	if (flag == 1):
		return lambda x,t: 0.0

dirichletConditions = {0:getDBC_Water}
advectiveFluxBoundaryConditions =  {0:noflux}
diffusiveFluxBoundaryConditions = {0:{0:noflux}}

import copy
class BoundaryDiffusiveFlux(AuxiliaryVariables.AV_base):
    def __init__(self):
        pass
    def attachModel(self,model,ar):
        """

        Get the model and setup the array F to store the flux on each part
        of the domain boundary that is shares a elementBoundaryFlag

        """
        self.model=model
        self.ar=ar
        self.writer = Archiver.XdmfWriter()
        self.nd = model.levelModelList[-1].nSpace_global
        m = self.model.levelModelList[-1]
        flagMax = max(m.mesh.elementBoundaryMaterialTypes)
        flagMin = min(m.mesh.elementBoundaryMaterialTypes)
        assert(flagMin == 0)
        assert(flagMax >= 0)
        self.nFluxes=flagMax+1
        self.levelFlist=[]
        for m in self.model.levelModelList:
            F = numpy.zeros((self.nFluxes,),'d')
            self.levelFlist.append(F)
        self.historyF=[]
        self.historyF.append(copy.deepcopy(self.levelFlist))
        return self
    def calculate(self):
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        for m,F in zip(self.model.levelModelList,self.levelFlist):
            F.flat[:]=0.0
            nQuadraturePoints_elementBoundary = m.ebqe[('u',0)].shape[1]
            for ebNE in range(m.mesh.nExteriorElementBoundaries_global):
                ebN = m.mesh.exteriorElementBoundariesArray[ebNE]                
                for k in range(nQuadraturePoints_elementBoundary):#number of
                    #compute -a.grad(u).n
                    diffusive_flux = 0.0
                    for I in range(m.nSpace_global): #loop over space dimension
                        for mm in range(m.coefficients.sdInfo[(0,0)][0][I],m.coefficients.sdInfo[(0,0)][0][I+1]): #loop over CSR representation of diffusion
                            J = m.coefficients.sdInfo[(0,0)][1][mm] #get column for nonzero entry
                            diffusive_flux -= m.ebqe[('a',0,0)][ebNE,k,mm]*m.ebqe[('grad(u)',0)][ebNE,k,J]*m.ebqe['n'][ebNE,k,I]
                        #                    
                    F[m.mesh.elementBoundaryMaterialTypes[ebN]] += diffusive_flux*m.ebqe[('dS_u',0)][ebNE,k]
            den_water = 1.e-21 #g/nm^3
            convert_nm_m = 1.e9
            if comm.rank==0:
            	# only processor 0 will actually get the data
            	totals = np.zeros_like(F)
            else:
            	totals = None
            comm.Reduce([F,MPI.DOUBLE], [totals, MPI.DOUBLE],op=MPI.SUM, root=0)
            if comm.rank == 0:
            	np.savetxt('test.txt', totals/Length/den_water/convert_nm_m)
            	log("Flux")
            	log(`totals/Length/den_water/convert_nm_m`)
        self.historyF.append(copy.deepcopy(self.levelFlist))




