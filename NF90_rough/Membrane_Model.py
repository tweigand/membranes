from math import *
import numpy
import proteus
from proteus.TransportCoefficients import *
from parameters import *
from diffusion import *

class Membrane_Model(TC_base):
    """
    (Mu)_t + \deld (Bu - A \grad phi) + C u = 0 
    """
    
    def __init__(self,nc=1,nd=2):
        self.Material_id = {}
        self.Material_id_qe = {}
        self.Material_id_ebq = {}
        self.Material_id_ebqe = {}
        self.Diff_w = {}
        self.Diff_qe_w = {}
        self.Diff_ebq_e_w = {}
        self.Diff_ebqe_e_w = {}
        mass={}
        advection={}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={}
        self.nd=nd
        mass = {0:{0:'constant'}}
        advection = {0:{0:'constant'}}
        diffusion = {0:{0:{0:'constant'}}}
        potential = {0:{0:'u'}}
        reaction = {0:{0:'constant'}}
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian)
        self.variableNames=['Water_Mass_Fraction']


    def initializeMesh(self,mesh):
        for eN in range(mesh.nElements_global):
             materialFlag = mesh.elementMaterialTypes[eN]
             #print materialFlag,eN
        self.elementMaterialTypes = mesh.elementMaterialTypes
        self.exteriorElementBoundaryTypes =  numpy.zeros((mesh.nExteriorElementBoundaries_global),'i')
        for ebNE in range(mesh.nExteriorElementBoundaries_global):
            ebN = mesh.exteriorElementBoundariesArray[ebNE]
            eN  = mesh.elementBoundaryElementsArray[ebN,0]
            self.exteriorElementBoundaryTypes[ebNE] = self.elementMaterialTypes[eN]
        self.elementBoundaryTypes = numpy.zeros((mesh.nElementBoundaries_global,2),'i')
        self.elementBoundariesArray = mesh.elementBoundariesArray
        for ebN in range(mesh.nElementBoundaries_global):
            eN_left = mesh.elementBoundaryElementsArray[ebN,0]
            eN_right= mesh.elementBoundaryElementsArray[ebN,1]
            self.elementBoundaryTypes[ebN,0] = self.elementMaterialTypes[eN_left]
            if eN_right >= 0:
                self.elementBoundaryTypes[ebN,1] = self.elementMaterialTypes[eN_right]
            else:
                self.elementBoundaryTypes[ebN,1] = self.elementMaterialTypes[eN_left]

    def initializeElementQuadrature(self,t,cq):
        self.q_e = cq[('u',0)]
        count = 0
        for ci in range(self.nc):
            cq[('f',ci)].flat[:] = 0.0
            for eN in range(cq['x'].shape[0]):
                material=self.elementMaterialTypes[eN]
                for k in range(cq['x'].shape[1]):
                  if (material == 100):
                    self.Material_id_qe[count] = material
                    self.Diff_qe_w[count] = D_ww
                  else:
                    self.Material_id_qe[count] = material
                    self.Diff_qe_w[count] = D_ws
                  count = count + 1

 
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        self.ebq_e = cebq[('u',0)]
        count = 0
        nd = self.nd
        for ci in range(self.nc):
            for ebN in range(cebq_global['x'].shape[0]):
                material_left = self.elementBoundaryTypes[ebN,0]
                material_right= self.elementBoundaryTypes[ebN,1]
                if (material_left == 100):
                  D_ml_w = D_ww
                else:
                  D_ml_w = D_ws
                if (material_right == 100):
                  D_mr_w = D_ww
                else:
                  D_mr_w = D_ws
                for k in range(cebq_global['x'].shape[1]):
                    self.Material_id_ebq[count] = material_right
                    self.Diff_ebq_e_w[count] = 2.0/( (1./D_ml_w) + (1./D_mr_w) )
                    count = count + 1
            for eN in range(cebq['x'].shape[0]):
                for ebN_local in range(cebq['x'].shape[1]):
                    ebN = self.elementBoundariesArray[eN,ebN_local]
                    material_left = self.elementBoundaryTypes[ebN,0]
                    material_right= self.elementBoundaryTypes[ebN,1]
                    if (material_left == 100):
                       D_ml_w = D_ww
                    else:
                       D_ml_w = D_ws
                    if (material_right == 100):
                       D_mr_w = D_ww
                    else:
                       D_mr_w = D_ws
                    for k in range(cebq['x'].shape[2]):
                       self.Material_id_ebq[count] = material_right
                       self.Diff_ebq_e_w[count] = 2.0/( (1./D_ml_w) + (1./D_mr_w) )
                       count = count + 1


    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        self.ebqe_e = cebqe[('u',0)]
        count = 0
        nd = self.nd
        for ci in range(self.nc):
            for ebNE in range(cebqe['x'].shape[0]):
                material = self.exteriorElementBoundaryTypes[ebNE]
                for k in range(cebqe['x'].shape[1]):
                  if (material == 100):
                    self.Material_id_ebqe[count] = material
                    self.Diff_ebqe_e_w[count] = D_ww
                  else:
                    self.Material_id_ebqe[count] = material
                    self.Diff_ebqe_e_w[count] = D_ws
                  count = count + 1                  

    def evaluate(self,t,c):

        if c[('u',0)].shape == self.q_e.shape:
           self.Material_id = self.Material_id_qe
           self.Diff_w = self.Diff_qe_w
        elif c[('u',0)].shape == self.ebq_e.shape:
            self.Material_id = self.Material_id_ebq
            self.Diff_w = self.Diff_ebq_e_w
        elif c[('u',0)].shape == self.ebqe_e.shape:
           self.Material_id = self.Material_id_ebqe
           self.Diff_w = self.Diff_ebqe_e_w





        space_dim = c[('f',0)].shape[-1]

        for k in range(len(c[('u',0)].flat)):
            w_w=c[('u',0)].flat[k]

            c[('m',0)].flat[k] = 0.0
            c[('dm',0,0)].flat[k] = 0.0

            for j in range(space_dim):
                c[('a',0,0)].flat[k*space_dim**2+j*space_dim+j] = den_wet_mem*self.Diff_w[k]
                c[('da',0,0,0)].flat[k*space_dim**2+j*space_dim+j] = 0.0



