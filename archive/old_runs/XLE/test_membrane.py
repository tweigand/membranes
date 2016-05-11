class Mem_Proteus:
	from proteus.iproteus import *
	from proteus import Comm
	from proteus import Archiver
	import numpy 
	import Membrane_p
	import Membrane_c0p1_n
	import parameters
	import Membrane_Model
	import diffusion


	def __init__(self):
		pass 


	def proteus_eval(self,k):

		f = open('diffusion.py','w')
		f.write("D_ws = ")
		f.write("%f" %k)
		f.close()

		reload(self.diffusion)
		reload(self.parameters)
		reload(self.Membrane_Model)
		reload(self.Membrane_p)
		reload(self.Membrane_c0p1_n)
		self.Comm.init()

		self.Profiling.logLevel=0
		self.Profiling.verbose=False
		self.pList = [self.Membrane_p]
		self.nList = [self.Membrane_c0p1_n]    
		self.so = self.default_so
		self.so.name = self.pList[0].name = "membrane_test"
		self.so.sList=[self.default_s]
		self.opts.logLevel=7
		self.opts.verbose=False
		self.opts.profile=True
		self.ns = self.NumericalSolution.NS_base(self.so,self.pList,self.nList,self.so.sList,self.opts)
		self.ns.calculateSolution('Membrane_test')

		flux_proteus = self.numpy.loadtxt('test.txt')

		print flux_proteus[3], k, "Values"
		flux_lab = 0.0000131544901065449
		print (flux_lab-flux_proteus[3])*(flux_lab-flux_proteus[3])
		return (flux_lab-flux_proteus[3])*(flux_lab-flux_proteus[3])