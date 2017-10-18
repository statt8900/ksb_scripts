from ase.db import connect


def fMaybe(f,x): return f(x) if x is not None else x
def maybeStr(x): return fMaybe(str,x)
def maybeInt(x): return fMaybe(int,x)
def maybeFloat(x): return fMaybe(float,x)

class Result(object):
	def __init__(self,jobtableid,launchdir,aseid=None,energy=None,forces=None
				,bulkmodulus=None,bfit=None,xccoeffs=None,viblist=None
				,dos=None,barrier=None,time=None,niter=None):

		self.jobtableid 	= jobtableid 		#int
		self.launchdir 		= launchdir 		#string
		self.aseid 			= aseid 			# int
		self.energy 		= energy 			#float
		self.forces 		= forces 			#[[Float]]
		self.bulkmodulus 	= bulkmodulus 		#float
		self.bfit 			= bfit 				#float
		self.xccoeffs 		= xccoeffs 			#[[float]]
		self.viblist 		= viblist 			#[float]
		self.dos 			= dos 				# binary
		self.barrier 		= barrier 			#float
		self.time 			= time 				#float
		self.niter 			= niter 			#int
		
	def sqlTable(self):  return 'result'

	def sqlCols(self):   return ['jobtableid'
								,'launchdir'
								,'aseid'
								,'energy'
								,'forces'
								,'bulkmodulus'
								,'bfit'
								,'xccoeffs'
								,'viblist'
								,'dos'
								,'barrier'
								,'time'
								,'niter']

	def sqlInsert(self): return [self.jobtableid
								,self.launchdir
								,self.aseid
								,maybeFloat(self.energy)
								,maybeStr(self.forces)
								,maybeFloat(self.bulkmodulus)
								,maybeFloat(self.bfit)
								,maybeStr(self.xccoeffs)
								,maybeStr(self.viblist)
								,maybeStr(self.dos)
								,maybeFloat(self.barrier)
								,maybeFloat(self.time)
								,maybeInt(self.niter)]

	def sqlEq(self): return ['jobtableid']

