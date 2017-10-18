from db import insertObject 

class Convergence(object):
	def __init__(self,dwrat,econv,mixing,nmix,maxstep,nbands,sigma,fmax,delta,climb):
		self.dwrat 		= dwrat
		self.econv 		= econv
		self.mixing	 	= mixing
		self.nmix 		= nmix
		self.maxstep 	= maxstep
		self.nbands 	= nbands
		self.sigma 		= sigma
		self.fmax 		= fmax
		self.delta 		= delta
		self.climb 		= climb
	
	def sqlTable(self):  return 'convergence'

	def sqlCols(self):   return ['dwrat'
								,'econv'
								,'mixing'
								,'nmix'
								,'maxstep'
								,'nbands'
								,'sigma'
								,'fmax'
								,'delta'
								,'climb']

	def sqlInsert(self): return [self.dwrat 
								,self.econv 	
								,self.mixing
								,self.nmix 		
								,self.maxstep	
								,self.nbands 	
								,self.sigma 	
								,self.fmax 		
								,self.delta	 	
								,self.climb]
	
	def sqlEq(self): return self.sqlCols()

conv1 = Convergence(10,1e-5,0.1,10,500,-12,0.1,0.05,0.03,False)  # BASELINE
conv2 = Convergence(10,1e-5,0.05,10,500,-12,0.2,0.05,0.03,False) # UP: Sigma, DOWN: Mixing
convs = [conv1,conv2]

def main():
	for conv in convs: insertObject(conv)

if __name__ == '__main__':
	main()
