# External modules
import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot 	as plt
# Internal Modules
import plotClasses
from plotCommands 		import *

#############################################################################
#############################################################################

# INPUT
plots = [ZnPdconvEng,ZnPdEngAbsDiff,ZnPdConvEng_d]

#############################################################################
#############################################################################

def axMaker(n):
	divs  = [i for i in range(1,n+1) if n%i==0]
	nrows = divs[len(divs)/2]
	ncols = n / nrows
	f,axs = plt.subplots(nrows=nrows,ncols=ncols)
	if   n == 1: return [axs]
	elif n <= 3: return axs
	else: return [item for sublist in axs for item in sublist]  

def maker(pfdict):
	"""Take plotFunc dictionary from plotCommands. Feed into plotfuncmaker from plotClasses"""
	f = plotClasses.makerDict[pfdict.pop('maker')]
	return f.makeFunc(**pfdict)

#############################################################################
#############################################################################

def main():	
	axs = axMaker(len(plots))
	for i,p in enumerate(plots):  maker(p).plot(axs[i])
	plt.show()

if __name__ == '__main__': 	main()