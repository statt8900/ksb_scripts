import data_solids_wPBE as bData
import numpy as np
from scipy.stats.mstats import gmean

"""
This module contains PlotFunc/SummaryFunc classes and functions for plotting
"""
########################################################

bulkFilterFunc={'bulk':     lambda x: x.bulk.stoichStruct
				,'dft':     lambda x: x.dftCode
				,'xc':      lambda x: x.calc.xc
				,'pw':      lambda x: x.calc.pw
				,'kpt':     lambda x: x.calc.kpt
				,'eConv':   lambda x: x.calc.eConv
				,'mix':     lambda x: x.calc.mixing
				,'nmix':    lambda x: x.calc.nmix
				,'maxstep': lambda x: x.calc.maxsteps
				,'nbands':  lambda x: x.calc.nbands
				,'sigma':   lambda x: x.calc.sigma
				,'xtol':    lambda x: x.calc.xtol
				,'mag':     lambda x: x.calc.magmom
				,'psp':     lambda x: x.calc.psp}

surfFilterFunc= {'facet':  lambda x: x.surf.facet
				,'scale':  lambda x: x.surf.scale
				,'sym'  :  lambda x: x.surf.symmetric
				,'fixed':  lambda x: x.surf.fixed
				,'vac'  :  lambda x: x.surf.vacuum
				,'ads'  :  lambda x: x.surf.adsorbates}.update(bulkFilterFunc)

########################################################
########################################################
class PlotFunc(object):
	def __init__(self,name,func,axisLabel,convTol=None,derivConvTol=None):
		self.name         = name
		self.func         = func            # Job -> Double
		self.axisLabel    = axisLabel       # String
		self.convTol      = convTol  
		self.derivConvTol = derivConvTol # {PlotFunc : (Double,Double)}

# Label Functions
def getStoich(pureBulkJob):	 return pureBulkJob.bulk.stoichStruct[0]
def getPW(job):              return str(job.calc.pw)
def getKPT(job):             return str(job.calc.kpt)

# Pure data functions
def getBulkEnergy(bulkJob):	 return bulkJob.getBulkResult().energy
def getBulkModulus(bulkJob): return bulkJob.getBulkResult().bulkModulus
def getTime(bulkJob):        return bulkJob.getBulkResult().time
def getLatticeA(pureBulkJob):
	s    =  pureBulkJob.bulk.stoichStruct[1]
	rawA =  np.linalg.norm(pureBulkJob.getBulkResult().cell[0])
	if s in ['fcc','rocksalt','diamond','zincblende']: return rawA*2**0.5 # primitive != std unit cell
	else: return rawA

def isConverged(xFunc,yFunc,filters):
	"""
	Uses three-point finite difference estimate of derivative: http://www.m-hikari.com/ijma/ijma-password-2009/ijma-password17-20-2009/bhadauriaIJMA17-20-2009.pdf
	Tests whether or not yFunc derivative has a low magnitude and is decreasing over a certain range
	"""
	# Test whether all the infrastructure is in place to determine convergence	
	try:
		dydxMax = yFunc.derivConvTol[xFunc.name][0] # threshold for max |dy/dx| 
		xRange  = yFunc.derivConvTol[xFunc.name][1] # range (from last data point) over which |dy/dx| must be decreasing and beneath threshold
	except KeyError: 
		raise AttributeError, 'No derivative convergence parameters specified for deriving {0} w/r/t {1}'.format(yFunc.name,xFunc.name)

	xy=xyResults(xFunc,yFunc,filters)
	if len(x) < 4: raise ValueError, "Not enough data points to take a single derivative!"
	if len(x) != len(set(x)): raise ValueError, 'Filters poorly chosen; there are multiple datapoints with the same x value'
	
	for i in range(1,len(x)): #all values except for first and last
		h1,h2= float(x[i]-x[i-1]), float(x[i+1]-x[i])
		sumH,diffH,prodH,quotH = h1+h2,h1-h2,h1*h2,h1/h2
		f0,f1,f2=y[i-1],y[i],y[i+1]
		dydx.append(-1/(sumH*quotH)*f0 - diffH/(prodH)*f1 + quotH/sumH*f2)
	filtdYdX = [abs(dYdX) for (x_i,dYdX) in zip(x[1:-1],dydx) if abs(x_i-max(x)) < xRange]
	if len(filtdYdX < 2): filtdYdX = dydx[-2:] # we could still have convergence if only data points are outside minimum xRange -- it's just a harder test to pass
	monotonicallyDecreasing = all(x>y for x, y in zip(filtdXdY, filtdXdY[1:]))	
	absoluteTest = all(x < dydxMax for x in filtdXdY)
	return monotonicallyDecreasing and absoluteTest

def convergence(xFunc,yFunc,filters):
	assert isConverged(xFunc,yFunc,filters), 'Trying to find convergence on unconverged function'
	xy = xyResults(xFunc,yFunc,filters)
	a = xy[-1][1] # last y value, assumed appx equal to inf limit
	alpha = yFunc.convTol
	def converged(xyPairs): return all(abs(yo-a)<alpha for xo,yo in xyPairs)
	for xyPair in xy:
		if converged(xyPair): return xyPair[0] 

#######################################
# Experimental Data - related functions
# Really designed only to parse Keld's particular scripts with bulk data
#######################################
def jobToDataEntryName(pureBulkJob): 
	trueNames = [(name.split('_')[0],name) for name in bData.data.keys()]
	jobStoich = pureBulkJob.bulk.stoich().keys()
	name1 = ''.join(jobStoich); name2 = ''.join(jobStoich[::-1])
	filteredNames = [namePair[1] for namePair in trueNames if (name1 == namePair[0] or name2 == namePair[0])]
	if len(filteredNames)!=1: raise ValueError, 'Strategy for finding bulk data names broken, %d names match %s'%len(filteredNames,pureBulkJob)
	else: return filteredNames[0]
def getExptLattice(pureBulkJob): return bData.data[jobToDataEntryName(pureBulkJob)]['lattice parameter']
def getExptBulkMod(pureBulkJob): 
	try:             return bData.data[jobToDataEntryName(pureBulkJob)]['bulk modulus']
	except KeyError: return bData.data[jobToDataEntryName(pureBulkJob)]['bulk modulus kittel']
##############################
# Calculator related functions
##############################
def prodKPT(job): return job.calc.kpt[0]*job.calc.kpt[1]*job.calc.kpt[2]

calcLatParamA  = PlotFunc(	'calcLatParamA'
							,getLatticeA
							,"Calculated Lattice Constant 'a', A")

exptLatParamA  = PlotFunc(	'exptLatParamA'
							,getExptLattice
							,"Experimental Lattice Constant 'a', A")

errorLatParamA = PlotFunc(	'errorLatParamA'
							,lambda x: getLatticeA(x)-getExptLattice(x)
							,"Error in Lattice Constant 'a', A")

calcBulkMod    = PlotFunc(	'calcBulkMod'
							,getBulkModulus
							,"Calculated Bulk Modulus, GPa")

exptBulkMod    = PlotFunc(	'exptBulkMod'
							,getExptBulkMod
							,"Experimental Bulk Modulus, GPa")

errorBulkMod   = PlotFunc(	'errorBulkMod'
							,lambda x: getBulkModulus(x)-getExptBulkMod(x)
							,"Error in Bulk Modulus, GPa")



jobPW          = PlotFunc('jobPW',lambda x: x.calc.pw,'Plane-wave cutoff, eV')
jobKPT         = PlotFunc('jobKPT',prodKPT, 'Total K-points (x * y * z)')
jobTime        = PlotFunc('jobTime',getTime,'Simulation time, h')
##########################################################################
##########################################################################

class SummaryFunc(object):
	def __init__(self,func,kind,name):
		self.func      = func       # Job -> Double
		self.kind      = kind       # [Double] -> Double
		self.name      = name  # String

# Summarizing functions
def meanAbsolute(xs): return sum(map(abs,xs))/float(len(xs))
def mean(xs): return sum(xs)/float(len(xs))
def gMeanAbsolute(xs): return gmean(map(abs,xs))
def variance(xs): raise NotImplementedError

maeLatParamA = SummaryFunc(errorLatParamA,gMeanAbsolute,"MAE Lattice Constant 'a', A")
maeBulkMod   = SummaryFunc(errorBulkMod,gMeanAbsolute,"MAE Bulk Modulus, GPa")
##########################################################################
def filterComplete(filters):
	from batchjobs import getStatus
	result=[]
	for j in getStatus('complete'):
		filterResults = [True] if filters is None else [v(bulkFilterFunc[k](j)) for k,v in filters.items()]
		if filterResults.count(True)==len(filterResults): result.append(j)
	return result


##########################################################################
# PLOTTING
##########################################################################

def xyPlot(xPlotFunc,yPlotFunc,filterList,ax=None,codeList=None,labelFunc=None,title=None):
	"""
	For a given x and y function, plots over different domains (specified by filterList)
	are shown together
	"""
	for i,f in enumerate(filterList):	
		x,y,label=[],[],[]
		for job in filterComplete(f.filters):
			x.append(xPlotFunc.func(job))
			y.append(yPlotFunc.func(job))
			if labelFunc is not None: label.append(labelFunc(job))
		
		try: 
			xy= sorted(zip(x,y))  #order the pairs
			x,y = zip(*xy)
			if ax is not None: 
				ax.plot(x,y,codeList[i],label=f.name)
				if labelFunc is not None: 
					for i in range(len(x)):	ax.annotate(label,xy=(x,y),fontsize=9)

		except ValueError: print "Warning, no data found for "+f.name
	if title is not None: ax.set_title(title)
	ax.set_xlabel(xPlotFunc.axisLabel)
	ax.set_ylabel(yPlotFunc.axisLabel)
	legend = ax.legend(loc='best', shadow=True)
	legend.get_frame().set_facecolor('#00FFCC')
	legend.draggable()


def xyResult(xPlotFunc,yPlotFunc,filters=None):	
	"""
	Returns vector of x and y values for given x,y functions and a Filter object
	"""
	x,y=[],[]
	for j in filterComplete(filters.filters):
		x.append(xPlotFunc.func(j))
		y.append(yPlotFunc.func(j))
	
	xy = sorted(zip(x,y))  #order the pairs
	return zip(*xy)

def barResults(summaryFunc,filterList,ax):
	y_pos   = np.arange(len(filterList))
	objects = [f.name for f in filterList]
	filts   = [f.filters for f in filterList]
	y = [summaryFunc.kind(map(summaryFunc.func.func,filterComplete(f))) for f in filts]
	ax.bar(y_pos,y,align='center', alpha=0.5)
	ax.set_ylabel(summaryFunc.name)
	ax.set_xticks(y_pos);ax.set_xticklabels(objects)

def plotEq(xmin,xmax,ax):	ax.plot(np.linspace(xmin,xmax),np.linspace(xmin,xmax)) # plot stright line
