import data_solids_wPBE 	as bData
import numpy 				as np
import matplotlib.pyplot 	as plt

from scipy.stats.mstats import gmean

from constraint import *
from db 		import plotQuery


identity = lambda x: x

"""
Functions for visualizing data from data.db (jobs and their results)

Experimental data functions depend currently on the structure of Keld's data_solids_wPBE

Constraints are used by plotting functions to select slices of the whole dataset

Two kinds of objects are introduced:
	- PlotFunc (effectively a relation between a job and a number (or a string, if used as a label function))
	- SummaryFunc (effectively a relation between a set of jobs and a number)

Ideally, elemental PlotFuncs and SummaryFuncs are defined in one module, another module has a bunch of interesting combinations, 
	and this module is modified to load different interesting queries and display to the screen
"""


####################################
# Experimental Data functions
# Designed only to parse Keld's data
####################################

def getExpt(inputFunc):
	def func(bulkname):
		stoich 	= bulkname.split('-')[0]
		lattice = bulkname.split('-')[1]
		try:
			key = ''.join(stoich)+'_'+lattice
			out = bData.data[key]
		except KeyError: #try other order
			key = ''.join(reversed(stoich))+'_'+lattice
			out = bData.data[key]
		return inputFunc(out)
	return func

def getExptA(dataEntry): 
	return dataEntry['lattice parameter']

def getExptBM(dataEntry): 
	try:             return dataEntry['bulk modulus']
	except KeyError: return dataEntry['bulk modulus kittel']

def getErrA(bulkname,a): 
	if 	'hcp' in bulkname: 	multFactor = 1
	elif 'bcc' in bulkname: multFactor = (3**(0.5)/2.)
	else:					multFactor = 2**(-0.5) 
	return a - multFactor*getExpt(getExptA)(bulkname)

def getErrBM(bulkname,bm): 
	print bulkname,bm,getExpt(getExptBM)(bulkname)
	return bm - getExpt(getExptBM)(bulkname)

##########################################
class PlotFunc(object):
	def __init__(self,name,cols,func,axisLabel = None,convTol=None,derivConvTol=None):
		if axisLabel is None: axisLabel = name
		self.name         = name         # String
		self.cols         = cols         # [String], JobID -> (Col#1,Col#2,...)     ] JobID -> Double
		self.func         = func         # (Col#1,Col#2,...) -> Double              ]
		self.axisLabel    = axisLabel    # String, default axis label
		self.convTol      = convTol      # Absolute magnitude difference from 'final' value considered converged
		self.derivConvTol = derivConvTol # {PlotFunc : (Double,Double)}


##########################################
##################
# Database Related
##################
exptA      = PlotFunc('Expt Lattice A, A', 		['bulkname'], 					lambda x: getExpt(getExptA)(x))
exptBM     = PlotFunc('Expt Bulk Mod, GPa',		['bulkname'], 					lambda x: getExpt(getExptBM)(x))
errA       = PlotFunc('Error in Lattice A, A', 	['bulkname','bulkresult.a'],	getErrA)
errBM      = PlotFunc('Error in BulkMod, GPa', 	['bulkname','bulkresult.bulkmodulus'],	getErrBM)
####################
# CALCULATOR RELATED
####################
pw          = PlotFunc('PW Cutoff, eV',           ['pw'],               identity)
kpt         = PlotFunc('KPT density, pt/A^-1',    ['kpt'],              identity)
pw15        = PlotFunc('PW Cutoff^(1.5)',         ['pw'],               lambda x: x**1.5)
################
# RESULT RELATED
################
energy      = PlotFunc('Energy, eV',              ['bulkresult.energy'],     identity)
calcBulkMod = PlotFunc('BulkModulus, GPa',        ['bulkresult.bulkmodulus'],identity)
bfit        = PlotFunc('Bulk Modulus Fit',        ['bulkresult.bfit'],       identity)            
calcA       = PlotFunc('Lattice parameter a, A',  ['bulkresult.a'],          identity) 
timePerStep = PlotFunc('Time Per Ionic Step, min',['bulkresult.time'],       identity)

################################
### LabelFunc :: (... -> String)
################################

labelName = PlotFunc('',['bulkname'],identity)
labelXC   = PlotFunc('',['xc'],identity)
labelCalc = PlotFunc('',['dftcode','xc','pw','kpt'],lambda x,y,z,q: '%s_%s_%d_%s'%(x,y,z,q))
##########################################################
class SummaryFunc(object):
	def __init__(self,name,func,kind):
		self.name 	= name       # String
		self.func 	= func       # PlotFunc 
		self.kind 	= kind       # [Double] -> Double
		"""NOT IMPLEMENTED: a general summaryfunc with a funcLIST, where KIND is [[Double]] -> Double"""
#########################################################

#########################
# Kinds of data summaries
#########################
def meanAbsolute(xs):  return sum(map(abs,xs))/float(len(xs))
def mean(xs):          return sum(xs)/float(len(xs))
def gMeanAbsolute(xs): return gmean(map(abs,xs))
def variance(xs):      raise NotImplementedError

###################
# Summary Functions
###################

maeLatParamA = SummaryFunc("MAE Lattice Constant 'a', A",errA,gMeanAbsolute)
maeBulkMod   = SummaryFunc("MAE Bulk Modulus, GPa",errBM,gMeanAbsolute)

########################################
def xyPlot(xPlotFunc
			,yPlotFunc
			,table
			,filterList
			,ax
			,legendLabel=None
			,labelFunc=None
			,title=None
			,commonConstraints=[completed]
			,codeList=['ro--','gx--','b^--','ms--','y*--','ko--','co--','ro:','gx:','b^:','ms:','y*:','ko:','co:','ro-','gx-','b^-','ms-','y*-','ko-','co-']):
	"""
	Plots a specified function (Job -> Double) vs another function, possibly over many sets of data (multiple lines) 
	"""
	xys=[]
	for i,constraintList in enumerate(filterList):	
		xs = [xPlotFunc.func(*x) for x in plotQuery(table,xPlotFunc.cols,constraintList+commonConstraints)]
		ys = [yPlotFunc.func(*y) for y in plotQuery(table,yPlotFunc.cols,constraintList+commonConstraints)]
		if labelFunc is not None: 
			label = [labelFunc.func(*l) for l in plotQuery(table,labelFunc.cols,constraintList+commonConstraints)]
		try: 
			xy= sorted(zip(xs,ys))  #order the pairs
			x,y = zip(*xy)
			ax.plot(x,y,codeList[i%len(codeList)],label='' if legendLabel is None else legendLabel[i])
			if labelFunc is not None: 
				for i in range(len(x)):	ax.annotate(label[i],xy=(x[i],y[i]),fontsize=9)
		except ValueError: print "Warning, no data found for constraint #"+str(i+1)
		xys.append(xy)
	if title is not None: ax.set_title(title)

	ax.set_xlabel(xPlotFunc.axisLabel)
	ax.set_ylabel(yPlotFunc.axisLabel)
	
	if legendLabel is not None: 
		legend = ax.legend(loc='best', shadow=True)
		legend.get_frame().set_facecolor('#00FFCC')
		legend.draggable()
	return xys

def barResults(summaryFunc,table,filterList,labelList,ax,commonConstraints=[completed]):
	"""
	Plots summaries of data as a bar graph
	"""
	y_pos   = np.arange(len(filterList)) # bar locations

	y = [summaryFunc.kind([summaryFunc.func.func(*x) for x in plotQuery(table,summaryFunc.func.cols,constraintList+commonConstraints)]) for constraintList in filterList] #this could've been clearer ....

	ax.bar(y_pos,y,align='center', alpha=0.5)
	ax.set_ylabel(summaryFunc.name)
	ax.set_xticks(y_pos);ax.set_xticklabels(labelList)
####################################################################################
####################################################################################
####################################################################################

def main():
	TWO = True
	
	if TWO:	fig, (ax1,ax2) 	= plt.subplots(nrows=1,ncols=2)
	else:	fig, (ax1) 		= plt.subplots(nrows=1,ncols=1)
	
	xyPlot(pw,timePerStep,'bulkresult',[[gpaw],[qe]],ax1
			#,singleElementLabels
			,title='Time per step vs PW'
			,commonConstraints=[li,kptden2,completed],labelFunc=None)
	if TWO:
		xyPlot(pw15,timePerStep,'bulkresult',[[gpaw],[qe]],ax2
				#,singleElementLabels
				,title='Time per step vs PW^1.5'
				,commonConstraints=[li,kptden2,completed],labelFunc=None)

	plt.show()


if __name__ == '__main__': 
	main()

########################
# Previous Plot Commands
#########################
"""
xyPlot(pw,errBM,'bulkresult',singleElementFilters,ax1
		,singleElementLabels
		,title='Gpaw'
		,commonFilters=[completed,rpbe,gpaw],labelFunc=labelCalc)

For some common calculator set (e.g. rbpe,gpaw), show B.M. error vs pw cutoff for all elements

		barResults(maeLatParamA,'bulkresult',[[gpaw],[qe]],['Gpaw ','qe'],ax2,[completed])

"""

##########################################################
##########################################################
##########################################################
"""
def isConverged(xFunc,yFunc,filters):
	
	#Uses three-point finite difference estimate of derivative: http://www.m-hikari.com/ijma/ijma-password-2009/ijma-password17-20-2009/bhadauriaIJMA17-20-2009.pdf
	#Tests whether or not yFunc derivative has a low magnitude and is decreasing over a certain range
	
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

"""
