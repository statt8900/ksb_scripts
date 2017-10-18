import data_solids_wPBE 	as bData
import numpy 				as np
import matplotlib.pyplot 	as plt

from scipy.stats.mstats import gmean

from constraint import *
from details 	import *
from db 		import plotQuery


identity = lambda x: x

"""
Functions for visualizing data from data.db (jobs and their results)

Constraints are used by plotting functions to select slices of the whole dataset

SummaryFunc (effectively a relation between a set of jobs and a number)

"""

##########################################################
class SummaryFunc(object):
	def __init__(self,name,details,kind):
		self.name 	 = name       # String
		self.details = details    # [Detail] 
		self.kind 	 = kind       # [[Double]] -> Double

################################################################
# Kinds of data summaries (on one list) :: [[Double]] -> Double
################################################################
def meanAbsolute(xs):  return sum(map(abs,xs[0]))/float(len(xs[0]))
def mean(xs):          return sum(xs[0])/float(len(xs[0]))
def gMeanAbsolute(xs): return gmean(map(abs,xs[0]))
def variance(xs):      raise NotImplementedError

###################
# Summary Functions
###################

maeLatParamA = SummaryFunc(r"MAE Lattice Constant 'a', $\AA$",[errA],gMeanAbsolute)
maeBulkMod   = SummaryFunc(r"MAE Bulk Modulus, GPa",[errBM],gMeanAbsolute)

########################################
def xyPlot(xDetail
			,yDetail
			,filterList
			,ax
			,legendLabel= None
			,labelDetail= None
			,title 		= None
			,commonConstraints=[COMPLETED]
			,codeList=['ro--','gx--','b^--','ms--','y*--','ko--','co--','ro:','gx:','b^:','ms:','y*:','ko:','co:','ro-','gx-','b^-','ms-','y*-','ko-','co-']):
	"""
	Plots a specified function (Job -> Double) vs another function, possibly over many sets of data (multiple lines) 
	"""
	xys=[]
	for i,constraintList in enumerate(filterList):	
		xs = [x[0] for x in plotQuery([xDetail.name],constraintList+commonConstraints)]
		ys = [y[0] for y in plotQuery([yDetail.name],constraintList+commonConstraints)]

		if labelDetail is not None: 
			label = [l[0] for l in plotQuery([labelDetail.name],constraintList+commonConstraints)]
		try: 
			xy= sorted(zip(xs,ys))  #order the pairs
			x,y = zip(*xy)
			ax.plot(x,y,codeList[i%len(codeList)],label='' if legendLabel is None else legendLabel[i])
			if labelDetail is not None: 
				for ii in range(len(x)):	ax.annotate(label[ii],xy=(x[ii],y[ii]),fontsize=9)
		except ValueError: print "Warning, no data found for constraint %s"%([str(x) for x in constraintList])
		xys.append(xy)

	if title is not None: ax.set_title(title)

	ax.set_xlabel(xDetail.axislabel)
	ax.set_ylabel(yDetail.axislabel)
	
	if legendLabel is not None: 
		legend = ax.legend(loc='best', shadow=True)
		legend.get_frame().set_facecolor('#00FFCC')
		legend.draggable()
	return xys

def barResults(summaryFunc,filterList,labelList,ax,commonConstraints=[COMPLETED]):
	#Plots summaries of data as a bar graph
	y_pos   = np.arange(len(filterList)) # bar locations
	y = [summaryFunc.kind(plotQuery([d.name for d in summaryFunc.details],constraintList+commonConstraints)) for constraintList in filterList]

	ax.bar(y_pos,y,align='center', alpha=0.5)
	ax.set_ylabel(summaryFunc.name)
	ax.set_xticks(y_pos);ax.set_xticklabels(labelList)

####################################################################################
####################################################################################

def main():
	TWO = False
	
	if TWO:	fig, (ax1,ax2) 	= plt.subplots(nrows=1,ncols=2)
	else:	fig, (ax1) 		= plt.subplots(nrows=1,ncols=1)
	
	xyPlot(pw,errA,[[PD,GPAW],[PD,QE],[AU,GPAW],[AU,QE]],ax1
			,["Pd,GPAW",'Pd QE',"Au,GPAW",'Au QE']
			,title=''
			,commonConstraints=[BEEF],labelDetail=None)

	if TWO:
		barResults(maeLatParamA,[[GPAW],[QE]],['Gpaw ','Q.E.'],ax2,[COMPLETED])
	
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
