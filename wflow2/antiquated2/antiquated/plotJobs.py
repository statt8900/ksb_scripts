# External modules
import matplotlib,itertools,thread
import scipy.optimize 		as opt
import scipy.stats.mstats 	as mstatt
import numpy 				as np
matplotlib.use('GTKAgg')
import matplotlib.pyplot 	as plt
# Internal Modules
import constraint
import details as d
import printParse 		as pp 
import dbase 			as db
from plotCommands 		import *

identity = lambda x: x

"""
Functions for visualizing data from data.db (jobs and their results)

"""

##########################################################


class PlotFuncMaker(object):
	def __init__(self,func,xlabelMaker,ylabelMaker,nameMaker):
		self.func 		= func
		self.xlabelMaker= xlabelMaker
		self.ylabelMaker= ylabelMaker
		self.nameMaker	= nameMaker

	def makePlotFunc(self,cols): return PlotFunc(self.makeName(cols),cols,self.makeX(cols),self.makeY(cols),self.func)
	def makeX(self,cols): 		return self.xlabelMaker(cols)
	def makeY(self,cols): 		return self.ylabelMaker(cols)
	def makeName(self,cols): 	return self.nameMaker(cols)

basicFunc =	PlotFuncMaker(lambda xlist,ylist: [(x,y) for x,y in zip(xlist,ylist)]
						,lambda cs: cs[0].colname
						,lambda cs: cs[1].colname
						,lambda cs: '_vs_'.join([c.colname for c in cs]))

def absDiffLast(xlist,ylist): 
		if ylist.count(None)==len(ylist): return zip(xlist,ylist)
		for x,y in sorted(zip(xlist,ylist)):
			if y is not None: maxy = y
		
		return zip(xlist,[None if y is None else abs(y-maxy) for y in ylist])

absDiffFunc=PlotFuncMaker(absDiffLast
						,lambda cs: cs[0].colname
						,lambda cs: 'Convergence of '+cs[1].colname
						,lambda cs: 'absdiff '+'_vs_'.join([c.colname for c in cs]))

class PlotFunc(object):
	def __init__(self,name,cols,xlabel,ylabel,func):
		self.name=name 							# STRING
		self.cols=[c.colname for c in cols]		# [STRING]
		self.xlabel=xlabel;self.ylabel=ylabel 	# STRING
		self.func=func 							# [[a]] -> [(FLOAT,FLOAT)]
	def addAxes(self,ax): ax.set_ylabel(self.ylabel);ax.set_xlabel(self.xlabel)

	def makeTriples(self,idgroups,labelName): 
		args 		= [db.query(self.cols,FWIDS(x)) for x in idgroups] #[[(ARG1,ARG2,...)]]
		labels 		= [db.queryCol(labelName,FWIDS(x)) for x in idgroups] #[[STRING]]
		xy 			= [self.func(*zip(*arg)) for arg in args]   #[[(FLOAT,FLOAT))
		xy_l 		= [(zip(xy,l)) for xy,l in zip(xy,labels)]
		filtered 	= [[(xyl[0][0],xyl[0][1],xyl[1]) for xyl in xyls if not (xyl[0][0] is None or xyl[0][1] is None)] for xyls in xy_l]
		sortd 		= [sorted(x) for x in filtered if len(x)>0]
		return [zip(*x) for x in sortd]


def eqPlot(	ax
			,plotFuncName
			,eqclasses
			,labelName='blank'
			,const='1'
			,legName=None
			,title=None
			,allowSingletons=True
			,fit=False
			,codeList=['ro--','gx--','b^--','ms--','y*--','ko--','co--','ro:','gx:','b^:','ms:','y*:','ko:','co:','ro-','gx-','b^-','ms-','y*-','ko-','co-']):
	if title is not None: ax.set_title(title)
	plotFunc = pDict[plotFuncName]
	plotFunc.addAxes(ax)

	classes = [eq.classes(extraConstraint=const) for eq in eqclasses] #[[[Int]]]
	allids 	= misc.allElems(eqclasses[0].output(extraConstraint=const))
	fingerp = {ID:tuple([misc.getFingerprint(ID,clss) for clss in classes]) for ID in allids}
	finger 	= misc.invertDict(fingerp)
	ids 	= finger.values()

	if legName is not None: legendLabel = [db.query1(legName,'fwid',i[0],'job') for i in ids] 											# ideally, legFunc is identical for all elements in class

	xyls = plotFunc.makeTriples(ids,labelName)

	for i,(xs,ys,ls) in enumerate(xyls):
		if allowSingletons or len(xs)>1:
			print len(xs),ls[0]
			if fit:
				xs=np.array(xs);ys = np.array(ys)
				def mdl(a): return a[0]*xs+a[1]
				def errVec(a): return mdl(a) - ys 
				try: 
					aHat, success = opt.leastsq(errVec, [1,0])
					yhat 	= mdl(aHat)
					r2 		= 1 - (np.sum((yhat-ys)**2) / float(np.sum((ys - np.sum(ys)/len(ys))**2))) # 1 - SSR / SST, SST = sum_i (yi-ybar)^2
					xfit 	= np.linspace(min(xs),max(xs),num=len(xs));yfit=np.linspace(min(yhat),max(yhat),num=len(xs)) #aHat*xfit
					newleg 	= '' if legName is None else str(legendLabel[i])+', fitted, r2=%f'%r2
					fitline = ax.plot(xfit,yfit,codeList[i%len(codeList)][:2]+':',linewidth=0.5,markersize=0.3,label=newleg)
				except: pass

			ax.plot(xs,ys,codeList[i%len(codeList)],label='' if legName is None else legendLabel[i])
			for ii in range(len(xs)): ax.text(xs[ii],ys[ii],ls[ii],size=9)
		
	try:
		legend = ax.legend(loc='upper left', shadow=True,fontsize='small')
		legend.get_frame().set_facecolor('#00FFCC')
		legend.draggable()

	except AttributeError: pass





###################
# Summary Functions
###################

maeLatParamA = BarFunc(r"MAE Lattice Constant 'a', $\AA$",['err_A'],gMeanAbsolute)
maeBulkMod   = BarFunc(r"MAE Bulk Modulus, GPa",['errBM'],gMeanAbsolute)

########################################

def barResults(ax,summaryFunc,labelFunc,eq,const='1'):
	#Plots summaries of data as a bar graph
	ids 	= [list(x) for x in eq.classes(extraConstraint=const)]
	xs  	= [db.query1(labelFunc,'fwid',i[0],'job') for i in ids] # [String]
	raw_ys  = [db.query(summaryFunc.details,FWIDS(i)) for i in ids] # [[a,b,c,...]]
	col_ys 	= zip(summaryFunc.details * len(raw_ys),raw_ys)
	ys 		= [summaryFunc.kind(summaryFunc.details,args) for cols,args in col_ys] 			# [Float]

	y_pos   = np.arange(len(ids)) # bar locations


	ax.bar(y_pos,ys,align='center', alpha=0.5)
	ax.set_ylabel(summaryFunc.name)
	ax.set_xticks(y_pos);ax.set_xticklabels(xs)

####################################################################################
####################################################################################

def plotDict(x,ax):
	keys=x.keys()
	if 'plotFunc' in keys:
		eqPlot(ax,x['plotFunc'],x['eq'],x['labelName'] if 'labelName' in keys else 'blank'
					,x['const'] if 'const' in keys else '1'
					,x['legName'] if 'legName' in keys else None
					,x['title'] if 'title' in keys else ''
					,x['singletons'] if 'singletons' in keys else True
					,x['fit'] if 'fit' in keys else False)
	elif 'summaryFunc' in keys:
		barResults(ax,x['summaryFunc'],x['labelFunc'],x['eq'],x['const'] if 'const' in keys else '1')

	else: raise ValueError

binaryFuncs = [basicFunc,absDiffFunc]
plotCols 	= [c for c in d.allDetails if c.axislabel is not None]
plotFuncs 	= [x.makePlotFunc([y,z]) for x,y,z in itertools.product(binaryFuncs,plotCols,plotCols)]
pDict 		= {pf.name:pf for pf in plotFuncs}

def main():	
	print "STARTING PLOT:"

	load = True ### Only do this once (if database doesn't change
	if load:
		d.load(d.staticcol,incomplete=False,pairs=True)

	fig, (ax) 	= plt.subplots(nrows=1,ncols=1)
	#(fig,ax) = plt.subplots(nrows=1,ncols=1)
	plotDict(interstitialEnergy,ax)

	plt.show()


if __name__ == '__main__': 	main()


