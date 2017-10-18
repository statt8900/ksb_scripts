#External Modules
import numpy as np
# Internal Modules
import dbase,constraint
import details as d
import dbase 			as db


########################################
true 	= lambda x: True
fst 	= lambda x: x[0]; mapfst = lambda xs: map(fst,xs)
snd 	= lambda x: x[1]; mapsnd = lambda xs: map(snd,xs)

codeList=['ro--','gx--','b^--','ms--','y*--','ko--','co--','ro:','gx:','b^:','ms:','y*:','ko:','co:','ro-','gx-','b^-','ms-','y*-','ko-','co-']

##############################################
class PlotFuncMaker(object):
	def __init__(self,name,xyfunc,colMaker,ordMaker,xlabelMaker,ylabelMaker,lFunc,legFunc,titleMaker):
		self.name 		= name 			# corresponds to 'maker' in plotCommands
		self.xyfunc 	= xyfunc 		# [[a]] -> [(FLOAT,FLOAT)]
		self.colMaker 	= colMaker	 	# String -> [String]
		self.ordMaker 	= ordMaker
		self.xlabelMaker= xlabelMaker
		self.ylabelMaker= ylabelMaker
		self.lFunc 		= lFunc
		self.legFunc 	= legFunc
		self.titleMaker = titleMaker

	def makeXlabel(self,cols): 		return self.xlabelMaker(cols)
	def makeYlabel(self,cols): 		return self.ylabelMaker(cols)
	def makeFunc(self
					,colstr
					,cnst 	= '1'
					,gb 	= 'xc'
					,h 		= true
					,lname 	= 'blank'
					,legname= 'blank'): 
		cols = colstr.split()
		return PlotFunc( cols 		= self.colMaker(cols,lname,legname)
						,title 		= self.titleMaker(cols)
						,xLabel 	= self.xlabelMaker(cols)
						,yLabel 	= self.ylabelMaker(cols)
						,constraint = cnst
						,groupby	= gb.split()
						,having		= h
						,order		= self.ordMaker(cols) #order by x
						,xyFunc 	= self.xyfunc
						,lFunc 		= self.lFunc
						,legFunc	= self.legFunc)
###########################################################################
class BarFuncMaker(object):
	def __init__(self,name,xfunc,yfunc,ordMaker,xlabelMaker,ylabelMaker,titleMaker):
		self.name 		= name
		self.xfunc 		= xfunc
		self.yfunc 		= yfunc
		self.ordMaker 	= ordMaker
		self.xlabelMaker= xlabelMaker
		self.ylabelMaker= ylabelMaker
		self.titleMaker = titleMaker

	def makeXlabel(self,cols): 		return self.xlabelMaker(cols)
	def makeYlabel(self,cols): 		return self.ylabelMaker(cols)
	def makeFunc(self
					,colstr
					,cnst 	= '1'
					,h 		= true): 
		cols = colstr.split()
		return BarFunc( cols 		= cols
						,title 		= self.titleMaker(cols)
						,xLabel 	= self.xlabelMaker(cols)
						,yLabel 	= self.ylabelMaker(cols)
						,constraint = cnst
						,groupby	= [cols[0]]
						,having		= h
						,order		= self.ordMaker(cols) #order by x
						,xFunc 		= self.xfunc
						,yFunc 		= self.yfunc)

class PlotFunc(object):
	def __init__(self,cols,title,xLabel,yLabel,constraint,groupby,having,order,xyFunc,lFunc,legFunc):
		self.cols 		= cols 			# [String]
		self.title 		= title 			# String
		self.xLabel		= xLabel 		# String
		self.yLabel 	= yLabel 		# String
		self.constraint = constraint	# String
		self.groupby	= groupby 		# [String], cols that must be equal
		self.having 	= having 		# [(a,b,c,...)] -> Bool
		self.order 		= order 		# String
		self.xyFunc 	= xyFunc 		# [[a]] -> [(FLOAT,FLOAT)]
		self.lFunc 		= lFunc			# [[a]] -> [String]
		self.legFunc	= legFunc 		# [[a]] -> [String] - ideally equal for all elements in the group

	def addAxes(self,ax): ax.set_ylabel(self.yLabel);ax.set_xlabel(self.xLabel)
	def setTitle(self,ax): ax.set_title(self.title)
	def setLegend(self,ax):
		try:
			legend = ax.legend(loc='upper left', shadow=True,fontsize='small')
			legend.get_frame().set_facecolor('#00FFCC')
			legend.draggable()
		except AttributeError: pass

	def addLines(self,ax):
		for i,arg in enumerate(self.query()):
			xs = [fst(x) for x in arg]	
			ys = [snd(x) for x in arg]
			ax.plot(xs,ys,codeList[i%len(codeList)],label=arg[0][3])
			for a in arg: ax.text(a[0],a[1],a[2],size=9)

	def query(self): 
		#print 'COLS = ',self.cols
		#print 'constraint = ',self.constraint
		#print 'lFunc = ',self.lFunc
		#print 'legFunc = ',self.legFunc
		output 		= dbase.queryDistinct(self.cols,self.groupby,self.constraint,self.order) #[[(a,b,c,..)]]
		filteroutput= [x for x in output if self.having(x)]
		lines=[]
		for arg in filteroutput:
			#print 'arg = ',arg
			xys = self.xyFunc(arg)
			ls 	= self.lFunc(arg)
			legs= self.legFunc(arg)
			line=[]
			for xy,l,leg in zip(xys,ls,legs):
				if not (xy[0] is None or xy[1] is None): line.append((xy[0],xy[1],l,leg))
			if len(line)>0: lines.append(line)
		return lines

	def plot(self,ax):
		self.setTitle(ax)
		self.addAxes(ax)
		self.addLines(ax)
		self.setLegend(ax)


#############################################################################
#############################################################################

class BarFunc(object):
	def __init__(self,cols,title,xLabel,yLabel,constraint,groupby,having,order,xFunc,yFunc):
		self.cols 		= cols 			# [String]
		self.title 		= title 		# String
		self.xLabel		= xLabel 		# String
		self.yLabel 	= yLabel 		# String
		self.constraint = constraint	# String
		self.groupby 	= groupby	# String
		self.having 	= having 		# [(a,b,c,...)] -> Bool
		self.order 		= order 		# String
		self.xFunc 		= xFunc 		# [[a]] -> [String]
		self.yFunc 		= yFunc 		# [[a]] -> [FLOAT]

	def query(self): 
		output 		= dbase.queryDistinct(self.cols,self.groupby,self.constraint,self.order) #[[(a,b,c,..)]]
		def process(x): return [a for a in x if a.count(None) == 0]
		filteroutput= [process(x) for x in output if self.having(x) and process(x)!=[]]
		out = []
		for arg in filteroutput:
			print arg
			out.append((self.xFunc(arg),self.yFunc(self.cols,arg)))
		return out #[(String,Float)]

	def addAxes(self,ax): ax.set_ylabel(self.yLabel);ax.set_xlabel(self.xLabel)
	def setTitle(self,ax): ax.set_title(self.title)
	def addBars(self,ax):
		xy 		= self.query()
		y_pos 	= np.arange(len(xy)) # bar locations
		xs,ys 	= mapfst(xy),mapsnd(xy)
		ax.bar(y_pos,ys,align='center', alpha=0.5)
		ax.set_xticks(y_pos);ax.set_xticklabels(xs)
	def plot(self,ax):
		self.setTitle(ax)
		self.addAxes(ax)
		self.addBars(ax)

#############################################################################
#############################################################################


##################
# Simple Functions
##################
fst 	= lambda x: x[0];			thrd 	= lambda x: x[2]
snd 	= lambda x: x[1];			frth 	= lambda x: x[3]
mapfst 	= lambda xs: map(fst,xs);	mapthrd = lambda xs: map(thrd,xs)
mapsnd 	= lambda xs: map(snd,xs);	mapfrth = lambda xs: map(frth,xs)

fstsnd 	= lambda xs: [(x[0],x[1]) for x in xs]
sndthrd	= lambda xs: [(x[1],x[2]) for x in xs]

def appnd(cols,lname,legname): return cols + [lname,legname]
def avg(xs): return sum(xs)/float(len(xs))

def absdiff(obj): 
		#print 'obj into absdiff ',obj
		xlist,ylist = mapfst(obj),mapsnd(obj)
		if ylist.count(None)==len(ylist): return zip(xlist,ylist)
		for x,y in sorted(zip(xlist,ylist)):
			if y is not None: maxy = y
		return zip(xlist,[None if y is None else abs(y - maxy) for y in ylist])

def derivabsdiff(obj): 
		xys = absdiff(obj)

		filtxy = [xy for xy in xys if not (xy[0] is None or xy[1] is None)]
		if len(filtxy)<3: 		return [(x, None,l,ll) for x,y,l,ll in obj] 
		x = zip(*filtxy)[0];y = zip(*filtxy)[1]
		dydx = {a:b for a,b in derivxy(x,y)}

		return [(x,dydx[x] if x in dydx.keys() else None,l,ll) for x,y,l,ll in obj]

def deriv(obj): 
		print 'obj ',obj
		xys = sorted(zip(mapfst(obj),mapsnd(obj)))
		filtxy = [xy for xy in xys if not (xy[0] is None or xy[1] is None)]
		x = zip(*filtxy)[0];y = zip(*filtxy)[1]
		dydx = {a:b for a,b in derivxy(x,y)}
		return [(x,dydx[x] if x in dydx.keys() else None,l,ll) for x,y,l,ll in obj]

def derivxy(x,y):
	#Uses three-point finite difference estimate of derivative: http://www.m-hikari.com/ijma/ijma-password-2009/ijma-password17-20-2009/bhadauriaIJMA17-20-2009.pdf
	dydx = []
	if len(x) != len(set(x)): raise ValueError, "EQ constraint poorly chosen for derivative: multiple y values per x value:"+str(x) ; return 0

	for i in range(1,len(x)-1): #all values except for first and last
		h1,h2= float(x[i]-x[i-1]), float(x[i+1]-x[i])
		sumH, diffH, prodH, quotH = h1+h2, h1-h2, h1*h2, h1/h2
		f0, f1, f2 = y[i-1], y[i], y[i+1]
		dydx.append(quotH/sumH*f2 - 1/(sumH*quotH)*f0 - diffH/(prodH)*f1)
	return zip(x[1:-1],dydx)

#############################################################################

def convergence(cols,vals):
	#Uses three-point finite difference estimate of derivative: http://www.m-hikari.com/ijma/ijma-password-2009/ijma-password17-20-2009/bhadauriaIJMA17-20-2009.pdf
	#Tests whether or not yFunc derivative has a low magnitude and is decreasing over a certain range
	assert len(cols)==3;assert len(vals[0])==3
	xc,yc,dydx = d.dDict[cols[1]],d.dDict[cols[2]],[]
	print 'xc name ',xc.colname,' yc name',yc.colname, yc.derivConv
	xy = sorted([(x,y) for l,x,y in vals if not (x is None or y is None)]) #filter pairs with either x or y None
	
	if len(xy) < 5: print "Not enough data points";return 0 # not converged
	
	try:
		dydxMax = yc.derivConv[xc.colname][0] # threshold for max |dy/dx| 
		xRange  = yc.derivConv[xc.colname][1] # range (from last data point) over which |dy/dx| must be decreasing and beneath threshold
		def dydxTest(xdydxlist): 
			for x,dydx in xdydxlist:
				above = -dydxMax <= dydx
				below = dydx <= 0.0001
				if not above or not below: return False
			return True
	except KeyError: raise AttributeError, 'No derivative convergence parameters specified for deriving {0} w/r/t {1}'.format(yc.colname,xc.colname)



	xs = zip(*xy)[0];ys = zip(*xy)[1]
	xmax = xs[-2]

	xdydxs = derivxy(xs,ys)
	for i,(x,dydx) in enumerate(xdydxs):
		if dydxTest(xdydxs[i:]): return x if xmax-x >= xRange else 0
	return 0


#############################################################################
def basicName(cs): return cs[0] + ' vs ' + cs[1]
###########################################################################

pfms = [PlotFuncMaker(name= 'basic'
								,xyfunc 		= fstsnd
								,colMaker 		= appnd
								,ordMaker 		= fst
								,xlabelMaker 	= lambda cs: d.dDict[cs[0]].axislabel
								,ylabelMaker 	= lambda cs: d.dDict[cs[1]].axislabel
								,lFunc 			= mapthrd, 	legFunc = mapfrth
								,titleMaker = basicName)
		,PlotFuncMaker(name='deriv'
								,xyfunc 		= deriv
								,colMaker 		= appnd
								,ordMaker 		= fst
								,xlabelMaker 	= lambda cs: d.dDict[cs[0]].axislabel
								,ylabelMaker 	= lambda cs: 'Derivative of '+d.dDict[cs[1]].axislabel
								,lFunc 			= mapthrd, 	legFunc = mapfrth
								,titleMaker 	= lambda cs: 'Derivative of '+basicName(cs))
		,PlotFuncMaker(name='absdiff'
								,xyfunc 		= absdiff
								,colMaker 		= appnd
								,ordMaker 		= fst
								,xlabelMaker 	= lambda cs: d.dDict[cs[0]].axislabel
								,ylabelMaker 	= lambda cs: 'Absolute convergence of '+d.dDict[cs[1]].axislabel
								,lFunc 			= mapthrd, 	legFunc = mapfrth
								,titleMaker 	= lambda cs: 'Convergence of '+basicName(cs))
		,PlotFuncMaker(name='derivabsdiff'
								,xyfunc 		= derivabsdiff
								,colMaker 		= appnd
								,ordMaker 		= fst
								,xlabelMaker 	= lambda cs: d.dDict[cs[0]].axislabel
								,ylabelMaker 	= lambda cs: 'Derivative of '+d.dDict[cs[1]].axislabel
								,lFunc 			= mapthrd, 	legFunc = mapfrth
								,titleMaker 	= lambda cs: 'Derivative of absolute convergence of '+basicName(cs))

		,BarFuncMaker(name= 'basicbar'
								,xfunc 			= lambda x: x[0][0]
								,yfunc 			= lambda c,x: avg(mapsnd(x))
								,ordMaker 		= fst
								,xlabelMaker 	= lambda cs: d.dDict[cs[0]].axislabel
								,ylabelMaker 	= lambda cs: d.dDict[cs[1]].axislabel
								,titleMaker 	= lambda cs: cs[0] + ' vs Average ' + cs[1])
		,BarFuncMaker(name= 'convergence'
								,xfunc 			= lambda x: x[0][0]
								,yfunc 			= convergence
								,ordMaker 		= fst
								,xlabelMaker 	= lambda cs: d.dDict[cs[1]].axislabel
								,ylabelMaker 	= lambda cs: d.dDict[cs[2]].axislabel
								,titleMaker 	= lambda cs: 'Convergence of '+ cs[1] + ' vs ' + cs[2])
		]

makerDict= {p.name:p for p in pfms}


