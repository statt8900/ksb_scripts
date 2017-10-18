from constraint import *
from printParse import *
import details as d
########
class BarFunc(object):
	def __init__(self,name,details,kind):
		self.name 	 = name       # String
		self.details = details    # [Detail] 
		self.kind 	 = kind       # ([String], [[Double]]) -> Double

################################################################
# Kinds of data summaries (on one list) :: [[Double]] -> Double
################################################################

def meanAbsolute(cols,xs):  return sum(map(abs,xs[0]))/float(len(xs[0]))
def mean(cols,xs):          return sum(xs[0])/float(len(xs[0]))
def gMeanAbsolute(cols,xs): return mstatt.gmean(map(abs,xs[0]))
def variance(cols,xs):      raise NotImplementedError
def convergence(cols,vals):
	#Uses three-point finite difference estimate of derivative: http://www.m-hikari.com/ijma/ijma-password-2009/ijma-password17-20-2009/bhadauriaIJMA17-20-2009.pdf
	#Tests whether or not yFunc derivative has a low magnitude and is decreasing over a certain range

	assert len(cols)==2;assert len(vals[0])==2
	print cols
	xc,yc,dydx = d.dDict[cols[0]],d.dDict[cols[1]],[]
	try:
		ytol	= yc.ytol
		dydxMax = yc.derivConv[xc.colname][0] # threshold for max |dy/dx| 
		xRange  = yc.derivConv[xc.colname][1] # range (from last data point) over which |dy/dx| must be decreasing and beneath threshold
	except KeyError: raise AttributeError, 'No derivative convergence parameters specified for deriving {0} w/r/t {1}'.format(yc.colname,xc.colname)

	xy = [(x,y) for x,y in vals if not (x is None or y is None)]
	
	if len(xy) < 4: print "Not enough data points";return 0 # not converged
	
	maxX,maxY = xy[-1][0], xy[-1][1]

	x = zip(*xy)[0];y = zip(*xy)[1]
	if len(x) != len(set(x)): print "EQ constraint poorly chosen: multiple y values per x value:"+str(x) ; return 0
	print 'x ',x
	for i in range(1,len(x)-1): #all values except for first and last
		print 'i ',i
		h1,h2= float(x[i]-x[i-1]), float(x[i+1]-x[i])
		sumH, diffH, prodH, quotH = h1+h2, h1-h2, h1*h2, h1/h2
		f0, f1, f2 = y[i-1], y[i], y[i+1]
		dydx.append(-1/(sumH*quotH)*f0 - diffH/(prodH)*f1 + quotH/sumH*f2)
	
	xDYDX = zip(x[1:-1],dydx)
	
	if any([x_i > maxX and DYDX > dydxMax for x_i,DYDX in xDYDX]): print "not converged due to dydxMax and xRange criteria";return 0
	
	for x,y in xy:
		if y-maxY < ytol: return x
	print 'Derivative meets convergence criteria, but absolute difference from last data point too high'
	return 0

# Bar Plots

errAvsStructure = {'summaryFunc': 	BarFunc('MAE Lattice Error',['err_A'],meanAbsolute)
					,'labelFunc': 	'structure'
					,'eq': 			structs(AND([LATTICEOPT,GPAW,PBE]))
					}

convErrA = {'summaryFunc': BarFunc('Convergence',['pw','err_A'],convergence)
			,'labelFunc': 'name'
			,'eq': names(AND([GPAW,PBE,LATTICEOPT,KPTDEN(2)]))}
################
## EQ Commands
###############
timeperstepVsPW = {'plotFunc': 		'pw_vs_timestep'
					,'legName':		'name'
					,'labelName':	'calclabel'
					,'eq':			names(AND([LATTICEOPT,ZEROFORCE]))
					,'const':		AND([LATTICEOPT,ZEROFORCE])
					,'title': 		'Electronic Time per step vs pw (zero force structures)'}

h2VsPw = 	{'plotFunc': 	'pw_vs_raw_energy'
			,'legName':		'calclabel'
			,'labelName':	'fwid'
			,'eq':			equalcalcs('1')
			,'const':		H2
			,'title': 		'Energy convergence for H2'}

liErrA = 	{'plotFunc':		'pw_vs_err_A'

			,'legName':		'calclabel'
			,'labelName':	'fwid'
			,'eq':			xcs('1')
			,'const':		LI
			,'title': 		'Energy convergence for Li error in lattice constant'}

alErrA = 	{'plotFunc':	'pw_vs_err_A'
			,'yName':		'err_A'
			,'legName':		'calclabel'
			,'labelName':	'fwid'
			,'eq':			xcs('1')
			,'const':		AL
			,'title': 		'Energy convergence for Al error in lattice constant'}

PBEerrA = {'plotFunc':		'pw_vs_err_A'
			#,'legName':	'name'
			,'labelName':	'name'
			,'eq':			names('1')
			,'const':		PBE
			,'title': 		'PBE pw convergence for error in lattice constant'}

mBEEFerrA = {'plotFunc':	'pw_vs_err_A'
			#,'legName':	'name'
			,'labelName':	'name'
			,'eq':			names('1')
			,'const':		MBEEF
			,'title': 		'mBEEF pw convergence for error in lattice constant'}

PBEerrAconv = {'plotFunc':	'absdiff pw_vs_err_A'
			#,'legName':	'name'
			,'labelName':	'name'
			,'eq':			names('1')
			,'const':		AND([PBE,OR([LI,AL,BE,CU])])
			,'title': 		'PBE pw convergence for error in lattice constant'}

mBEEFerrAconv = {'plotFunc':'absdiff pw_vs_err_A'
			#,'legName':		'name'
			,'labelName':	'name'
			,'eq':			[names('1')]
			,'const':		AND([MBEEF,KPTDEN(2),LATTICEOPT])#,OR([LI,AL,BE,CU])])
			,'title': 		'mBEEF pw convergence for error in lattice constant'}

mBEEFerrBM = {'plotFunc':	'absdiff pw_vs_err_BM'
			,'legName':		None#'name'
			,'labelName':	'name'
			,'eq':			[names('1')]
			,'const':		AND([BULKMOD,KPTDEN(4),MBEEF])
			,'title': 		'mBEEF pw convergence for error in bulk modulus'}

mBEEFbfit = {'plotFunc':	'pw_vs_bfit'
			#,'legName':	'name'
			,'labelName':	'name'
			,'eq':			[names('1')]
			,'const':		AND([BULKMOD,MBEEF])
			,'title': 		'mBEEF pw convergence for bulk modulus quadratic fit'}

interstitialEnergy = {'plotFunc': 'hfrac_vs_eform'
			,'yName':		'raw_energy'
			,'legName':		'name'
			,'labelName':	'symbols_str'
			,'eq':			[metalspecies('1')]
			,'const':		AND([PW(500),BULK,KPTDEN(2),QE,RELAXORLAT])
			,'title': 		'interstitial energy'
			,'singletons': False}
