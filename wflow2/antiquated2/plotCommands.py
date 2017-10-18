from constraint import *
from printParse import *
import details as d
########
def xyX(lst): return not (lst[0] is None or lst[1] is None)
def len2(xs): return 1 <  len(set(([x[0] for x in xs if xyX(x)])))


################################################################
# Kinds of data summaries (on one list) :: [[Double]] -> Double
################################################################

def meanAbsolute(cols,xs):  return sum(map(abs,xs[0]))/float(len(xs[0]))
def mean(cols,xs):          return sum(xs[0])/float(len(xs[0]))
def gMeanAbsolute(cols,xs): return mstatt.gmean(map(abs,xs[0]))
def variance(cols,xs):      raise NotImplementedError


# Bar Plots

avgTime = {'maker':'basicbar'
		,'colstr':'xc timeperstep'
		,'cnst': 	AND([BEEF,PW(500),KPTDEN(2)])}


ZnPdconvEng = {'maker': 'convergence'
			,'colstr':'name pw err_A'
			,'cnst': AND([MBEEF,KPTDEN(2),LATTICEOPT,OR([CU,BE])])}

ZnPdEngAbsDiff = {'maker': 		'absdiff'
			,'colstr':		'pw err_A'
			,'lname': 		'fwid'
			,'legname':		'name'
			,'gb':			'name'
			,'cnst':		AND([MBEEF,KPTDEN(2),OR([CU,BE]),LATTICEOPT])}

ZnPdConvEng_d = {'maker': 		'derivabsdiff'
			,'colstr':		'pw err_A'
			,'lname': 		'fwid'
			,'legname':		'name'
			,'gb':			'name'
			,'cnst':		AND([MBEEF,KPTDEN(2),OR([CU,BE]),LATTICEOPT])}

#maeLatParamA = BarFunc(r"MAE Lattice Constant 'a', $\AA$",['err_A'],gMeanAbsolute)
#maeBulkMod   = BarFunc(r"MAE Bulk Modulus, GPa",['errBM'],gMeanAbsolute)

#errAvsStructure = {'summaryFunc': 	BarFunc('MAE Lattice Error',['err_A'],meanAbsolute)
#					,'labelFunc': 	'structure'
#					,'eq': 			structs(AND([LATTICEOPT,GPAW,PBE]))
#					}

################
## EQ Commands
###############

exampleDict = {'maker':'absdiff'
				,'colstr':'pw err_A fwid name'
				,'cnst':'xc=\'mBEEF\' and kptden=2'
				,'gb':'name'}

timeperstepVsPW = {'maker':'basic'
				,'colstr':'pw timeperstep blank blank'
				,'cnst':"xc='mBEEF' and kptden=2 and jobkind='latticeopt'"
				,'gb':'name'}

h2VsPw = 	{'maker': 	'absdiff'
			,'colstr':	'pw raw_energy blank blank'
			,'gb':		'dftcode'
			,'cnst':	H2}


liErrBM = {'maker': 		'absdiff'
			,'colstr':		'pw err_BM blank calclabel'
			,'gb':			'xc dftcode kptden'
			,'cnst':		AND([LI,BULKMOD])}


mbeefErrA = {'maker': 		'absdiff'
			,'colstr':		'pw err_A fwid name'
			,'gb':			'name'
			,'cnst':		AND([KPTDEN(2),BULKMOD,GPAW,MBEEF])}

mbeefErrAderiv = {'maker': 		'derivabsdiff'
			,'colstr':		'pw err_A'
			,'lname': 'fwid'
			,'legname': 'name'
			,'gb':			'name'
			,'cnst':		AND([KPTDEN(2),BULKMOD,GPAW,MBEEF])}

convErrBM = {'maker': 'convergence'
			,'colstr':'name pw err_BM'
			,'cnst': AND([KPTDEN(2),BULKMOD,GPAW,PBE])}

mBEEFbfit = {'maker': 		'basic'
			,'colstr':		'pw bfit fwid name'
			,'gb':			'name'
			,'cnst':		AND([KPTDEN(2),BULKMOD,GPAW,MBEEF])}

hydride = {'maker':'basic'
			,'colstr':	'hfrac eform fwid metalstoich_str'
			,'gb':		'metalcomp_str kptden'
			,'h': len2
			,'cnst':		AND([BEEF,PW(500),BULK,QE,RELAXORLAT])}


molecule = {'maker':'absdiff'
			,'colstr':	'pw raw_energy fwid name'
			,'gb':		'name xc dftcode'
			#,'h': len2
			,'cnst':		AND([COMPLETED,RELAX,MOLECULE])}
