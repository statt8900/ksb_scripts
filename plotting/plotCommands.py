from sqlShortcuts import *

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
		,'cnst': 	AND(BEEF,PW(500),KPTDEN(2))}


##convEng = {'maker': 'convergence'
#			,'colstr':'job_name planewave_cutoff error_lattice_A'
#			,'cnst': AND([BEEF,KPTDEN(4),LATTICEOPT,OR(AL,IR)])}

#########

test = {'maker':'absdiff'
				,'colstr':'planewave_cutoff error_BM fwid job_name'
				,'cnst':AND(KPTDEN(2),BULKMOD,MBEEF,PAWPSP)
				,'gb':'job_name psp_ksb'}





test2 = {'maker':'absdiff'
				,'colstr':'planewave_cutoff error_lattice_A job_name psp_ksb'
				,'cnst':AND(KPTDEN(2),PAWELEM,SG15,LATTICEOPT,MBEEF)
				,'gb':'job_name psp_ksb'}




engAbsDiff = {'maker': 		'absdiff'
			,'colstr':		'planewave_cutoff error_lattice_A'
			,'lname': 		'fwid'
			,'legname':		'job_name'
			,'gb':			'job_name'
			,'cnst':		AND(BEEF,KPTDEN(4),OR(AL,IR),LATTICEOPT)}



convEng_d = {'maker': 		'derivabsdiff'
			,'colstr':		'planewave_cutoff error_lattice_A'
			,'lname': 		'fwid'
			,'legname':		'job_name'
			,'gb':			'job_name'
			,'cnst':		AND(BEEF,KPTDEN(4),OR(AL,IR),LATTICEOPT)}






#maeLatParamA = BarFunc(r"MAE Lattice Constant 'a', $\AA$",['error_lattice_A'],gMeanAbsolute)
#maeBulkMod   = BarFunc(r"MAE Bulk Modulus, GPa",['errBM'],gMeanAbsolute)

#errAvsStructure = {'summaryFunc': 	BarFunc('MAE Lattice Error',['error_lattice_A'],meanAbsolute)
#					,'labelFunc': 	'structure'
#					,'eq': 			structs(AND([LATTICEOPT,GPAW,PBE]))
#					}

################
## EQ Commands
###############

exampleDict = {'maker':'absdiff'
				,'colstr':'planewave_cutoff error_lattice_A fwid job_name'
				,'cnst':'xc=\'mBEEF\' and kptden_ksb=2'
				,'gb':'job_name'}

timeperstepVsPW = {'maker':'basic'
				,'colstr':'planewave_cutoff timeperstep blank blank'
				,'cnst':"xc='mBEEF' and kptden_ksb=2 and jobkind='latticeopt'"
				,'gb':'job_name'}

h2VsPw = 	{'maker': 	'absdiff'
			,'colstr':	'planewave_cutoff raw_energy blank blank'
			,'gb':		'dftcode'
			,'cnst':	H2}


liErrBM = {'maker': 		'absdiff'
			,'colstr':		'planewave_cutoff err_BM blank calclabel'
			,'gb':			'xc dftcode kptden_ksb'
			,'cnst':		AND(LI,BULKMOD)}


mbeefErrA = {'maker': 		'absdiff'
			,'colstr':		'planewave_cutoff error_lattice_A fwid job_name'
			,'gb':			'job_name'
			,'cnst':		AND(KPTDEN(2),BULKMOD,GPAW,MBEEF)}

mbeefErrAderiv = {'maker': 		'derivabsdiff'
			,'colstr':		'planewave_cutoff error_lattice_A'
			,'lname': 'fwid'
			,'legname': 'job_name'
			,'gb':			'job_name'
			,'cnst':		AND(KPTDEN(2),BULKMOD,GPAW,MBEEF)}

convErrBM = {'maker': 'convergence'
			,'colstr':'job_name planewave_cutoff err_BM'
			,'cnst': AND(KPTDEN(2),BULKMOD,GPAW,PBE)}

mBEEFbfit = {'maker': 		'basic'
			,'colstr':		'planewave_cutoff bulk_modulus_quadfit fwid job_name'
			,'gb':			'job_name'
			,'cnst':		AND(KPTDEN(2),BULKMOD,GPAW,MBEEF)}

hydride = {'maker':'basic'
			,'colstr':	'hfrac eform fwid metalstoich_str'
			,'gb':		'metalcomp_str kptden_ksb'
			,'h': len2
			,'cnst':		AND(BEEF,PW(500),BULK,QE,RELAXORLAT)}


molecule = {'maker':'absdiff'
			,'colstr':	'planewave_cutoff raw_energy fwid job_name'
			,'gb':		'job_name xc dftcode'
			,'cnst':		AND(RELAX,MOLECULE)}
