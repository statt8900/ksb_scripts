import re,sys
from ast import literal_eval
import numpy as np

"""
Auxillary functions that are occasionally needed to translate between strings and things
"""

def parseStoich(x): 
	if not x: return {}
	return {alpha(y) : int(digit(y)) for y in cleanSplit(x,'-')}
def digit(x):    return re.sub("[^0-9]","",x)
def alpha(x):    return re.sub("[^a-zA-Z]", "", x)

def printTime(floatHours): 
	intHours = int(floatHours)
	return "%02d:%02d" % (intHours,(floatHours-intHours)*60)


def alphaNumSplit(x): return (alpha(x),digit(x))

def parseChemicalFormula(x): return {ele:(1 if n is '' else int(n)) for ele,n in re.findall(r'([A-Z][a-z]*)(\d*)', x)}

def printAds(x): 
	return '__'.join([k+'_'+'-'.join(v) for k,v in x.items()])
def parseAds(x):
	items	= x.split('__')
	keys 	= [i.split('_')[0] for i in items]
	vals	= [i.split('_')[1].split('-') for i in items]
	return {k:v for k,v in zip(keys,vals)}

def parseAds(x): raise NotImplementedError

def parseNp(x):
	if x is None: return None
	stripped = ''.join(str(x).split('\n'))
	strIn = stripped.replace('] [',' ; ').replace('[','').replace(']','').replace('(','').replace(')','')
	mOut= np.matrix(strIn,dtype=np.float_)	
	return mOut

