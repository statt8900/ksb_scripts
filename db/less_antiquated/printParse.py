import re

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

def printAds(x): raise NotImplementedError
def parseAds(x): raise NotImplementedError