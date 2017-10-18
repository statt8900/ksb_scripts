import os,datetime
from  objectClass import Object

"""
Classes for the traj data that is used to create a Bulk or Surface structure.
They are distinguished from 'pure' Bulk and Surface objects that are completely 
dictated by their internal parameters.
"""

class TrajDataBulk(Object):
	def __init__(self,name,path,bravais,createdBy='Me',created=None,tags=[]): 
		if created is None: created=getTimeStr(path)
			
		self.name      = name        # String
		self.bravais   = bravais     # Need to know what cell dimensions are free to optimize
		self.path      = path        # FilePath (should reference ~, not absolute, in order to work on clusters)
		self.created   = created     # String, YYMMDDHHMMSS
		self.createdBy = createdBy   # String
		self.tags      = tags        # shortcuts for filtering


class TrajDataSurf(Object):
	def __init__(self,name,path,scale=(1,1,1),vac=10,created=None,createdBy='Me',tags=[]): 
		if created is None: created=getTimeStr(path)
		self.name      = name        # String
		self.path      = path        # FilePath (should reference ~, not absolute, in order to work on clusters)
		self.created   = created     # String, YYMMDDHHMMSS
		self.createdBy = createdBy   # String
		self.scale     = scale     
		self.vac       = vac
		self.tags      = tags        # shortcuts for filtering


def getTimeStr(path):
	utc=  os.stat(path).st_mtime
	t = datetime.datetime.fromtimestamp(utc)
	return "{0:0>2}{1:0>2}{2:0>2}{3:0>2}{4:0>2}{5:0>2}".format((t.year)%100
								,t.month
								,t.day
								,t.hour
								,t.minute
								,t.second)