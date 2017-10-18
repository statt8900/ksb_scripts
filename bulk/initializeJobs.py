import sqlite3
from db_insert import query1
# Keywords for Atoms.info that may be used:
# name    (assumed identical to path if not specified)
# bravais (assumed triclinic if not specified)
#for d in os.listdir('/scratch/users/ksb/db/traj'): 	insertObject(BulkTraj(d))


magList = [('Fe2','bcc'),('Ni','fcc'),('Co2','hcp'),('MnO','rocksalt')
			,('MnS','rocksalt'),('MnN','rocksalt'),('MnC','rocksalt')
			,('FeC','rocksalt'),('FeAl','cscl'),('CrC','rocksalt'),('CrN','rocksalt')]

magElems = ['Fe','Mn','Cr','Co','Ni']

############################
#Aux functions to help rules
############################
def stoich(j): return parseStoich(query1('bulk','stoich','pth',j[0]))
######################################################
#Rules: (trajPath,CalcInd,CalcInd,dftName,Int) -> Bool
######################################################

def noJobsOver10h(j): return j[4] <= 10

def magGuard(j):  
	stoich = stoich, in magElems



rules=[noJobsOver10h]
###################################
timeLimDomain = [1,10,20]

combs = product(trajDomain,nC,nC,['gpaw','quantumespresso'],timelimDomain)



for c in combs:
	comb = [time.time()]+['kris']+[time.time()]+['bulkrelax']+['']+c
	insertObject(Job(*comb))

