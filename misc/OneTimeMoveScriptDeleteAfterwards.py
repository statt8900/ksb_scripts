import ase.io as aseio
#OUTCAR 2.588
#INIT  2.624
scale = 2.588/2.624

atoms = aseio.read('initial.traj')

moveDict={57:[59],61:[58]
		,50:[48,49]
		,56:[55,54]
		,53:[52,51]
		,23:[21,24,22]
		,26:[25,60]
		,47:[45,46]
		,20:[18,19]}

invertedDict = {value: key for key in moveDict for value in moveDict[key]}

vectorDict = {k:atoms[k].position-atoms[v].position for k,v in invertedDict.items()}

unScaledAtoms = invertedDict.keys()
print unScaledAtoms
scaledAtoms = [i for i in range(1,62) if i not in unScaledAtoms]

def scaleXY(coord): 	
	coord[0]*=scale ; coord[1]*=scale
	return coord

atoms.set_cell(scaleXY(atoms.get_cell()))

for a in atoms:
	if a.index in scaledAtoms: 
			print 'unscaled atom ',a.index,a.symbol, ' at ',a.position[:2]
			a.position=scaleXY(a.position)
			print 'new position ',a.position[:2]

for a in atoms:
	if a.index in unScaledAtoms: 
		a.position=atoms[invertedDict[a.index]].position + vectorDict[a.index]

aseio.write('scaled.traj',atoms)

