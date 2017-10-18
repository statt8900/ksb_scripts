#######################################################
## Based on the tutorial by Hongliang Xin
## Topic: Plotting the density of states data from ase-espresso
## dependencies: pickle, matplotlib, numpy
#######################################################
#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import os, sys, pickle
##########################################
##########################################
#INPUTS
##########################################
##########################################
num = [2,3]                 #Atom numbers to sum over
orb = ['s','p','d']       #Orbitals to include projection
xl,xh = -10,30 #Axis limits
yl,yh =   0,5  #Axis limits
filename = 'dos_PdAu'
figLabel = 'PdAu full DoS'
##########################################
##########################################
axLim = [xl,xh,yl,yh]      #Axis parameters
xtick = np.linspace(xl,xh,10,endpoint=False).tolist()
ytick = np.linspace(yl,yh,10,endpoint=False).tolist()

def read_dos(dir,tetrahedra=False):
    import pickle
    try:
        if tetrahedra==True:
            f = open(dir + '/dos_tetra.pickle')
            energies, dos, pdos = pickle.load(f)
            f.close()
        else:
            f = open(dir + '/dos.pickle')
            energies, dos, pdos = pickle.load(f)
            f.close()
    except:
        print "No Density of States DATA Found."
        sys.exit(1)
    return energies, dos, pdos

rcParams['figure.figsize'] = 6*1.67323,4*1.67323
rcParams['ps.useafm'] = True
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

rcParams['pdf.fonttype'] = 42
matplotlib.rc('xtick.major', size=6)
matplotlib.rc('xtick.minor', size=3)
matplotlib.rc('ytick.major', size=6)
matplotlib.rc('ytick.minor', size=3)
matplotlib.rc('lines', markeredgewidth=0.5*2)
matplotlib.rc('font', size=12*2)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('PDOS (Arb. Unit)')

energies,dos,pdos = read_dos('./')

resov_dos = np.zeros((pdos[num[0]][orb[0]][0]).size) #Initialize vector

for i in num:
    for o in orb:
        resov_dos = np.sum([resov_dos,pdos[i][o][0]],axis=0)

ax.plot(energies,resov_dos,
        color='k',
        linestyle='-',
        label=figLabel)

ax.axis(axLim)
plt.xticks(xtick); plt.yticks(ytick)
ax.minorticks_on()

leg=plt.legend(loc=1,prop={'size':24},numpoints=1)
leg.draw_frame(False)

left   = 0.125  # the left side of the subplots of the figure
right  = 0.95    # the right side of the subplots of the figure
bottom = 0.2   # the bottom of the subplots of the figure
top    = 0.95      # the top of the subplots of the figure
wspace = 0.35   # the amount of width reserved for blank space between subplots
hspace = 0.35   # the amount of height reserved for white space between subplots
subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

plt.savefig(filename+'.png', format='png')
plt.show()