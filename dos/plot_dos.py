#!/usr/bin/env python276


import numpy as np
import matplotlib.pyplot as plt
from pylab import *

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, MaxNLocator
import matplotlib
from scipy.fftpack import hilbert, ihilbert
from ase.io import read, write
from os import listdir
from os.path import isfile, join

import os, sys, pickle

path   = os.getcwd() + '/'
files  = [f for f in listdir(path) if isfile(join(path, f))] ; files.sort()


for f in files:
    if f[-5:]=='.traj': 
        trajname=f




Len = len(sys.argv)



fsize = 10

plotUB =False
plotmax=False

font = {'family' : 'serif', #options: 'serif' (e.g. Times), 'sans-serif' (e.g. Helvetica), 'cursive' (e.g. Zapf-Chancery), 'fantasy' (e.g. Western), and 'monospace' (e.g. Courier)
        'style'  : 'normal', #options: normal (or roman), italic or oblique
        'weight' : 'normal', #options: 13 values: normal, bold, bolder, lighter, 100, 200, 300, ...900.
        'size'   : fsize
#   'serif'   : ['Times New Roman'],
#   'sans-serif' : ['Helvetica Neue'],
    }

#fig = plt.figure(figsize=(3.5,7.5))
#fig.set_size_inches(3.25,5)

matplotlib.rc('font', **font)
matplotlib.rc("text", usetex=True)
#matplotlib.rc("text.latex", preamble=["\\usepackage{helvet}\\usepackage{sansmath}\\sansmath"])

axes_font = {'labelsize' : fsize, #: medium  # fontsize of the x any y labels
             'titlesize' : fsize #fontsize of the axes title
           }

mathtext_prop = {
         #'default' : 'regular',
         'fontset' : 'stix',
        }

matplotlib.rc('axes', **axes_font)
matplotlib.rc('mathtext', **mathtext_prop)

plt.rc('legend',**{'fontsize':fsize})
plt.rc('font',**font)
plt.rc('axes',**axes_font)

npts=2000
xlabel = '$\mathit{\epsilon}-\mathit{\epsilon}_\mathrm{F}$ (eV)'
ylabel = '$\mathsf{d}$-projected DOS on Mo'
xlim = [-10,10]#[-5,5]
ylim = [0,5]#[0,5]
xlim2 = [-1.5,1.5]
ylim2 = [0,0.15]
MaxN = 2    #number of ticks
fontsize = fsize    #font size

sym = 0
scale = 1
hscale = 2.3

color = '#4cae4c'
color2 = '#3276b1'
shade = '#bfbfbf'
ubcolor = '#bfbfbf'
pastelred='#e14f4f'
maxp='#275b89'
hcolor = pastelred

zmax=1.0 # most likely need to change this.
zO=zmax
zH=zmax

Len = len(sys.argv)

exclude=[]  # H and O atoms to exclude, e.g. H and O that are part of adsorbates and not the water layer. enter in system arguments last.
idxH = []   # idx of H that needs to be included
save=True

outname=''
dosname = 'dos.pickle'
if Len>1:
    for i in range(1,Len):
        if sys.argv[i-1]=='-z':
            zmax=float(sys.argv[i])
            zO=zmax
            zH=zmax
        if sys.argv[i-1]=='-zO':
            zO=float(sys.argv[i])
        if sys.argv[i-1]=='-zH':
                        zH=float(sys.argv[i])
        if sys.argv[i-1]=='-e':
            for j in range(i,Len):
                exclude.append(int(sys.argv[j]))
        if sys.argv[i-1]=='-o':
            outname=sys.argv[i]
        if sys.argv[i-1]=='-t':
            trajname=sys.argv[i]
        if sys.argv[i-1]=='-d':
            dosname=sys.argv[i]
        if sys.argv[i]=='-ns':
            save=False
        if sys.argv[i-1]=='-u':
            dosname = 'dos_'+sys.argv[i]+'.pickle'
            outname = sys.argv[i]+'_'
        if sys.argv[i-1]=='-iH':
            for j in range(i,Len):
                idxH.append(int(sys.argv[j]))

if save:
    outfilename = 'pdos_'+outname+'.txt'
    outfile=open(outfilename, 'w+')
    print >>outfile, 'trajectory file:'+trajname
    print >>outfile, 'zO:'+str(zO)
    print >>outfile,'zH:'+str(zH)
else:
    print  'trajectory file:'+trajname
    print  'zO:'+str(zO)
    print 'zH:'+str(zH)



def plot_pdos( path, title,fignum, figsave):


    n=len(title)

    fig = plt.figure(fignum)

    for imagenum in range(0,n):
        f = open(path[imagenum]+dosname)
        energies, dos, pdos = pickle.load(f)
        f.close()
        slab = read(path[imagenum]+trajname)

        atoms_O=[]
        atoms_H=[]
        if not idxH:
            for a in slab :
                if a.symbol== 'O' and a.z > zO:
                    atoms_O.append(a.index)
                elif a.symbol== 'H' and a.z > zH:
                    atoms_H.append(a.index)

            for k in range(0,len(exclude)):   #delete from O and H lists what we want excluded
                list_atoms_O = range(0,len(atoms_O))
                list_atoms_O.reverse()
                list_atoms_H = range(0,len(atoms_H))
                list_atoms_H.reverse()
                #print 'atoms_O', atoms_O
                #print 'atoms_H', atoms_H
                for i in list_atoms_O:
                    if atoms_O[i]==exclude[k]:
                        atoms_O.pop(i)

                for i in list_atoms_H:
                    if atoms_H[i]==exclude[k]:
                        atoms_H.pop(i)
        else:
            atoms_H = idxH

        if save:
            print >>outfile,'excluded atoms: '+str(exclude)
            print >>outfile,'atoms_O: '+str(atoms_O)
            print >>outfile,'atoms_H: '+str(atoms_H)

        p1 = plt.subplot(n,1,imagenum+1)

        plotdos=[]
        for i in range(0,len(atoms_H)):
            natom = atoms_H[i]
            dos_H = pdos[natom]['s'][sym]
            if i==0:
                plotdos = dos_H
                #print plotdos
            else:
                plotdos = np.add(plotdos , dos_H)
                # if i==1:
#                   print dos_H
#                   print plotdos

        plt.plot(energies,plotdos,label='$s$-dos on H ',color=pastelred)

        f = open('H_dos.txt', 'w')
        f.write(','.join([str(e) for e in energies]))
        f.write('\n')
        f.write(','.join([str(d) for d in plotdos]))
        f.close()
        del plotdos

        if atoms_O:
            for i in range(0,len(atoms_O)):
                natom = atoms_O[i]
                dos_O = pdos[natom]['s'][sym]
                if i == 0:
                    plotdos = dos_O
                else:
                    plotdos = np.add(plotdos , dos_O)


            plt.plot(energies,plotdos,label='$s$-dos on O ',color=color2)
            del plotdos


            for i in range(0,len(atoms_O)):
                natom = atoms_O[i]
                dos_O = pdos[natom]['p'][sym]
                if i == 0:
                    plotdos = dos_O
                else:
                    plotdos = np.add(plotdos , dos_O)


            plt.plot(energies,plotdos,label='$p$-dos on O ',color=color)
            del plotdos

        p1.yaxis.set_major_locator(MaxNLocator(MaxN))
        plt.ylabel(title[imagenum],fontdict={'fontsize':fontsize})
        plt.setp(p1.get_yticklabels(), visible=False)
        #plt.setp(p1.get_xticklabels(), visible=False)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.plot([0, 0], [ylim[0],ylim[1]], 'k:', lw=0.7)
        if imagenum==0:
            leg=plt.legend(loc=1,prop={'size':10},numpoints=1)

    plt.xlabel(xlabel,fontdict={'fontsize':fontsize})
    if save:
        fig.savefig(figsave)
    return atoms_O, atoms_H


atoms_O, atoms_H = plot_pdos(['./'],['arb. units'],1,'pdos_'+outname+'.pdf')

if save:
    outfile.close()
    outputfile=open(outfilename, "r")
    printout=outputfile.readlines()
    for line in printout:
            print line

plt.show()