#!/usr/bin/env python

from ase.io import read,write
import sys,os,random

default = 'qn.traj'
alist = {}
if sys.argv[1:] == []:
    try:
        alist[default] = read(default)
    except:
        print 'Usage: python screenshot.py traj1 traj2 ...'
else:
    for arg in sys.argv[1:]:
        try:
            alist[arg] = read(arg)
        except:
            print 'Invalid traj file: ' + arg


for i,file in enumerate(alist):
    textures = ['simple' for i in range(len(alist[file]))]
    name=file[0:-5]
    write(name+'_ortho.pov',alist[file],rotation='-60x,-45y,-20z',run_povray=True,textures=textures,canvas_height=1000)
    os.remove(name+'_ortho.pov')
    os.remove(name+'_ortho.ini')

    write(name+'_top.pov',alist[file],rotation='0x,0y,0z',run_povray=True,textures=textures,canvas_height=1000)
    os.remove(name+'_top.pov')
    os.remove(name+'_top.ini')

    write(name+'_side.pov',alist[file],rotation='-90y',run_povray=True,textures=textures,canvas_height=1000)
    os.remove(name+'_side.pov')
    os.remove(name+'_side.ini')

    my_randoms = [str(x) for x in random.sample(xrange(90), 9)]
    
    randrot = []
    
    for i in range(3):
        randrot.append(','.join([a+b for a,b in zip(my_randoms[3*i:3*i+3],"xyz")]))

    for j in range(3):
        write(name+'_random%d.pov'%j,alist[file],rotation=randrot[j],run_povray=True,textures=textures,canvas_height=1000)
        os.remove(name+'_random%d.pov'%j)
        os.remove(name+'_random%d.ini'%j)