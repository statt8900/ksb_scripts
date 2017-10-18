#!/usr/bin/env python


import pickle
import sys

pickle_file = sys.argv[1] 
display = 'electronic energy'
Len = len(sys.argv)
if Len > 2:
	for i in range(2,Len): 
		if sys.argv[i] == "script":
			display = "calculation script"
		if sys.argv[i] == "path":
			display = "path"
		if sys.argv[i] == "vibrations":
			display = "vibrations"				

print pickle.load(open(pickle_file))[display]

