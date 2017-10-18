#!/usr/bin/python 
import sys,os

if len(sys.argv) < 2: user = 'ewh'
else: user = sys.argv[1]

os.system('squeue | grep ' + user + ' | grep iric | grep R | grep -v Resources > tmp.txt')
f = open('tmp.txt')
lines = f.readlines()
f.close()
os.system('rm tmp.txt')

nodes = 0
nodelist = []
for line in lines:
	node = line.split(' ')[-1].rstrip()
	num = int(line.split(' ')[-2])
	if num > 1:
		nodes += num
	else:
		if node in nodelist:
			continue
		else:
			nodes += 1
			nodelist.append(node)

print 'User '+ user + ' is using %i nodes on iric'%nodes