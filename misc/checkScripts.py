import os
from importlib import import_module
root = '/Users/kbrown/scripts/'
for folder in ['fireworks','db']:
	if os.path.isdir(root+folder):
		for d in os.listdir(root+folder):
			try:
				print 'splitting ',d
				fs=d.split('.')
				f,ext = fs[0],fs[1]
				if ext == 'py': py = True
				else: py = False
			except IndexError:
				py = False
			if py:
				print 'checking ',f
				import_module(f)
