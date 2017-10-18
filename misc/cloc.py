import sys,os
from prettytable import PrettyTable

"""
Reads through scripts folder, makes prettytable files for each subdirectory

Requires that cloc is installed.
"""

def parseCloc(output):
	output_2 = output.readlines()[-2].split()
	try: return [output_2[i] for i in [0,-2,-1]]
	except IndexError: pass

def main():

	scriptsPath = '/Users/kbrown/scripts/' if len(sys.argv)<2 else os.path.abspath(sys.argv[1])+'/'

	for d in os.listdir(scriptsPath):
		if os.path.isdir(scriptsPath+d):
			with open(scriptsPath+d+'/CodeLengthSummary.txt','w') as f:
				
				table = PrettyTable(['File','Language','Comments','Code','Definitions','Size,kB'])
				totCom,totCode,totDef,totSize=0,0,0,0

				for x in os.listdir(scriptsPath+d):
					pth = scriptsPath+d+'/'+x
					if os.path.isfile(pth):
					
						parsed = parseCloc(os.popen('cloc '+pth))  
						if parsed is None or parsed[0] not in ['Python','Bourne']: pass
						else: 
							statinfo = os.stat(pth)
							size = int(round(statinfo.st_size/1000.0))

							with open(pth,'r') as ff: 
								defs = ff.read().count('def ')

							table.add_row([x]+parsed+[defs,size])
							totCom += int(parsed[1]) ;	totCode+= int(parsed[2])
							totDef += defs ; 			totSize += size
	
				table.add_row(['Total','-',totCom,totCode,totDef,totSize])
				f.write(str(table))						

if __name__ == '__main__':
	main()