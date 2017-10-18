#!/home/vossj/suncat/bin/python
#above line selects special python interpreter which knows all the paths
#SBATCH -p iric,owners,normal
#SBATCH -x gpu-14-1,sh-20-35
#################
#set a job name
#SBATCH --job-name=myjob
#################
#a file for job output, you can check job progress
#SBATCH --output=myjob.out
#################
# a file for errors from the job
#SBATCH --error=myjob.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH --time=02:50:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB per CPU
#SBATCH --mem-per-cpu=4000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END,FAIL
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=ksb@stanford.edu
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks-per-node=16
#################

################## INPUT PARAMETERS ##############

#Convergence
ftol=.01 #eV, Difference between energy required for convergence
maxiter=500
x0 = [1,1,1] #initial guess for strain relative to initial structure

############ NO NEED TO CHANGE ANYTHING BELOW HERE #############
import sys,os,ase,shutil
import glob,subprocess
from ase.all import *
from scipy.optimize import *

path    = os.getcwd() + '/'
traj    = glob.glob('*.traj')[0]
initial = ase.io.read(traj)


os.makedirs('input')
shutil.copy(traj,'input')


initial_cell=initial.get_cell()
print "Initial cell: "
print initial_cell



def get_energy(strains):
    "Function to create an atoms object with lattice parameters specified by x and obtain its energy using calc."
    
    dirs = [f for f in os.listdir('.') if os.path.isdir(f)]
    for d in dirs:
        if str(d)[:5] == 'input' and len(str(d))>5: shutil.rmtree(d)

    sys.stdout.flush()
    traj = glob.glob('*.traj')[0]

    atoms=ase.io.read(traj)
    #Scale x,y,z coords by constant if bulk, scale only x,y coords if slab
    atoms.set_cell(np.multiply(initial_cell,strains),scale_atoms=True)

    ase.io.write(traj,atoms)
    #Relax cell with new cell lengths
    subprocess.call(['python', 'RELAX.py'])

    #If job converges, read energy from output file, store in logfile, return energy as output
    with open ("converged.txt", "r") as output:    energy = float(output.read())
    with open ('lattice_opt.log','a+') as logfile: logfile.write('%s\t%s\n' %(energy,strains))
    return energy

x = fmin(get_energy,x0,maxiter=maxiter,ftol=ftol) #this can be changed to minimize(...method = method...) in scipy > v0.1
#x = minimize(get_energy,x0)#,method='Nelder-Mead',bounds=(.5,1.5))

print "Optimized multiplier: "
print x
print "Final cell: "
print initial_cell*x
initial.set_cell(np.multiply(initial_cell,strains),scale_atoms=True)
ase.io.write(traj,initial)
