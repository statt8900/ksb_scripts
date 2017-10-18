# .bashrc


#########################
### Environment Variables
#########################

export HISTCONTROL=ignoreboth:erasedups  # no duplicate entries
export PYTHONPATH=~/scripts/datalog/:~/scripts/data:~/scripts/misc:$PYTHONPATH
export PYTHONPATH=~/scripts/quantumespresso:~/scripts/wflow3:$PYTHONPATH
export PYTHONPATH=~/scripts/newJobCreation:~/scripts/plotting:$PYTHONPATH
export PYTHONPATH=~/scripts/fireworks:$PYTHONPATH

export PATH=~/Zeo++/zeo/trunk:~/Zeo++/zeo/trunk/cython_wrapper:~/Voro++/voro/trunk/src:$PATH
export PYTHONPATH=~/Zeo++/zeo/trunk:~/Zeo++/zeo/trunk/cython_wrapper:~/Voro++/voro/trunk/src:$PYTHONPATH
export PATH=/home/ksb/enumlib/src:$PATH

##################
# Quantum espresso
##################
export OMP_NUM_THREADS=1  
export NO_STOP_MESSAGE=yes #eliminates spurious warning from quantum espresso
export GPAW_SETUP_PATH=/scratch/users/ksb/gpaw/gpaw_sg15/norm_conserving_setups
export PATH=$HOME/stack-1.4.0-linux-x86_64-static:$PATH
export PATH=$SCRATCH/fireworks:$PATH


############
### Espresso
############
export PATH=$HOME/espresso/q-e/bin:$PATH
export PATH=$HOME/scripts/misc:$PATH
export PYTHONPATH=$HOME/espresso/lib/python2.7/site-packages:$PYTHONPATH
export PYTHONPATH=$HOME/espresso:$PYTHONPATH #is this superfluous?
export ESP_PSP_PATH=/scratch/PI/suncat/dacapopsp

#####
# ASE
#####
export PYTHONPATH=/scratch/users/ksb/asethree14/lib/python2.7/site-packages:$PYTHONPATH
export PATH=/scratch/users/ksb/asethree14/bin:$PATH

########
# CUSTOM
########
export FW_PATH=$SCRATCH/fireworks/jobs/
export ASE_PATH=$SCRATCH/db/ase.db
export ALL_FWS=$SCRATCH/fireworks/alljobs/
export DATA_PATH=$SCRATCH/db/data.db
export IMG_PATH=$SCRATCH/img/

###########
### Aliases
###########


# Terminal Navigation
alias cds='cd /scratch/users/ksb/'
alias fire='cd /scratch/users/ksb/fireworks'
alias sunc='ssh -Y ksb@suncatls1.slac.stanford.edu'

# SLURM
alias idle='sinfo |grep idle'
alias nrun='squeue -u ksb | grep ksb | grep -w R | wc -l'
alias nq='squeue -u ksb | grep ksb | wc -l'
alias sq='squeue -u ksb -o "%.18i %.9P %.18j %.8u %.2t %.10M %.6D %Z"'
alias sInfo='scontrol show job '
alias jobs='$HOME/scripts/misc/jobSHERLOCK.py'
# FIREWORKS
alias ready='lpad get_fws -s READY'
alias waiting='lpad get_fws -s WAITING'
alias reserved='lpad get_fws -s RESERVED'
alias fizzled='lpad get_fws -s FIZZLED'
alias archived='lpad get_fws -s ARCHIVED'
alias running='lpad get_fws -s RUNNING'
alias completed='lpad get_fws -s COMPLETED'
alias info='fire;lpad get_fws -i'
alias lpadc='lpad -l $SCRATCH/fireworks_queue_dir/my_launchpad.yaml'

alias launcher='~/scripts/fireworks/launcher.sh'
alias maintain='~/scripts/fireworks/daemon_sherlock.sh'
#Misc
alias ewh='/home/ksb/scripts/misc/ewh.py'
alias img='/home/ksb/scripts/misc/imgcat.sh'
alias piScratch='source $PI_SCRATCH/sw/env.bash' 
alias upgradeASE='export PYTHONPATH=/scratch/users/ksb/ase-3.14.1:$PYTHONPATH'

#Function shorthand
alias mod='python -m'
alias vi='vim -p'
alias mem='du -sh *'
alias ltr='ll -tr'
alias count='tree -L 1 | tail -1'

alias agneb="ase gui initial.traj neb*.traj final.traj -n -1"
alias ag='ase gui '
alias actualRM='/bin/rm'

alias all='/home/ksb/scripts/misc/all.sh'
alias cpall='/home/ksb/scripts/misc/cpall.sh'
alias rmx='/home/ksb/scripts/misc/rmx.sh'

###########
# Functions
##########

py() {
    python -c "from sqlShortcuts import *;import $1;$1.$2"
}
py2() {
    python -c "from sqlShortcuts import *;from $1 import *;from $2 import *;$3"
}

pyprint() {
    python -c "from sqlShortcuts import *;import $1;print $1.$2"
}


sub() {
    eval 'sbatch -J $PWD "$@"'
}

# File transfer
tosunc ()
{
    DIR=${2:-/afs/slac.stanford.edu/u/if/ksb/scp/};
    rsync -C -a --stats $1 ksb@suncatls1.slac.stanford.edu:/$DIR
}

fromsunc ()
{
    DIR=${2:-/home/ksb/scp_tmp};
    rsync -C -a ksb@suncatls1.slac.stanford.edu:$1 $DIR
}

tocori() {
   DIR=${2:-/global/u1/k/krisb/scp_tmp/}
    scp -C -r $1 krisb@cori.nersc.gov:/$DIR
}

fromcori() {
    DIR=${2:-/home/ksb/scp/}
    scp -C -r krisb@cori.nersc.gov:$1 $DIR
}

auto_update() {
      processes=$(ps -ef | grep fswatch | wc -l | xargs echo)
      echo $processes
      if [ ! $processes -gt 1 ]
      then
          fswatch -o /scratch/users/ksb/wflow/ | xargs -n1 /bin/sh -c /home/ksb/wflow/updateSuncatDB.sh &
      fi
  }

# LS 
LS_COLORS="rs=0:di=01;34:ln=01;36:mh=00:pi=40;33:so=01;35:do=01;35:bd=40;33;01:cd=40;33;01:or=40;31;01:mi=01;05;37;41:su=37;41:sg=30;43:ca=30;41:tw=30;42:ow=37;42:st=37;44:ex=01;32:*.tar=01;31:*.tgz=01;31:*.arj=01;31:*.taz=01;31:*.lzh=01;31:*.lzma=01;31:*.tlz=01;31:*.txz=01;31:*.zip=01;31:*.z=01;31:*.Z=01;31:*.dz=01;31:*.gz=01;31:*.lz=01;31:*.xz=01;31:*.bz2=01;31:*.tbz=01;31:*.tbz2=01;31:*.bz=01;31:*.tz=01;31:*.deb=01;31:*.rpm=01;31:*.jar=01;31:*.rar=01;31:*.ace=01;31:*.zoo=01;31:*.cpio=01;31:*.7z=01;31:*.rz=01;31:*.jpg=01;35:*.jpeg=01;35:*.gif=01;35:*.bmp=01;35:*.pbm=01;35:*.pgm=01;35:*.ppm=01;35:*.tga=01;35:*.xbm=01;35:*.xpm=01;35:*.tif=01;35:*.tiff=01;35:*.png=01;35:*.svg=01;35:*.svgz=01;35:*.mng=01;35:*.pcx=01;35:*.mov=01;35:*.mpg=01;35:*.mpeg=01;35:*.m2v=01;35:*.mkv=01;35:*.ogm=01;35:*.mp4=01;35:*.m4v=01;35:*.mp4v=01;35:*.vob=01;35:*.qt=01;35:*.nuv=01;35:*.wmv=01;35:*.asf=01;35:*.rm=01;35:*.rmvb=01;35:*.flc=01;35:*.avi=01;35:*.fli=01;35:*.flv=01;35:*.gl=01;35:*.dl=01;35:*.xcf=01;35:*.xwd=01;35:*.yuv=01;35:*.cgm=01;5:*.emf=01;35:*.axv=01;35:*.anx=01;35:*.ogv=01;35:*.ogx=01;35:*.aac=01;36:*.au=01;36:*.flac=01;36:*.mid=01;36:*.midi=01;36:*.mka=01;36:*.mp3=01;36:*.mpc=01;36:*.ogg=01;36:*.ra=01;36:*.wav=01;36:*.axa=01;36:*.oga=01;36:*.spx=01;36:*.xspf=01;36:*.py=01;91:*.txt=01;37:*.traj=01;33:*.cif=01;33:*.err=01;41:"
export LS_COLORS
# Reverse history search
if [ -t 1 ]
then

    bind '"\e[A": history-search-backward'
    bind '"\e[B": history-search-forward'
fi

source $PI_SCRATCH/sw/env.bash

######
######
######
######
