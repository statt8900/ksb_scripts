#
# this file gets executed by every shell you start
#

setenv PRINTER InvalidPrinter
setenv WORK /nfs/slac/g/suncatfs/ksb
setenv SCRATCH $WORK
setenv CLUSTER slac
setenv LD_LIBRARY_PATH /nfs/slac/g/suncatfs/sw/espv17/intellib:/nfs/slac/g/suncatfs/sw/external/openmpi/1.6.3/install/lib:/nfs/slac/g/suncatfs/sw/vasp/v7/intellibs:/nfs/slac/g/suncatfs2/sw2/openmpi1.8.5rhel5/lib

setenv PATH /u/if/ksb/scripts/misc:$PATH

setenv PATH /nfs/slac/g/suncatfs/ksb/q-e/bin:$PATH
setenv PYTHONPATH /u/if/ksb/ase-qe/lib/python2.7/site-packages:$PYTHONPATH
setenv ESP_PSP_PATH /nfs/slac/g/suncatfs/sw/external/esp-psp/gbrv1.5pbe

setenv PYTHONPATH ~/scripts/auto/:$PYTHONPATH
setenv PYTHONPATH /u/if/ksb/scripts/data:$PYTHONPATH
setenv PYTHONPATH /u/if/ksb/scripts/datalog:$PYTHONPATH
setenv PYTHONPATH /u/if/ksb/scripts/misc:$PYTHONPATH
setenv PYTHONPATH /u/if/ksb/scripts/quantumespresso:$PYTHONPATH
setenv PYTHONPATH /u/if/ksb/scripts/wflow3:$PYTHONPATH
setenv PYTHONPATH /u/if/ksb/scripts/fireworks:$PYTHONPATH

setenv FW_PATH /nfs/slac/g/suncatfs/ksb/fireworks/jobs

# Set the prompt with ueful color scheme
set prompt="%{\033[1;31m%}%n@suncat%{\033[1;34m%} %~ >%{\033[0m%} "

if( ${?prompt} ) then
  # put aliases and other things down here
  alias ls 'ls --color=auto'
  alias e 'emacs'
  alias ebash 'emacs ~/.cshrc'
  alias sbash 'source ~/.cshrc'
  alias cls 'clear;ls'
  alias wd 'cd $WORK'
  alias cds 'cd $WORK'
  alias ag 'ase gui'
  alias grep 'grep --color=auto'
  alias jobs '/u/if/ksb/scripts/misc/jobsSLAC.py'
  alias launcher 'ps -ef | grep rapidfire'
  alias testSub '/u/if/ksb/scripts/rc/testSub.sh'
  alias sub110 '/u/if/ksb/scripts/rc/sub110.sh'
  alias	sub120 '/u/if/ksb/scripts/rc/sub120.sh'
  alias	sub210 '/u/if/ksb/scripts/rc/sub210.sh'
  alias	sub220 '/u/if/ksb/scripts/rc/sub220.sh'
  alias ready 'cd $SCRATCH/fireworks;source startFireworksSUNCAT.sh;$SCRATCH/fireworks/fireworks_virtualenv/bin/lpad get_fws -s READY'
  alias killpend 'bkill ` bjobs -u ksb |grep PEND |cut -f1 -d" "`'
  alias all '/u/if/ksb/scripts/misc/all.sh'
  alias cpall '/u/if/ksb/scripts/misc/cpall.sh'
  alias rmx '/u/if/ksb/scripts/misc/rmx.sh'
  alias hori '/u/if/ksb/scripts/misc/get_hori.py'
endif


