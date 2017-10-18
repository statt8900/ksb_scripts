#!/bin/bash

# Prints a Slurm cluster status with 1 line per node/partition and job info
# Author: Ole.H.Nielsen@fysik.dtu.dk
# Homepage: https://github.com/OleHolmNielsen/Slurm_tools/tree/master/pestat
my_version="pestat version 0.61. Date: 07 Jul 2017"

#####################################################################################
#
# Environment configuration lines
#

# CONFIGURE the paths of commands:
# Directory where Slurm commands live:
export prefix=/usr/bin
# The awk command
my_awk=$prefix/awk
# Variable width of the hostname column (default=8)
export hostnamelength="8"
# Global configuration file for pestat
export PESTAT_GLOBAL_CONFIG=/etc/pestat.conf
# Per-user configuration file for pestat
export PESTAT_CONFIG=$HOME/.pestat.conf
# End of CONFIGURE lines

#####################################################################################
#
# Command usage:
#
function usage()
{
	cat <<EOF
Usage: pestat [-p partition(s)] [-u username] [-g groupname]
	[-q qoslist] [-s statelist] [-n/-w hostlist] [-j joblist]
	[-f | -F | -m free_mem | -M free_mem ] [-1] [-E] [-C/-c] [-V] [-h]
where:
	-p partition: Select only partion <partition>
	-u username: Print only user <username> 
	-g groupname: Print only users in UNIX group <groupname>
	-q qoslist: Print only QOS in the qoslist <qoslist>
	-s statelist: Print only nodes with state in <statelist> 
	-n/-w hostlist: Print only nodes in hostlist
	-j joblist: Print only nodes in job <joblist>
	-f: Print only nodes that are flagged by * (unexpected load etc.)
	-F: Like -f, but only nodes flagged in RED are printed.
	-m free_mem: Print only nodes with free memory LESS than free_mem MB
	-M free_mem: Print only nodes with free memory GREATER than free_mem MB (under-utilized)
	-1: Only 1 line per node (unique nodes in multiple partitions are printed once only)
	-E: Job EndTime is printed after each jobid/user
	-C: Color output is forced ON
	-c: Color output is forced OFF
	-h: Print this help information
	-V: Version information

Global configuration file for pestat: $PESTAT_GLOBAL_CONFIG
Per-user configuration file for pestat: $PESTAT_CONFIG
EOF
}

#####################################################################################
#
# Configuration of default pestat parameters.
# Default values are overridden by $PESTAT_GLOBAL_CONFIG and $PESTAT_CONFIG
#

# Flagging nodes
export flaglevel=-1
# Free memory threshold value (off by default)
export free_mem_less=-1
export free_mem_more=-1

# Thresholds for flagging nodes
export cpuload_delta1=2.0	# CPU load delta from ideal load (RED)
export cpuload_delta2=0.5	# CPU load delta from ideal load (MAGENTA)
export memory_thres1=0.1	# Fraction of memory which is free (RED)
export memory_thres2=0.2	# Fraction of memory which is free (MAGENTA)

# Colored output by default:
export colors=1

# Print all nodes by default
export printallnodes=1
# Unique nodes only (nodes in multiple partitions are printed once only)
export uniquenodes=0

# Controls whether job EndTime will be printed
export printjobendtime=0

# Select UNIX group
export groupname=""

export partition=""
export hostlist=""

# Check if output is NOT a terminal: Turn colors off (can be overruled by "-C" flag).
FD=1	# File Descriptor no. 1 = stdout
if test ! -t $FD
then
	export colors=0
fi

# First read the global configuration file for pestat
if test -s $PESTAT_GLOBAL_CONFIG
then
	source $PESTAT_GLOBAL_CONFIG
fi

# Next read the per-user configuration file for pestat
if test -s $PESTAT_CONFIG
then
	source $PESTAT_CONFIG
fi

#####################################################################################
#
# Parse command options
#

while getopts "p:u:g:q:s:n:j:w:m:M:CchVFf1E" options; do
	case $options in
		p )	export partition="-p $OPTARG"
			echo Print only nodes in partition $OPTARG
			;;
		u )	export username=$OPTARG
			if test "`$prefix/sacctmgr -p -n show assoc where users=$username`"
			then
				echo Select only user $username
			else
				echo Error selecting Slurm username $username 
				exit -1
			fi
			;;
		g )	export groupname="$OPTARG"
			if test "`getent group $OPTARG`"
			then
				echo Print only users in UNIX group $OPTARG
			else
				echo Error selecting UNIX group $OPTARG
				exit -1
			fi
			;;
		q )	export qoslist=$OPTARG
			if test "`$prefix/sacctmgr -n -p show qos $qoslist`"
			then
				echo Select only QOS=$qoslist
			else
				echo Error selecting QOS $qoslist
				echo Print all available QOSes by: sacctmgr show qos
				exit -1
			fi
			;;
		s )	statelist="--states $OPTARG"
			echo Select only nodes with state=$OPTARG
			;;
		n|w )	hostlist="-n $OPTARG"
			echo Select only nodes in hostlist=$OPTARG
			;;
		j )	hostlist="-n `$prefix/squeue -j $OPTARG -h -o %N`"
			echo Select only nodes with jobs in joblist=$OPTARG
			;;
		f|F )	export flaglevel=1
			flag_color="(RED and MAGENTA nodes)"
			if test $options = "F"
			then
				export flaglevel=2	# Flag RED nodes only
				flag_color="(RED nodes)"
			fi
			if test $free_mem_less -ge 0 -o $free_mem_more -ge 0
			then
				echo ERROR: The -f -m -M flags are mutually exclusive
				exit -1
			fi
			export printallnodes=0
			echo Print only nodes that are flagged by \* $flag_color
			;;
		m )	export free_mem_less=$OPTARG
			if test $flaglevel -gt 0 -o $free_mem_more -ge 0
			then
				echo ERROR: The -f -m -M flags are mutually exclusive
				exit -1
			fi
			export printallnodes=0
			echo Select only nodes with free memory LESS than $free_mem_less MB
			;;
		M )	export free_mem_more=$OPTARG
			if test $flaglevel -gt 0 -o $free_mem_less -ge 0
			then
				echo ERROR: The -f -m -M flags are mutually exclusive
				exit -1
			fi
			export printallnodes=0
			echo Select only nodes with free memory GREATER than $free_mem_more MB
			;;
		1 )	export uniquenodes=1
			echo "Only 1 line per node (Unique nodes in multiple partitions are printed once only)"
			;;
		E )	export printjobendtime=1
			echo "Job EndTime is printed after each jobid/user"
			;;
		C )	export colors=1
			echo Force colors ON in output
			;;
		c )	export colors=0
			echo Force colors OFF in output
			;;
		V ) echo $my_version
			exit 1;;
		h|? ) usage
			exit 1;;
		* ) usage
			exit 1;;
	esac
done

# Test for extraneous command line arguments
if test $# -gt $(($OPTIND-1))
then
	echo ERROR: Too many command line arguments: $*
	usage
	exit 1
fi

#####################################################################################
#
# Main pestat function: Execute Slurm commands and print output.
#

# Print all nodes: NODELIST PARTITION CPU CPU_LOAD MEMORY FREE_MEM STATE
$prefix/sinfo -h -N $partition $hostlist $statelist -o "%N %P %C %O %m %e %t" | $my_awk '
BEGIN {
	#####################################################################################
	# Initialization

	# Read the environment variables configuring actions of pestat:
	prefix=ENVIRON["prefix"]
	username=ENVIRON["username"]
	groupname=ENVIRON["groupname"]
	qoslist=ENVIRON["qoslist"]
	free_mem_less=ENVIRON["free_mem_less"]
	free_mem_more=ENVIRON["free_mem_more"]
	cpuload_delta1=ENVIRON["cpuload_delta1"]
	cpuload_delta2=ENVIRON["cpuload_delta2"]
	memory_thres1=ENVIRON["memory_thres1"]
	memory_thres2=ENVIRON["memory_thres2"]
	printallnodes=ENVIRON["printallnodes"]
	printjobendtime=ENVIRON["printjobendtime"]
	uniquenodes=ENVIRON["uniquenodes"]
	flaglevel=ENVIRON["flaglevel"]
	hostnamelength=ENVIRON["hostnamelength"]
	colors=ENVIRON["colors"]
	# Define terminal colors for the output if requested
	if (colors != 0) {
		# See http://en.wikipedia.org/wiki/ANSI_escape_code#Colors
		RED="\033[1;4;31m"
		GREEN="\033[1;32m"
		MAGENTA="\033[1;35m"
		NORMAL="\033[0m"
	}

	if (qoslist != "") selection = selection " -q " qoslist
	# The "scontrol show hostnames" command is used to expand NodeList expressions
	HOSTLIST=prefix "/scontrol show hostnames "

	#####################################################################################
	# Gather the list of running jobs with squeue

	# Running jobs info: JobState JobId User group NodeList EndTime
	JOBLIST = prefix "/squeue -t RUNNING -h -o \"%T %A %u %g %N %e\" " selection
	while ((JOBLIST | getline) > 0) {
		JobState=$1
		# Replaced by -t RUNNING flag: Skip jobs if not in RUNNING state
		# if (JobState != "RUNNING") continue
		JobId=$2
		User=$3
		Group=$4
		NodeList=$5
		EndTime=$6
		# Select job information to be printed for this job
		JobInfo = JobId " " User " "
		if (printjobendtime == 1) JobInfo = JobInfo EndTime " "
		# May select a UNIX group name
		if (groupname != "" && Group != groupname) continue
		# Create the list of nodes for this job, expand the list if necessary
		if (index(NodeList,"[") == 0) {
			jobnodes[1] = NodeList
		} else {
			# Put hostname lines into an array jobnodes[]
			cmd = HOSTLIST NodeList
			i=0
			while ((cmd | getline) > 0) jobnodes[++i] = $1
			close (cmd)
		}
		# Populate the node arrays with "JobId User" (multiple jobs may exist, EndTime may be added)
		for (i in jobnodes) {
			n = jobnodes[i]
			hostname[n] = n
			jobs[n] = jobs[n] JobInfo
			numjobs[n]++
			# If username has been selected and node "n" runs job belonging to username:
			if (User == username) selecteduser[n] = User
			if (Group == groupname) selectedgroup[n] = Group
		}
		delete jobnodes
	}
	close (JOBLIST)
	# The job information may include EndTime
	JobInfo = "JobId User"
	if (printjobendtime == 1) JobInfo = JobInfo " EndTime"
	JobInfo = JobInfo " ..."
	# Print a header line (variable hostnamelength)
	printf("%*s %15s %8s %7s %8s %8s %8s  %s\n", hostnamelength, "Hostname", "Partition", "Node", "Num_CPU", "CPUload", "Memsize", "Freemem", "Joblist")
	printf("%*s %15s %8s %7s %8s %8s %8s  %s\n", hostnamelength, "", "", "State", "Use/Tot", "", "(MB)", "(MB)", JobInfo)
}
{
	#####################################################################################
	# Main section: Process lines from sinfo

	node=$1
	# Selection of subset of nodes
	if (selection != "" && jobs[node] == "") next
	if (username != "" && selecteduser[node] == "") next
	if (groupname != "" && selectedgroup[node] == "") next

	partition=$2
	# sinfo -o %C gives number of CPUs by state in the format "allocated/idle/other/total"
	split($3,cpulist,"/")
	cpuload=$4
	memory=$5
	freemem=$6
	state=$7
	gsub("*", "", state)	# Strip "*" from nodename in default partition

	# Select only subset of nodes with certain values/states
	listnode = printallnodes

	if (free_mem_less > 0) {
		# Free memory on the node LESS than free_mem_less
		if (freemem < free_mem_less) listnode++
	} else if (free_mem_more > 0) {
		# Free memory on the node GREATER than free_mem_more
		if (freemem > free_mem_more) listnode++
	} else {
		if (state == "drain" || state == "drng" || state == "resv" || state == "down" || state == "error") {
			# Flag nodes with status down, drain etc.
			stateflag="*"
			statecolor=RED
			listnode++
		} else {
			stateflag=" "
			statecolor=NORMAL
		}
		# Flag unexpected CPU load average
		loaddiff = cpuload - cpulist[1]
		if (loaddiff > cpuload_delta1 || loaddiff < -cpuload_delta1) {
			loadflag="*"
			loadcolor=RED
			cpucolor=GREEN
			listnode++
		} else if (loaddiff > cpuload_delta2 || loaddiff < -cpuload_delta2) {
			loadflag="*"
			loadcolor=MAGENTA
			cpucolor=GREEN
			if (flaglevel == 1) listnode++
		} else {
			loadflag=" "
			loadcolor=NORMAL
			cpucolor=NORMAL
		}
 		# Flag unexpected number of jobs
		if (numjobs[node] > cpulist[1]) {	# Should be at least 1 task per job
			jobflag="*"
			jobcolor=RED
			listnode++
		} else {
			jobflag=" "
			jobcolor=NORMAL
		}

		# Free memory on the node
		if (freemem < memory*memory_thres1) {	# Very high memory usage (RED)
			memflag="*"
			freememcolor=RED
			memcolor=GREEN
			listnode++
		} else if (freemem < memory*memory_thres2) {	# High memory usage (MAGENTA)
			memflag="*"
			freememcolor=MAGENTA
			memcolor=GREEN
			if (flaglevel == 1) listnode++
		} else {
			memflag=" "
			freememcolor=NORMAL
			memcolor=NORMAL
		}
	}
		
	if (listnode > 0 && nodewasprinted[node] == 0) {
		if (uniquenodes > 0) nodewasprinted[node]++	# Count this node (may be in multiple partitions)
		printf("%*s %15s ", hostnamelength, node, partition)
		printf("%s%8s%1s%s", statecolor, state, stateflag, NORMAL)
		printf("%s%3d%s %3d ", cpucolor, cpulist[1], NORMAL, cpulist[4])
		printf("%s%7.2f%1s%s ", loadcolor, cpuload, loadflag, NORMAL)
		printf("%s%8d %s%8d%1s%s ", memcolor, memory, freememcolor, freemem, memflag, NORMAL)
		printf("%s%s%1s%s", jobcolor, jobs[node], jobflag, NORMAL)
		printf("\n")
	}
	delete cpulist
}'

