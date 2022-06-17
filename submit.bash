#!/usr/bin/env bash

# Aim: submit a simple job to HTCondor (https://htcondor.org/)
# Copyright (C) 2022 INRAE
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: https://github.com/timflutre/quantgen

progVersion="0.1.0" # http://semver.org/

# Display the help on stdout.
# The format complies with help2man (http://www.gnu.org/s/help2man)
function help () {
  msg="\`${0##*/}' submits a simple job to HTCondor (https://htcondor.org/).\n"
  msg+="\n"
  msg+="Usage: ${0##*/} [OPTIONS] ...\n"
  msg+="\n"
  msg+="Options:\n"
  msg+="  -h, --help\tdisplay the help and exit\n"
  msg+="  -V, --version\toutput version information and exit\n"
  msg+="  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
  msg+="  -e, --exe\tpath to the executable\n"
  msg+="  -a, --args\targuments to the executable (optional)\n"
  msg+="  -o, --out\tprefix of output files for the job (default=${out})\n"
  msg+="  -m, --mem\trequired memory (default=${mem})\n"
  msg+="  -c, --cpu\trequired CPUs (default=${cpu})\n"
  msg+="  -j, --job\tfile name for the job (default=${job})\n"
  msg+="\n"
  msg+="Examples:\n"
  msg+="  ${0##*/} -e '/bin/echo' -a 'Hello, world'\n"
  msg+="\n"
  msg+="Report bugs to <timothee.flutre@inrae.fr>."
  echo -e "$msg"
}

# Display version and license information on stdout.
# The person roles complies with R's guidelines (The R Journal Vol. 4/1, June 2012).
function version () {
  msg="${0##*/} ${progVersion}\n"
  msg+="\n"
  msg+="Copyright (C) 2022 INRAE.\n"
  msg+="License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
  msg+="This is free software; see the source for copying conditions.  There is NO\n"
  msg+="warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
  msg+="\n"
  msg+="Written by Timothée Flutre [cre,aut]."
  echo -e "$msg"
}

# http://www.linuxjournal.com/content/use-date-command-measure-elapsed-time
function timer () {
  if [[ $# -eq 0 ]]; then
    echo $(date '+%s')
  else
    local startRawTime=$1
    endRawTime=$(date '+%s')
    if [[ -z "$startRawTime" ]]; then startRawTime=$endRawTime; fi
    elapsed=$((endRawTime - startRawTime)) # in sec
    nbDays=$((elapsed / 86400))
    nbHours=$(((elapsed / 3600) % 24))
    nbMins=$(((elapsed / 60) % 60))
    nbSecs=$((elapsed % 60))
    printf "%01dd %01dh %01dm %01ds" $nbDays $nbHours $nbMins $nbSecs
  fi
}

# Parse the command-line arguments.
# http://stackoverflow.com/a/4300224/597069
function parseCmdLine () {
  getopt -T > /dev/null # portability check (say, Linux or Mac OS?)
  if [ $? -eq 4 ]; then # GNU enhanced getopt is available
	  TEMP=`getopt -o hVv:e:a:o:m:c:j: -l help,version,verbose:,exe:,args:out:mem:cpu:job:, \
        -n "$0" -- "$@"`
  else # original getopt is available (no long options, whitespace, sorting)
	  TEMP=`getopt hVv:e:a:o:m:c:j: "$@"`
  fi
  if [ $? -ne 0 ]; then
	  echo "ERROR: "$(which getopt)" failed" 1>&2
	  getopt -T > /dev/null
	  if [ $? -ne 4 ]; then
	    echo "did you use long options? they are not handled \
on your system, use -h for help"
	  fi
	  exit 2
  fi
  eval set -- "$TEMP"
  while [ $# -gt 0 ]; do
    case "$1" in
      -h | --help) help; exit 0; shift;;
      -V | --version) version; exit 0; shift;;
      -v | --verbose) verbose=$2; shift 2;;
      -e | --exe) exe=$2; shift 2;;
      -a | --args) args=$2; shift 2;;
      -o | --out) out=$2; shift 2;;
      -m | --mem) mem=$2; shift 2;;
      -c | --cpu) cpu=$2; shift 2;;
      -j | --job) job=$2; shift 2;;
      --) shift; break;;
      *) echo "ERROR: options parsing failed, use -h for help" 1>&2; exit 1;;
    esac
  done

  hash condor_submit 2>/dev/null || \
    { echo >&2 "ERROR: condor_submit is not in your PATH"; exit 1; }
  if [ -z "${exe}" ]; then
    echo -e "ERROR: missing compulsory option -e\n" 1>&2
    help
    exit 1
  fi
}

function run () {
  if [ $verbose -gt "0" ]; then
    msg="write config file"
    echo -e $msg
  fi
  txt="Universe = vanilla"
  txt+="\nExecutable = "${exe}
  if [ ! -z "${args}" ]; then
    txt+="\nArguments = "${args}
  fi
  txt+="\n"
  txt+="\nshould_transfer_files = no"
  txt+="\n"
  txt+="\ninput = /dev/null"
  txt+="\noutput = ${out}.o\$(Cluster)"
  txt+="\nerror = ${out}.e\$(Cluster)"
  txt+="\nlog = ${out}.l\$(Cluster)"
  txt+="\n"
  txt+="\nrequest_memory = ${mem}"
  txt+="\nrequest_cpus = ${cpu}"
  txt+="\nJobLeaseDuration = 30"
  txt+="\nrequirements = ( HAS_ASREML =?= False )"
  txt+="\ngetenv = true"
  txt+="\n"
  txt+="\nQueue"
  txt+="\n"
  echo -e ${txt} > ${job}

  if [ $verbose -gt "0" ]; then
    msg="submit job"
    echo -e $msg
  fi
  condor_submit ${job}
}

verbose=0
exe=""
args=""
out="out_condor"
mem="4G"
cpu="1"
job="job_file"
parseCmdLine "$@"

if [ $verbose -gt "0" ]; then
  startTime=$(timer)
  msg="START ${0##*/} ${progVersion} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
  # msg+="\ncmd-line: $0 "$@ # comment if an option takes a glob as argument
  msg+="\ncwd: $(pwd)"
  echo -e $msg
fi

run exe args

if [ $verbose -gt "0" ]; then
  msg="END ${0##*/} ${progVersion} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
  msg+=" ($(timer startTime))"
  echo $msg
fi
