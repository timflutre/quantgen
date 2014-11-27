#!/usr/bin/env bash

# Aim: does this and that
# ---
# choose between:
# Not copyrighted -- provided to the public domain
# Author: Timothée Flutre
# or:
# Copyright (C) 2011-2013 Timothée Flutre
# License: GPL-3+
# Author: Timothée Flutre
# ---
# Versioning: https://github.com/timflutre/...

progVersion="1.0.0" # http://semver.org/

# Display the help on stdout.
# The format complies with help2man (http://www.gnu.org/s/help2man)
function help () {
    msg="\`${0##*/}' does this and that.\n"
    msg+="\n"
    msg+="Usage: ${0##*/} [OPTIONS] ...\n"
    msg+="\n"
    msg+="Options:\n"
    msg+="  -h, --help\tdisplay the help and exit\n"
    msg+="  -V, --version\toutput version information and exit\n"
    msg+="  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
    msg+="  -i, --input\tpath to the input file\n"
    msg+="\n"
    msg+="Examples:\n"
    msg+="  ${0##*/} -i <input>\n"
    msg+="\n"
    msg+="Report bugs to <>."
    echo -e "$msg"
}

# Display version and license information on stdout.
function version () {
    msg="${0##*/} ${progVersion}\n"
    msg+="\n"
    msg+="Copyright (C) 2011-2013 Timothée Flutre.\n"
    msg+="License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
    msg+="This is free software; see the source for copying conditions.  There is NO\n"
    msg+="warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
    msg+="\n"
    msg+="Written by Timothée Flutre."
    echo -e "$msg"
# or choose "Not copyrighted -- provided to the public domain\n"
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
	TEMP=`getopt -o hVv:i: -l help,version,verbose:,input: \
        -n "$0" -- "$@"`
    else # original getopt is available (no long options, whitespace, sorting)
	TEMP=`getopt hVv:i: "$@"`
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
            -i | --input) inFile=$2; shift 2;;
            --) shift; break;;
            *) echo "ERROR: options parsing failed, use -h for help" 1>&2; exit 1;;
        esac
    done
    if [ -z "${inFile}" ]; then
        echo -e "ERROR: missing compulsory option --input\n" 1>&2
        help
        exit 1
    fi
    if [ ! -f "${inFile}" ]; then
        echo -e "ERROR: can't find file ${inFile}\n" 1>&2
        help
        exit 1
    fi
}

function run () {
    
    # specific code ...
    true # to avoid an error when the function is empty: http://stackoverflow.com/a/2421637/597069
    
}

verbose=1
inFile=""
parseCmdLine "$@"

if [ $verbose -gt "0" ]; then
    startTime=$(timer)
    msg="START ${0##*/} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
    msg+="\ncmd-line: $0 "$@ # comment if an option takes a glob as argument
    msg+="\ncwd: $(pwd)"
    echo -e $msg
fi

run inFile verbose

if [ $verbose -gt "0" ]; then
    msg="END ${0##*/} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
    msg+=" ($(timer startTime))"
    echo $msg
fi
