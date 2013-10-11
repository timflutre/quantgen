#!/usr/bin/env bash

# Aim: does this and that
# choose between:
# Author: Timothee Flutre
# Not copyrighted -- provided to the public domain
# or:
# Copyright (C) 2011-2013 Timothee Flutre
# License: GPLv3+

function help () {
    msg="\`${0##*/}' does this and that.\n"
    msg+="\n"
    msg+="Usage: ${0##*/} [OPTIONS] ...\n"
    msg+="\n"
    msg+="Options:\n"
    msg+="  -h, --help\tdisplay the help and exit\n"
    msg+="  -V, --version\toutput version information and exit\n"
    msg+="  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
    msg+="  -i, --in\tinput\n"
    msg+="\n"
    msg+="Examples:\n"
    msg+="  ${0##*/} -i <input>\n"
    echo -e "$msg"
}

function version () {
    msg="${0##*/} 1.0\n"
    msg+="\n"
    msg+="Written by Timothee Flutre.\n"
    msg+="\n"
# choose between:
    msg += "Not copyrighted -- provided to the public domain\n"
# or:
    msg+="Copyright (C) 2011-2013 Timothee Flutre.\n"
    msg+="License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
    msg+="This is free software; see the source for copying conditions.  There is NO\n"
    msg+="warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
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

function parseArgs () {
    TEMP=`getopt -o hVv:i: -l help,version,verbose:,in: \
        -n "$0" -- "$@"`
    if [ $? != 0 ] ; then echo "ERROR: getopt failed" >&2 ; exit 1 ; fi
    eval set -- "$TEMP"
    while true; do
        case "$1" in
            -h|--help) help; exit 0; shift;;
            -V|--version) version; exit 0; shift;;
            -v|--verbose) verbose=$2; shift 2;;
            -i|--in) input=$2; shift 2;;
            --) shift; break;;
            *) echo "ERROR: options parsing failed"; exit 1;;
        esac
    done
    if [ -z "${input}" ]; then
        echo -e "ERROR: missing compulsory option -i\n"
        help
        exit 1
    fi
    if [ ! -f "${input}" ]; then
        echo -e "ERROR: can't find '${input}'\n"
        help
        exit 1
    fi
}

verbose=1
input=""
parseArgs "$@"

if [ $verbose -gt "0" ]; then
    startTime=$(timer)
    msg="START ${0##*/} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
    msg+="\ncmd-line: $0 "$@ # comment if an option takes a glob as argument
    echo -e $msg
fi

# ... specific code ...

if [ $verbose -gt "0" ]; then
    msg="END ${0##*/} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
    msg+=" ($(timer startTime))"
    echo $msg
fi
