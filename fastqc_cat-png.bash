#!/usr/bin/env bash

# Aim: concatenate PNG files from FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) into a single PDF file
# Copyright (C) 2014 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: https://github.com/timflutre/quantgen

progVersion="1.1.0" # http://semver.org/

# Display the help on stdout.
# The format complies with help2man (http://www.gnu.org/s/help2man)
function help () {
  msg="\`${0##*/}' concatenates PNG files from FastQC into a single PDF file.\n"
  msg+="\n"
  msg+="Usage: ${0##*/} [OPTIONS] ...\n"
  msg+="\n"
  msg+="Options:\n"
  msg+="  -h, --help\tdisplay the help and exit\n"
  msg+="  -V, --version\toutput version information and exit\n"
  msg+="  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
  msg+="  -I, --zip\tglob for the ZIP archives\n"
  msg+="  -i, --png\tglob for the PNG files (if already extracted from the archives)\n"
  msg+="  -p, --pre\tprefix of the PNG file inside the ZIP archives\n"
  msg+="  -o, --out\tpath to the PDF file\n"
  msg+="  -n, --name\tadd file name to each PNG file\n"
  msg+="  -c, --clean\tremove PNG files\n"
  msg+="\n"
  msg+="Examples:\n"
  msg+="  fastqc -o ./ *.fastq.gz # run FastQC on several fastq files\n"
  msg+="  unzip -l lane1_fastqc.zip # list one ZIP archive to check file names\n"
  msg+="  ${0##*/} -I \"*_fastqc.zip\" -p per_base_quality -o plots_per_base_quality.pdf -c\n"
  msg+="\n"
  msg+="Report bugs to <timothee.flutre@supagro.inra.fr>."
  echo -e "$msg"
}

# Display version and license information on stdout.
# The person roles complies with R's guidelines (The R Journal Vol. 4/1, June 2012).
function version () {
  msg="${0##*/} ${progVersion}\n"
  msg+="\n"
  msg+="Copyright (C) 2014 Institut National de la Recherche Agronomique (INRA).\n"
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
	  TEMP=`getopt -o hVv:I:i:p:o:nc -l help,version,verbose:,zip:,png:,prefix:,out:,name,clean, \
        -n "$0" -- "$@"`
  else # original getopt is available (no long options, whitespace, sorting)
	  TEMP=`getopt hVv:I:i:p:o:nc "$@"`
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
      -I | --zip) globZip=$2; shift 2;;
      -i | --png) globPng=$2; shift 2;;
      -p | --pre) fastqcPrefix=$2; shift 2;;
      -o | --out) outPdf=$2; shift 2;;
      -n | --name) addFileName=true; shift;;
      -c | --clean) clean=true; shift;;
      --) shift; break;;
      *) echo "ERROR: options parsing failed, use -h for help" 1>&2; exit 1;;
    esac
  done

  hash convert 2>/dev/null || \
    { echo >&2 "ERROR: convert (from ImageMagick) is not in your PATH"; exit 1; }
  if [ "${addFileName}" == true ]; then
    hash gs 2>/dev/null || \
      { echo >&2 "ERROR: gs (from Ghostscript) is not in your PATH"; exit 1; }
  fi
  if [ -z "${globZip}" -a -z "${globPng}" ]; then
    echo -e "ERROR: missing compulsory option -I or -i\n" 1>&2
    help
    exit 1
  fi
  if [ -z "${fastqcPrefix}" ]; then
    echo -e "ERROR: missing compulsory option -p\n" 1>&2
    help
    exit 1
  fi
  if [ -z "${outPdf}" ]; then
    echo -e "ERROR: missing compulsory option -o\n" 1>&2
    help
    exit 1
  fi
}

function run () {
  if [ ! -z "${globZip}" ]; then
    if [ $verbose -gt "0" ]; then
      echo -n "extract '${fastqcPrefix}' PNG files"
      echo " from $(ls ${globZip} | wc -l) ZIP archives ..."
    fi
    ls ${globZip} | while read f; do \
      i=$(echo $f | awk '{split($0,a,"."); print a[1]}'); \
      unzip -p ${f} ${i}/Images/${fastqcPrefix}.png > ${i}_${fastqcPrefix}.png; \
      done
    globPng="*_${fastqcPrefix}.png"
  fi

  globPng2=""
  if [ "${addFileName}" == true ]; then
    if [ $verbose -gt "0" ]; then
      echo "add file name to each PNG file ..."
    fi
    for f in ${globPng}; do
      f2=$(basename "$f")
      f2="${f2%.*}"
      newSuffix="pdf"
      convert "${f}" -gravity south -annotate +0+100 "%f" "tmp_${f2}.${newSuffix}"
    done
    globPng2="tmp_*.${newSuffix}"
  fi

  if [ $verbose -gt "0" ]; then
    echo "concatenate $(ls ${globPng} | wc -l) PNG files into a single PDF file ..."
  fi
  rm -f ${outPdf}
  if [ "${addFileName}" == true ]; then
    # convert "tmp_*_${fastqcPrefix}.png" ${outPdf}
    gs -dNOPAUSE -dQUIET -sDEVICE=pdfwrite -dPDFSETTINGS=/screen \
      -sOUTPUTFILE=${outPdf} -dBATCH ${globPng2}
  else
    convert ${globPng} ${outPdf}
  fi

  if [ "${clean}" == true ]; then
    rm -f ${globPng} ${globPng2}
  fi
}

verbose=1
globZip=""
globPng=""
fastqcPrefix=""
outPdf=""
addFileName=false
clean=false
parseCmdLine "$@"

if [ $verbose -gt "0" ]; then
  startTime=$(timer)
  msg="START ${0##*/} ${progVersion} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
  # msg+="\ncmd-line: $0 "$@ # comment if an option takes a glob as argument
  msg+="\ncwd: $(pwd)"
  echo -e $msg
fi

run globZip globPng fastqcPrefix outPdf addFileName clean verbose

if [ $verbose -gt "0" ]; then
  msg="END ${0##*/} ${progVersion} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
  msg+=" ($(timer startTime))"
  echo $msg
fi
