#!/usr/bin/env Rscript

## Aim: does this and that
## choose between:
## Author: Timothee Flutre
## Not copyrighted -- provided to the public domain
## or:
## Copyright (C) 2011-2013 Timothee Flutre
## License: GPLv3+

rm(list=ls())
prog.name <- "myprogram.R"
prog.version <- "1.0"

## Display the help on stdout.
## The format complies with help2man (http://www.gnu.org/s/help2man)
help <- function(){
  txt <- paste0("`", prog.name, "' does this and that.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Usage: ", prog.name, " [OPTIONS] ...\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Options:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "  -i, --input\tpath to the input file\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Examples:\n")
  txt <- paste0(txt, "  ", prog.name, " -i <input>\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Remarks:\n")
  txt <- paste0(txt, "  This is my typical template file for R.")
  message(txt)
}

## Display version and license information on stdout.
version <- function(){
  txt <- paste0(prog.name, " ", prog.version, "\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Written by Timothee Flutre.\n")
  txt <- paste0(txt, "\n")
  ## choose between:
  txt <- paste0(txt, "Not copyrighted -- provided to the public domain\n")
  ## or:
  txt <- paste0(txt, "Copyright (C) 2011-2013 Timothee Flutre.\n")
  txt <- paste0(txt, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n")
  txt <- paste0(txt, "This is free software; see the source for copying conditions.  There is NO\n")
  txt <- paste0(txt, "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n")
  message(txt)
}

## Parse the command-line arguments.
parseCmdLine <- function(params){
  args <- commandArgs(trailingOnly=TRUE)
  ## print(args)
  
  i <- 0
  while(i < length(args)){ # use "while" loop for options with no argument
    i <- i + 1
    if(args[i] == "-h" || args[i] == "--help"){
      help()
      quit("no", status=0)
    }
    else if(args[i] == "-V" || args[i] == "--version"){
      version()
      quit("no", status=0)
    }
    else if(args[i] == "-v" || args[i] == "--verbose"){
      params$verbose <- as.numeric(args[i+1])
      i <- i + 1
    }
    else if(args[i] == "-i" || args[i] == "--input"){
      params$in.file <- args[i+1]
      i <- i + 1
    }
    else{
      write(paste0(prog.name, ": invalid option -- ", args[i], "\n"), stderr())
      help()
      quit("no", status=1)
    }
  }
  
  return(params)
}

## Check the values of the command-line parameters.
checkParams <- function(params){
  if(is.null(params$in.file)){
    write("ERROR: missing compulsory option --input\n", stderr())
    help()
    quit("no", status=1)
  }
  if(! file.exists(params$in.file)){
    write(paste0("ERROR: can't find file ", params$in.file, "\n"), stderr())
    help()
    quit("no", status=1)
  }
}

run <- function(params){
  
  ## specific code ...
  
}

main <- function(){
  params <- list(verbose=1,
                 in.file=NULL)
  
  params <- parseCmdLine(params)
  
  checkParams(params)
  
  if(params$verbose > 0){
    message(paste0("START ", prog.name, " ",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    message(paste0("cwd: ", getwd()))
  }
  
  system.time(run(params))
  
  if(params$verbose > 0){
    message(paste0("END ", prog.name, " ",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    ## print(object.size(x=lapply(ls(), get)), units="Kb") # return an error I don't understand
  }
}

main()
