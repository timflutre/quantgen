#!/usr/bin/env Rscript

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

rm(list=ls())
prog.name <- "myprogram.R"
prog.version <- "1.0.0" # http://semver.org/

R.v.maj <- as.numeric(R.version$major)
R.v.min.1 <- as.numeric(strsplit(R.version$minor, "\\.")[[1]][1])
if(R.v.maj < 2 || (R.v.maj == 2 && R.v.min.1 < 15))
    stop("require R >= 2.15 (for paste0)", call.=FALSE)

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Help
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
  txt <- paste0(txt, "Report bugs to <>.")
  message(txt)
}

##' Display version and license information on stdout
##'
##' To comply with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Version
version <- function(){
  txt <- paste0(prog.name, " ", prog.version, "\n")
  txt <- paste0(txt, "\n")
  ## choose between:
  ## txt <- paste0(txt, "Not copyrighted -- provided to the public domain\n")
  ## or:
  txt <- paste0(txt, "Copyright (C) 2011-2014 Timothée Flutre.\n")
  txt <- paste0(txt, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n")
  txt <- paste0(txt, "This is free software; see the source for copying conditions.  There is NO\n")
  txt <- paste0(txt, "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Written by Timothée Flutre.")
  message(txt)
}

##' Parse the command-line arguments
##'
##' Allow short and long options
##' @title Command-line
##' @param params list of parameters initialized with default values
##' @return List of parameters
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

##' Check the values of the command-line parameters
##'
##' @param params list of parameters
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
    start.time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(paste0("START ", prog.name, " ", start.time))
    message(paste0("cwd: ", getwd()))
  }

  system.time(run(params))

  if(params$verbose > 0){
    end.time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    difft <- as.numeric(
        difftime(as.POSIXct(end.time, format="%Y-%m-%d %H:%M:%S"),
                 as.POSIXct(start.time, format="%Y-%m-%d %H:%M:%S"),
                 units="days"))
    ## difft <- 1 + 2/24 + 3/(24*60) + 3/(24*3600) # 1d 2h 3m 4s in days
    difft.d <- floor(difft)
    difft.h <- floor(((difft - difft.d) * 24) %% 24)
    difft.m <- floor(((difft - difft.d - difft.h/24) * 24*60) %% (24 * 60))
    difft.s <- floor(((difft - difft.d - difft.h/24 - difft.m/(24*60)) *
                      24*60*60) %% (24 * 60 * 60))
    run.length <- sprintf("%02i:%02i:%02i", difft.h, difft.m, difft.s)
    message(paste0("END ", prog.name, " ", end.time, " (", run.length, ")"))
    ## print(object.size(x=lapply(ls(), get)), units="Kb") # return an error I don't understand
  }
}

if(! interactive())
    main()
