#!/usr/bin/env Rscript

## Aim: does this and that
## choose between:
## Author: Timothee Flutre
## Not copyrighted -- provided to the public domain
## or:
## Copyright (C) 2011-2013 Timothee Flutre
## License: GPLv3+

rm(list=ls())
prog.name <- "myscript.R"

help <- function(){
  txt <- paste0("`", prog.name, "' does this and that.\n")
  txt <- paste0(txt, "\nUsage: ", prog.name, " [OPTIONS] ...\n")
  txt <- paste0(txt, "\nOptions:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "  -i\t\tinput\n")
  message(txt)
}

version <- function(){
  txt <- paste0(prog.name, " 1.0\n")
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

parseArgs <- function(params){
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
    else if(args[i] == "-i"){
      params$in.file <- args[i+1]
      i <- i + 1
    }
    else
      stop(paste0("unknown option ", args[i]))
  }
  
  if(params$verbose > 0){
    message("parameters:")
    print(params)
  }
  
  return(params)
}

checkParams <- function(params){
  stopifnot(! is.null(params$input),
            file.exists(params$input))
}

main <- function(){
  params <- list(verbose=1,
                 input=NULL)
  params <- parseArgs(params)
  checkParams(params)
  if(params$verbose > 0)
    message(paste0("START ", prog.name, " (", date(), ")"))
  
  ## ... specific code ...
  
  if(params$verbose > 0)
    message(paste0("END ", prog.name, " (", date(), ")"))
}

system.time(main())
print(object.size(x=lapply(ls(), get)), units="Kb")
