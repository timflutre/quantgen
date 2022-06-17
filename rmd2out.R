#!/usr/bin/env Rscript

# Aim: convert an Rmd file into html or pdf from the command-line
# Copyright (C) 2016,2018,2021 Timothée Flutre
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: https://github.com/timflutre/quantgen

rm(list=ls())
prog.name <- "rmd2out.R"
prog.version <- "0.5.1" # http://semver.org/

R.v.maj <- as.numeric(R.version$major)
R.v.min.1 <- as.numeric(strsplit(R.version$minor, "\\.")[[1]][1])
if(R.v.maj < 2 || (R.v.maj == 2 && R.v.min.1 < 15))
    stop("require R >= 2.15 (for paste0)", call.=FALSE)

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Help
help <- function(){
  txt <- paste0("`", prog.name, "' converts an Rmd file into html or pdf.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Usage: ", prog.name, " [OPTIONS] ...\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Options:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "  -i, --input\tpath to the input file (Rmd format)\n")
  txt <- paste0(txt, "  -O, --outf\toutput format (default=html/pdf)\n")
  txt <- paste0(txt, "  -o, --output\tpath to the output file\n")
  txt <- paste0(txt, "\t\tdefault=<input prefix>.<format>\n")
  txt <- paste0(txt, "  -r, --root\tworking directory in which to render the document\n")
  txt <- paste0(txt, "\t\tdefault=use the directory of the input document\n")
  txt <- paste0(txt, "  -p, --param\tparameters to be passed on as a named list, separated by one space\n")
  txt <- paste0(txt, "\t\tformat: 'name1=value name2=value etc'\n")
  txt <- paste0(txt, "\t\tvalues made only of digits will be converted as numeric\n")
  txt <- paste0(txt, "\t\texample 1: -p 'n=3'\n")
  txt <- paste0(txt, "\t\texample 2: -p 'input=~/dat.txt n=3'\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Examples:\n")
  txt <- paste0(txt, "  ", prog.name, " -i myreport.Rmd -O html\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Report bugs to <timothee.flutre@inrae.fr>.")
  write(txt, stdout())
}

##' Display version and license information on stdout
##'
##' The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012). To comply with help2man (http://www.gnu.org/s/help2man), use --no-discard-stderr.
##' @title Version
version <- function(){
  txt <- paste0(prog.name, " ", prog.version, "\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Copyright (C) 2016,2018,2021 Institut national de la recherche agronomique.\n")
  txt <- paste0(txt, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Written by Timothée Flutre [cre,aut].")
  write(txt, stdout())
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
    else if(args[i] == "-O" || args[i] == "--outf"){
      params$out.format <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "-o" || args[i] == "--output"){
      params$out.file <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "-r" || args[i] == "--root"){
      params$root.dir <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "-p" || args[i] == "--param"){
      tmp <- args[i+1]
      tmp <- strsplit(tmp, " ")[[1]]
      if(length(tmp) > 0){
        params$rmd.params <- list()
        for(j in seq_along(tmp)){
          if(tmp[j] == "") # case happens when there are several successive spaces
            next
          if(! grepl("=", tmp[j])){
            msg <- "WARNING: custom parameters should be formatted as name=value"
            write(msg, stderr())
            print(tmp[j])
            next
          }
          var.name <- strsplit(tmp[j], "=")[[1]][1]
          var.value <- strsplit(tmp[j], "=")[[1]][2]
          if(grepl("^[[:digit:]]+$", var.value))
             var.value <- as.numeric(var.value)
          params$rmd.params[[var.name]] <- var.value
        }
      }
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
  if(! params$out.format %in% c("html", "pdf")){
    write(paste0("ERROR: unrecognized output format ", params$out.format, "\n"),
          stderr())
    help()
    quit("no", status=1)
  }
  library(rmarkdown)
}

params <- list(verbose=1,
               in.file=NULL,
               out.format="html",
               out.file=NULL,
               root.dir=NULL,
               rmd.params=NULL)

params <- parseCmdLine(params)

checkParams(params)

verbose <- params$verbose

if(verbose > 0){
  start.time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  write(paste0("START ", prog.name, " ", prog.version, " ", start.time),
        stdout())
  args <- commandArgs(trailingOnly=TRUE)
  write(paste("cmd-line:", prog.name, paste(args, collapse=" ")), stdout())
  write(paste0("cwd: ", getwd()), stdout())
}
if(verbose > 1){
  print(params$rmd.params)
}

system.time(
    rmarkdown::render(input=params$in.file,
                      output_format=paste0(params$out.format, "_document"),
                      output_file=params$out.file,
                      intermediates_dir=tempfile(pattern="dir"),
                      knit_root_dir=params$root.dir,
                      clean=TRUE,
                      params=params$rmd.params,
                      quiet=TRUE)
)

if(verbose > 0){
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
  write(paste0("END ", prog.name, " ", prog.version, " ", end.time,
               " (", run.length, ")"), stdout())
  ## print(object.size(x=lapply(ls(), get)), units="Kb") # return an error I don't understand
}
