## R script for the game master
## https://github.com/timflutre/atelier-prediction-genomique
## Authors: Timothée Flutre, Jacques David
## Copyright 2016-2017 INRA, Montpellier SupAgro
## License: AGPL-3+

## TO BE CHANGED FOR EACH BREEDER
breeder <- "test"
fin <- "todo.txt"
year <- 2015 # TO BE CHANGED ALONG THE GAME

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3){
  print(args)
  stop("usage: Rscript game-master_cross.R <breeder> <file_name> <year>")
} else{
  breeder <- args[1]
  fin <- args[2]
  year <- as.numeric(args[3])
}

library(RSQLite)
library(rutilstimflutre)

## root.dir <- "~/work2/atelier-prog-selection-2017"
root.dir <- getwd()
setup <- getBreedingGameSetup(root.dir)
fin <- paste0(setup$breeder.dirs[[breeder]], "/", fin)
print(breeder)
print(fin)
print(year)
stopifnot(breeder %in% setup$breeders)
stopifnot(file.exists(fin))
db <- dbConnect(SQLite(), dbname=setup$dbname)

message("1. read the input file from the students")
flush.console()
crosses.todo <- readCheckBreedPlantFile(fin)
(cross.types <- countRequestedBreedTypes(crosses.todo))

message("2. check the presence of new individuals in the set of existing individuals")
flush.console()
parent.ids <- unique(c(crosses.todo$parent1, crosses.todo$parent2))
parent.ids <- parent.ids[! is.na(parent.ids)]
child.ids <- crosses.todo$child
tbl <- paste0("plant_material_", breeder)
stopifnot(tbl %in% dbListTables(db))
query <- paste0("SELECT child FROM ", tbl)
res <- dbGetQuery(conn=db, query)
stopifnot(all(parent.ids %in% res$child))
stopifnot(all(! child.ids %in% res$child))

message("3. load the haplotypes of all parents")
flush.console()
parents <- list(haplos=list())
for(parent.id in parent.ids){
  if("ind" %in% ls())
    rm(ind)
  f <- paste0(setup$truth.dir, "/", breeder, "/", parent.id, "_haplos.RData")
  if(! file.exists(f))
    stop(paste0(f, " doesn't exist"))
  load(f)
  if(length(parents$haplos) == 0){ # first to insert
    parents$haplos <- ind$haplos
  } else{
    for(chr.id in names(parents$haplos))
      parents$haplos[[chr.id]] <- rbind(parents$haplos[[chr.id]],
                                        ind$haplos[[chr.id]])
  }
}
stopifnot(sapply(parents$haplos, nrow) / 2 == length(parent.ids))

message("4. perform the requested crosses")
flush.console()
new.inds <- list()
loc.crossovers <- drawLocCrossovers(crosses=crosses.todo,
                                    nb.snps=sapply(parents$haplos, ncol))
new.inds$haplos <- makeCrosses(haplos=parents$haplos,
                               crosses=crosses.todo,
                               loc.crossovers=loc.crossovers, verbose=2)

message("5. save the haplotypes of the new individuals")
flush.console()
for(new.ind.id in getIndNamesFromHaplos(new.inds$haplos)){
  message(new.ind.id)
  ind <- list(haplos=getHaplosInd(new.inds$haplos, new.ind.id))
  f <- paste0(setup$truth.dir, "/", breeder, "/", new.ind.id, "_haplos.RData")
  save(ind, file=f)
}

message("6. insert the requested crosses into their table")
flush.console()
nrow(res <- dbGetQuery(db, paste0("SELECT * FROM ", tbl)))
for(i in 1:nrow(crosses.todo)){
  message(paste0(i, "/", nrow(crosses.todo)))
  query <- paste0("INSERT INTO ", tbl, " VALUES",
                  "('", crosses.todo$parent1[i],
                  "', '", crosses.todo$parent2[i],
                  "', '", crosses.todo$child[i], "')")
  res <- dbGetQuery(conn=db, query)
}
nrow(res <- dbGetQuery(db, paste0("SELECT * FROM ", tbl)))

message("7. log")
flush.console()
for(type in names(cross.types)){
  if(cross.types[type] > 0){
    query <- paste0("INSERT INTO log(breeder,year,task,quantity)",
                    " VALUES ('", breeder,
                    "', '", year,
                    "', '", type,
                    "', '", cross.types[type], "')")
    res <- dbGetQuery(db, query)
  }
}
