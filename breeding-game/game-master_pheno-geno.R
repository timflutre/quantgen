## R script for the game master
## https://github.com/timflutre/atelier-prediction-genomique
## Authors: Timothée Flutre, Jacques David
## Copyright 2016-2017 INRA, Montpellier SupAgro
## License: AGPL-3+

## TO BE CHANGED FOR EACH BREEDER
breeder <- "test"
fin <- "data_todo.txt"
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
library(MASS)
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
pre.fin <- strsplit(basename(fin), "\\.")[[1]][1]
db <- dbConnect(SQLite(), dbname=setup$dbname)

message("0. load required data")
flush.console()
f <- paste0(setup$truth.dir, "/p0.RData")
load(f)
subset.snps <- list()
f <- paste0(setup$init.dir, "/snp_coords_hd.txt.gz")
subset.snps[["hd"]] <- rownames(read.table(f))
f <- paste0(setup$init.dir, "/snp_coords_ld.txt.gz")
subset.snps[["ld"]] <- rownames(read.table(f))

message("1. read the input file from the students")
flush.console()
inds.todo <- readCheckBreedDataFile(fin, max.nb.plots=p0$K,
                                    subset.snps=subset.snps)
(data.types <- countRequestedBreedTypes(inds.todo))

message("2. check that the requested individuals already exist")
flush.console()
tbl <- paste0("plant_material_", breeder)
stopifnot(tbl %in% dbListTables(db))
query <- paste0("SELECT child FROM ", tbl)
res <- dbGetQuery(conn=db, query)
stopifnot(all(inds.todo$ind %in% res$child))

message("3. load the haplotypes and convert to genotypes")
flush.console()
X <- NULL # TODO: allocate whole matrix at this stage
for(i in 1:nrow(inds.todo)){
  ind.id <- inds.todo$ind[i]
  if(ind.id %in% rownames(X))
    next
  message(paste0(i, "/", nrow(inds.todo), " ", ind.id))

  f <- paste0(setup$truth.dir, "/", breeder, "/", ind.id, "_haplos.RData")
  if(! file.exists(f))
    stop(paste0(f, " doesn't exist"))
  load(f)

  ind$genos <- segSites2allDoses(seg.sites=ind$haplos, ind.ids=ind.id,
                                 rnd.choice.ref.all=FALSE)

  if(is.null(X)){
    X <- ind$genos
  } else
    X <- rbind(X, ind$genos)
}
## TODO: remove empty rows, if any
dim(X)

message("4. handle the 'pheno' tasks for the requested individuals")
flush.console()
idx <- which(inds.todo$task == "pheno")
length(idx)
if(length(idx) > 0){
  phenos.df <- makeDfPhenos(ind.ids=inds.todo$ind[idx],
                            nb.plots=as.numeric(inds.todo$details[idx]),
                            year=year,
                            pathogen=ifelse((year - 2005) %% 3 == 0,
                                            TRUE, FALSE))
  phenos <- simulTraits12(dat=phenos.df,
                          mu=p0$mu,
                          sigma.alpha2=p0$sigma.alpha2,
                          X=X[levels(phenos.df$ind),],
                          Beta=p0$Beta,
                          h2=p0$h2)
  phenos$trait3 <- simulTrait3(dat=phenos.df,
                               X=X[levels(phenos.df$ind),],
                               qtn.id=p0$trait3$qtn.id,
                               resist.genos=p0$trait3$resist.genos,
                               prob.resist.no.qtl=p0$trait3$prob.resist.no.qtl)

  phenos.df$trait1.raw <- phenos$Y[,1]
  phenos.df$trait2 <- phenos$Y[,2]
  phenos.df$trait3 <- phenos$trait3$y
  phenos.df$trait1 <- phenos.df$trait1.raw
  tmp <- (phenos.df$pathogen & as.logical(phenos.df$trait3))
  if(any(tmp))
    phenos.df$trait1[tmp] <- (1 - p0$prop.yield.loss) * phenos.df$trait1[tmp]

  ## write the phenotypes (all inds into the same file)
  fout <- paste0(setup$breeder.dirs[[breeder]], "/", pre.fin,
                 "_phenos-", year, ".txt.gz")
  if(file.exists(fout))
    stop(paste0(fout, " already exists"))
  write.table(x=phenos.df[, -grep("raw", colnames(phenos.df))],
              file=gzfile(fout), quote=FALSE,
              sep="\t", row.names=FALSE, col.names=TRUE)
}

message("5. handle the 'geno' tasks for the requested individuals")
flush.console()
for(dty in c("ld", "hd")){
  idx <- which(inds.todo$task == "geno" & inds.todo$details == dty &
               inds.todo$ind %in% rownames(X))
  message(paste0(dty, ": ", length(idx)))
  if(length(idx) > 0){

    ## write the genotypes (all inds into the same file)
    fout <- paste0(setup$breeder.dirs[[breeder]], "/", pre.fin,
                   "_genos-", dty, ".txt.gz")
    if(file.exists(fout))
      stop(paste0(fout, " already exists"))
    write.table(x=X[inds.todo$ind[idx], subset.snps[[dty]], drop=FALSE],
                file=gzfile(fout), quote=FALSE,
                sep="\t", row.names=TRUE, col.names=TRUE)
  }
}

message("6. handle the 'snp' tasks for the requested individuals")
flush.console()
idx <- which(inds.todo$task == "geno" & ! inds.todo$details %in% c("ld","hd"))
length(idx)
if(length(idx) > 0){
  all.genos <- data.frame(ind=inds.todo$ind[idx],
                          snp=inds.todo$details[idx],
                          geno=NA,
                          stringsAsFactors=FALSE)
  for(i in idx){
    ind.id <- inds.todo$ind[i]
    snp <- inds.todo$details[i]
    all.genos$geno[all.genos$ind == ind.id &
                   all.genos$snp == snp] <- X[ind.id, snp]
  }

  ## write the genotypes (all inds into the same file)
  fout <- paste0(setup$breeder.dirs[[breeder]], "/", pre.fin,
                 "-single-snps", ".txt.gz")
  if(file.exists(fout))
    stop(paste0(fout, " already exists"))
  write.table(x=all.genos,
              file=gzfile(fout), quote=FALSE,
              sep="\t", row.names=FALSE, col.names=TRUE)
}

message("7. log")
flush.console()
for(type in names(data.types)){
  if(data.types[type] > 0){
    query <- paste0("INSERT INTO log(breeder,year,task,quantity)",
                    " VALUES ('", breeder,
                    "', '", year,
                    "', '", type, "', '",
                    data.types[type], "')")
    res <- dbGetQuery(db, query)
  }
}
