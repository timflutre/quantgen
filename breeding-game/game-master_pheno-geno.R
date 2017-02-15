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

root.dir <- "~/work2/atelier-prog-selection-2017"
setup <- getBreedingGameSetup(root.dir)
fin <- paste0(setup$breeder.dirs[[breeder]], "/", fin)
print(breeder)
print(fin)
print(year)
stopifnot(breeder %in% setup$breeders)
stopifnot(file.exists(fin))
pre.fin <- strsplit(basename(fin), "\\.")[[1]][1]

## 0. load required data
f <- paste0(setup$truth.dir, "/p0.RData")
load(f)
subset.snps <- list()
f <- paste0(setup$init.dir, "/snp_coords_hd.txt.gz")
subset.snps[["hd"]] <- rownames(read.table(f))
f <- paste0(setup$init.dir, "/snp_coords_ld.txt.gz")
subset.snps[["ld"]] <- rownames(read.table(f))

## 1. read the input file from the students
inds.todo <- readCheckBreedDataFile(fin, max.nb.plots=p0$K,
                                    subset.snps=subset.snps)
(data.types <- countRequestedBreedTypes(inds.todo))

## 2. check that the requested individuals already exist
db <- dbConnect(SQLite(), dbname=setup$dbname)
tbl <- paste0("crosses_", breeder)
stopifnot(tbl %in% dbListTables(db))
query <- paste0("SELECT child FROM ", tbl)
res <- dbGetQuery(conn=db, query)
stopifnot(all(inds.todo$ind %in% res$child))

## 3. load the haplotypes and convert to genotypes
X <- NULL
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
dim(X)

## 4. handle the "pheno" tasks for the requested individuals
idx <- which(inds.todo$task == "pheno")
length(idx)
if(length(idx) > 0){
  nb.plots <- as.numeric(inds.todo$details[idx])
  phenos.df <- data.frame(ind=rep(inds.todo$ind[idx], nb.plots),
                          year=as.factor(rep(year, sum(nb.plots))),
                          plot=as.character(1:sum(nb.plots)),
                          pathogen=ifelse((year - 2015) %% 3 == 0, TRUE, FALSE),
                          trait1.raw=NA,
                          trait1=NA,
                          trait2=NA,
                          trait3=NA)

  ## draw phenotypes
  phenos <- list(N=nrow(phenos.df),
                 I=nlevels(phenos.df$ind),
                 J=nlevels(phenos.df$year),
                 K=nlevels(phenos.df$plot),
                 P=ncol(X))
  if(phenos$I == 1){
    phenos$Z.I <- matrix(1, phenos$N)
  } else
    phenos$Z.I <- model.matrix(~ phenos.df$ind - 1)
  phenos$Z.J <- matrix(1, nrow=phenos$N)
  phenos$Z.K <- model.matrix(~ phenos.df$plot - 1)

  phenos$alpha <- setNames(rnorm(n=phenos$J, mean=0, sd=sqrt(p0$sigma.alpha2)),
                           levels(phenos.df$year))
  phenos$gamma <- setNames(rnorm(n=phenos$K, mean=0, sd=sqrt(p0$sigma.gamma2)),
                           levels(phenos.df$plot))

  ## as in simulTraits12Rnd()
  X.tmp <- X[inds.todo$ind[idx], , drop=FALSE]
  G.A <- (X.tmp - 1) %*% cbind(p0$trait1$beta[colnames(X.tmp)],
                               p0$trait2$beta[colnames(X.tmp)])
  Sigma <- matrix(c(p0$trait1$sigma2, 0, 0, p0$trait2$sigma2), nrow=2, ncol=2)
  Epsilon <- mvrnorm(n=nrow(X.tmp), mu=c(0,0), Sigma=Sigma)
  Y <- G.A + Epsilon
  for(t in 1:2){
    trait <- paste0("trait", t)
    phenos[[trait]]$y <- p0[[trait]]$mu + phenos$Z.J %*% phenos$alpha +
      phenos$Z.K %*% phenos$gamma + phenos$Z.I %*% Y[,t]
  }

  phenos$trait3 <- list(y=rep(0, phenos$N)) # "0" means "no symptoms"
  phenos$trait3$y[phenos.df$pathogen] <- 1
  idx <- which(phenos.df$pathogen & (X.tmp[, p0$trait3$qtn.id] %in%
                                     p0$trait3$resist.genos))
  if(length(idx) > 0)
    phenos$trait3$y[idx] <- 1 - rbinom(n=length(idx), size=1,
                                       prob=p0$trait3$h2)

  phenos.df$trait1.raw <- phenos$trait1$y[,1]
  phenos.df$trait2 <- phenos$trait2$y[,1]
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

## 5. handle the "geno" tasks for the requested individuals
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

## 6. handle the "snp" tasks for the requested individuals
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

## 7. log
for(type in names(data.types)){
  if(data.types[type] > 0){
    query <- paste0("INSERT INTO log(breeder,task,quantity)",
                    " VALUES ('", breeder, "', '", type, "', '",
                    data.types[type], "')")
    res <- dbGetQuery(db, query)
  }
}
