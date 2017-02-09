## R script for the game master
## https://github.com/timflutre/atelier-prediction-genomique
## Authors: Timothée Flutre, Jacques David
## Copyright 2016-2017 INRA, Montpellier SupAgro
## License: AGPL-3+

library(RSQLite)
library(rutilstimflutre)

root.dir <- "~/work2/atelier-prog-selection-2017"
setup <- getBreedingDirs(root.dir)

## 1. determine if the candidates are registrable (i.e. better than controls)
candidates <- list()
for(breeder in setup$breeders){
  message(breeder)
  candidates[[breeder]] <- data.frame(ind=rep(NA,5),
                                      trait1=NA,
                                      trait2=NA,
                                      trait3=NA,
                                      stringsAsFactors=FALSE)
  f <- paste0(setup$breeder.dirs[[breeder]], "/inscription.txt")
  if(! file.exists(f)){
    warning(paste0(breeder, " skipped because no 'inscription.txt'"))
  } else{
    ind.ids <- read.table(f, header=TRUE, stringsAsFactors=FALSE)
    if("ind" %in% colnames(ind.ids) &
       nrow(ind.ids) == nrow(candidates[[breeder]])){
      candidates[[breeder]]$ind <- ind.ids
      for(i in 1:nrow(candidates[[breeder]])){
        f <- paste0(setup$truth.dir, "/", breeder, "/",
                    candidates[[breeder]][i,"ind"], "_haplos.RData")
        if(file.exists(f)){
          load(f)
          genos <- segSites2allDoses(ind$haplos)
          for(trait in paste0("trait", 1:2))
            candidates[[breeder]][i, trait] <- p0[[trait]]$mu +
              genos %*% p0[[trait]]$beta
          candidates[[breeder]][i, "trait3"] <-
            ! (candidates[[breeder]][i, "ind"] %in%
               names(p0$trait3$inds.resist)[p0$trait3$inds.resist])
        }
      }
    }
  }
}

candidates

lapply(candidates, function(x){
  idx <- which.max(x$trait1[x$trait2 >= 14])
  c(x$trait1[idx], x$trait2[idx], x$trait3[idx])
})

## 2. sum the cost of each breeder
db <- dbConnect(SQLite(), dbname=setup$dbname)
tasks <- c("allofecundation", "autofecundation", "haplodiploidization",
           "pheno", "geno-ld", "geno-hd", "geno-single-snp")
costs <- setNames(c(1/2, 1/25, 1/5, 1, 1/3, 1, 1/50), tasks)

budgets <- list()
for(breeder in setup$breeders){
  message(breeder)
  files <- Sys.glob(paths=paste0(setup$breeder.dirs[[breeder]], "/*.txt"))
  if(length(files) > 0){
    budgets[[breeder]] <- setNames(rep(0, 7), tasks)
    for(f in files){
      if(basename(f) == "inscription.txt")
        next
      message(f)
      dat <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
      if(colnames(dat)[1] != "ind"){
        tbl <- "log"
        query <- paste0("SELECT * FROM ", tbl, " WHERE breeder='", breeder, "'")
        res <- dbGetQuery(db, query)
        for(task in tasks[1:4])
          budgets[[breeder]][task] <- sum(res$task == task)
      } else if("pheno" %in% dat$task){
        for(i in which(dat$task == "pheno"))
          budgets[[breeder]]["pheno"] <- budgets[[breeder]]["pheno"] +
            as.numeric(dat$details[i])
      } else if("geno" %in% dat$task){
        for(i in which(dat$task == "geno")){
          if(dat$details[i] == "hd"){
            budgets[[breeder]]["geno-hd"] <- budgets[[breeder]]["geno-hd"] + 1
          } else if(dat$details[i] == "ld"){
            budgets[[breeder]]["geno-ld"] <- budgets[[breeder]]["geno-ld"] + 1
          } else
            budgets[[breeder]]["geno-single-snp"] <- budgets[[breeder]]["geno-single-snp"] + 1
        }
      }
    }
  }
}

budgets

total <- sapply(budgets, function(budget){budget[names(costs)] * costs})
(total <- rbind(total, colSums(total)))
## breeder1=2323 PE breeder2=1962 PE breeder3=6781
