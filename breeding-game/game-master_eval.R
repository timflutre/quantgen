## R script for the game master
## https://github.com/timflutre/atelier-prediction-genomique
## Authors: Timothée Flutre, Jacques David
## Copyright 2016-2017 INRA, Montpellier SupAgro
## License: AGPL-3+

library(RSQLite)
library(rutilstimflutre)

root.dir <- "~/work2/atelier-prog-selection-2017"
setup <- getBreedingGameSetup(root.dir)

## 1. determine if the candidates are registrable (i.e. better than controls)
candidates <- list()
f <- paste0(setup$truth.dir, "/p0.RData")
load(f)
for(breeder in c("controls", setup$breeders)){
  message(breeder)
  if(breeder == "controls"){
    f <- paste0(setup$init.dir, "/temoins.txt")
  } else
    f <- paste0(setup$breeder.dirs[[breeder]], "/inscription.txt")
  if(! file.exists(f)){
    warning(paste0(breeder, " skipped because input file doesn't exist"))
  } else{
    ind.ids <- read.table(f, header=TRUE, stringsAsFactors=FALSE)[,1]
    candidates[[breeder]] <- data.frame(ind=ind.ids,
                                        trait1=NA, trait2=NA, trait3=NA,
                                        stringsAsFactors=FALSE)
    for(i in 1:nrow(candidates[[breeder]])){
      if(breeder == "controls"){
        f <- paste0(setup$truth.dir, "/",
                    candidates[[breeder]][i,"ind"], "_haplos.RData")
      } else
        f <- paste0(setup$truth.dir, "/", breeder, "/",
                    candidates[[breeder]][i,"ind"], "_haplos.RData")
      if(! file.exists(f)){
        warning(paste0("skip ", ind.ids[i], " as it doesn't exist"))
      } else{
        load(f)
        X <- segSites2allDoses(seg.sites=ind$haplos, rnd.choice.ref.all=FALSE)
        if(breeder == "controls")
          print(paste0(candidates[[breeder]][i,"ind"], " ",
                       X[, p0$trait3$qtn.id]))
        for(trait in paste0("trait", 1:2)){
          stopifnot(all(colnames(X) == names(p0[[trait]]$beta)))
          candidates[[breeder]][i, trait] <- p0[[trait]]$mu +
            (X - 1) %*% p0[[trait]]$beta
        }
        candidates[[breeder]][i, "trait3"] <-
          as.numeric(! X[, p0$trait3$qtn.id] %in% p0$trait3$resist.genos)
        candidates[[breeder]] <- candidates[[breeder]][
            order(candidates[[breeder]]$trait1,
                  candidates[[breeder]]$trait2,
                  decreasing=TRUE),]
      }
    }
  }
}

lapply(candidates, head)
lapply(candidates, nrow)

(min.trait1 <- 1.01 * mean(candidates$controls$trait1))
min.trait2 <- 14
lapply(candidates, function(x){
  best.candidates <- x$ind[x$trait1 >= min.trait1 & x$trait2 >= min.trait2]
  x[x$ind %in% best.candidates,]
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
    budgets[[breeder]] <- setNames(rep(0, length(tasks)), tasks)
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
## breeder1= breeder2= breeder3=
