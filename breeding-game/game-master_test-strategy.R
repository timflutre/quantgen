## R script for the game master
## https://github.com/timflutre/atelier-prediction-genomique
## Authors: Timothée Flutre, Jacques David
## Copyright 2016-2017 INRA, Montpellier SupAgro
## License: AGPL-3+

library(rutilstimflutre)

root.dir <- "~/work2/atelier-prog-selection-2017"
setup <- getBreedingGameSetup(root.dir)

breeder <- "test"
stopifnot(breeder %in% setup$breeders)

## 1. cross the two controls once
f <- paste0(setup$init.dir, "/temoins.txt")
controls <- read.table(f, stringsAsFactors=FALSE)[,1]
F1.todo <- data.frame(parent1=controls[1],
                      parent2=controls[2],
                      child="F1",
                      stringsAsFactors=FALSE)
f <- paste0(setup$breeder.dirs[[breeder]], "/F1.txt")
write.table(x=F1.todo, file=f, append=FALSE, quote=FALSE, sep="\t",
            na="", row.names=FALSE, col.names=TRUE)
## use 'game-master_cross.R'

## 2. make several haplodiploidizations
nb.HDs <- 200
HD.todo <- data.frame(parent1=NA,
                      parent2=NA,
                      child=NA,
                      stringsAsFactors=FALSE)
for(i in 1:nb.HDs)
  HD.todo <- rbind(HD.todo,
                   c("F1", NA, paste0("F1.", i)))
HD.todo <- HD.todo[-1,]
str(HD.todo)
f <- paste0(setup$breeder.dirs[[breeder]], "/HD.txt")
write.table(x=HD.todo, file=f, append=FALSE, quote=FALSE, sep="\t",
            na="", row.names=FALSE, col.names=TRUE)
## use 'game-master_cross.R'

## 3. assess the breeding value of these HDs versus their parents
f <- paste0(setup$truth.dir, "/p0.RData")
load(f)
candidates <- list(test=data.frame(ind=c(controls[1], controls[2],
                                         paste0("F1.", 1:nb.HDs)),
                                   trait1=NA, trait2=NA, trait3=NA,
                                   stringsAsFactors=FALSE))
for(i in 1:nrow(candidates[[breeder]])){
  ind.id <- candidates[[breeder]]$ind[i]
  message(ind.id)
  f <- paste0(setup$truth.dir, "/", breeder, "/",
              candidates[[breeder]][i,"ind"], "_haplos.RData")
  load(f)
  genos <- segSites2allDoses(ind$haplos)
  for(trait in paste0("trait", 1:2)){
    stopifnot(all(colnames(genos) == names(p0[[trait]]$beta)))
    candidates[[breeder]][i, trait] <- p0[[trait]]$mu +
      (genos - 1) %*% p0[[trait]]$beta
  }
  candidates[[breeder]][i, "trait3"] <- genos[, p0$trait3$qtn.id] %in%
    p0$trait3$resistant.genos
}
candidates$test <- candidates$test[order(candidates$test$trait1,
                                         decreasing=TRUE),]
table(candidates$test$trait3)
nrow(candidates$test[candidates$test$trait2 >= 14,])
rbind(candidates$test[candidates$test$ind %in% c(controls[1], controls[2]),],
      candidates$test[candidates$test$trait2 >= 14,][1:3,])
## the best offspring has trait1=41.47710 and trait2=14.35004
## whereas the best parent has trait1=39.93625 and trait2=14.27491

## 4. generate genotypes and phenotypes to test 'game-master_pheno-geno.R'
data.todo <- data.frame(ind=c("F1", "F1", "F1.1", "F1.2"),
                        task=c("pheno", "geno", "geno", "geno"),
                        details=c(100, "hd", "ld", p0$trait3$qtn.id))
f <- paste0(setup$breeder.dirs[[breeder]], "/data_todo.txt")
write.table(x=data.todo, file=f, append=FALSE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=TRUE)
## use 'game-master_pheno-geno.R'
