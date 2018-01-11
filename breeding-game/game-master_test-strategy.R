## R script for the game master
## https://github.com/timflutre/atelier-prediction-genomique
## Authors: Timothée Flutre, Jacques David
## Copyright 2016-2018 INRA, Montpellier SupAgro
## License: AGPL-3+

library(rutilstimflutre)
stopifnot(compareVersion("0.156.3",
                         as.character(packageVersion("rutilstimflutre")))
          != 1)

root.dir <- "~/src/PlantSelBreedGame/data"
setup <- getBreedingGameSetup(root.dir)

breeder <- "test"
stopifnot(breeder %in% setup$breeders)

## 1. cross the two controls once
f <- paste0(setup$init.dir, "/controls.txt")
controls <- read.table(f, stringsAsFactors=FALSE)[,1]
F1.todo <- data.frame(parent1=controls[1],
                      parent2=controls[2],
                      child="F1",
                      stringsAsFactors=FALSE)
f <- paste0(setup$breeder.dirs[[breeder]], "/F1.txt")
write.table(x=F1.todo, file=f, append=FALSE, quote=FALSE, sep="\t",
            na="", row.names=FALSE, col.names=TRUE)
## Shiny app -> request plant material

## 2. make several haplodiploidizations
nb.HDs <- 500
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
## Shiny app -> request plant material

## 3. assess the breeding value of these HDs versus the controls
f <- paste0(setup$breeder.dirs[[breeder]], "/registration.txt")
write.table(x=c("F1", paste0("F1.", 1:nb.HDs)),
            file=f, append=FALSE, quote=FALSE, row.names=FALSE,
            col.names=FALSE)
## use 'game-master_eval.R'

## 4. generate genotypes and phenotypes to test 'game-master_pheno-geno.R'
f <- paste0(setup$truth.dir, "/p0.RData")
load(f)
data.todo <- data.frame(ind=c("F1", "F1", "F1.1", "F1.2"),
                        task=c("pheno", "geno", "geno", "geno"),
                        details=c(100, "hd", "ld", p0$trait3$qtn.id))
f <- paste0(setup$breeder.dirs[[breeder]], "/data_todo.txt")
write.table(x=data.todo, file=f, append=FALSE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=TRUE)
## Shiny app -> request phenotyping and genotyping
