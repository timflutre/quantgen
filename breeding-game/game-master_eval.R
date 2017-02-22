## R script for the game master
## https://github.com/timflutre/atelier-prediction-genomique
## Authors: Timothée Flutre, Jacques David
## Copyright 2016-2017 INRA, Montpellier SupAgro
## License: AGPL-3+

library(RSQLite)
library(rutilstimflutre)

## root.dir <- "~/work2/atelier-prog-selection-2017"
root.dir <- getwd()
setup <- getBreedingGameSetup(root.dir)

## 1. determine if the candidates are registrable (i.e. better than controls)
## by phenotyping them, along with controls, WITH NO "YEAR" EFFECT!
fins <- c("test"="temoins_phenos-2025.txt.gz",
          "selectionneur1"="inscriptions_phenos-2025.txt.gz",
          "selectionneur2"="inscriptions_phenos-2025.txt.gz",
          "selectionneur3"="inscriptions_phenos-2025.txt.gz",
          "selectionneur4"="inscriptions_phenos-2025.txt.gz")
dat <- list()
for(i in 1:length(fins)){
  breeder <- names(fins)[i]
  fin <- paste0(setup$breeder.dirs[[breeder]], "/", fins[i])
  dat[[breeder]] <- read.table(fin, header=TRUE)
  dat[[breeder]] <- cbind("breeder"=rep(breeder, nrow(dat[[breeder]])),
                          dat[[breeder]])
}
dat <- do.call(rbind, dat)
dat$year <- as.factor(dat$year)
dat$plot <- as.factor(dat$plot)
dat$yield.ha <- dat$trait1 * dat$trait2 / 100
str(dat)

(res <- cbind(tapply(dat$trait1, dat$ind, mean),
              tapply(dat$trait2, dat$ind, mean)))
res <- cbind(res, res[,1] * res[,2] / 100)

res[order(res[,1], decreasing=TRUE),]
res[order(res[,2], decreasing=TRUE),]
res[order(res[,3], decreasing=TRUE),]

par(mar=c(7, 4, 0.5, 0.5))
tmp <- boxplot(trait1 ~ ind, data=dat, xaxt="n", las=1, ylab="trait1")
text(x=1:nlevels(dat$ind), y=par("usr")[3] - 0.5, labels=levels(dat$ind),
     xpd=TRUE, srt=45, pos=1)

1.03 * mean(dat$trait1[dat$breeder == "test"])
mean(dat$trait2[dat$breeder == "test"])
mean(dat$yield.ha[dat$breeder == "test"])

## 2. sum the cost of each breeder
db <- dbConnect(SQLite(), dbname=setup$dbname)
costs <- c("allofecundation"=1/2,
           "autofecundation"=1/25,
           "haplodiploidization"=1/5,
           "pheno"=1,
           "geno-ld"=1/3,
           "geno-hd"=1,
           "geno-single-snp"=1/50)
tasks <- names(costs)

budgets <- list()
for(breeder in setup$breeders){
  message(breeder)
  tbl <- "log"
  query <- paste0("SELECT * FROM ", tbl, " WHERE breeder='", breeder, "'")
  res <- dbGetQuery(db, query)
  budgets[[breeder]] <- setNames(rep(0, length(tasks)), tasks)
  for(task in tasks)
    budgets[[breeder]][task] <- budgets[[breeder]][task] +
      sum(res$quantity[res$task == task])
}
budgets

total <- sapply(budgets, function(budget){budget[names(costs)] * costs})
(total <- rbind(total, colSums(total)))
colSums(total)
cost.PE <- 50 # in Mendels
(total.Mendels <- total * cost.PE)
colSums(total.Mendels)
