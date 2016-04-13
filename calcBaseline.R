####################################################################
# (C) University of Texas 
#     Health Science Center at Houston
#     School of Biomedical Informatics
# 
#     @author: Scott Malec
#     @last_updated: 11/14/2015
#     DISLAIMER:   {insert boilerplate prose}
#     BSD license: {insert boilerplate prose}
#####################################################################
#
#     USAGE: Rscript ./METATRONnoConf.R # adjust aki vs ali vs gib vs mi
#
#####################################################################



library(data.table)
library(foreach)
library(doParallel)
registerDoParallel(cores=4)
#setwd("/home/kingfish/testdata/testdata/pv_condition0")
#setwd("/home/kingfish/data/dat") # dat")
#setwd("/Users/smalec/Documents/data")
###setwd("/home/kingfish/NEWDATA")
#setwd("/home/smalec/datafarm")
setwd("/db/datafarm")
#setwd("/data/scott/datafarm/datafarm")
home <- getwd()
home
#require(pROC)
#require(Deducer)
library(pROC)
library(Deducer)
setwd(getwd())
docs <- seq(1:2197496)
conditions <- c("ami") #c("aki", "ali", "gib", "ami") #, "ali", "gib", "mi") #"aki") #, "ali", "gib", "mi") #, "ali", "gib", "mi") #, "ali", "gib", "mi") #, "ali", "gib", "mi") #, "aki", "mi", "gib") # , "ali", "mi", "gib") # c("mi") #
#c("gib", "ali", "mi", "aki")
levels <- c("1", "0")
generatePathForEntity <- function(ade, lev, entityType) {
  if (entityType == "c") {
    resultPath <- paste(ade, "/", "confounding", lev, sep="")
  }
  if (entityType == "d") {
    resultPath <- paste(ade, "/", "med", lev, sep="")
  }
  if (entityType == "a") {
    resultPath <- paste(ade, "/", "condition", sep="")
  }
  return(resultPath)
}

getOddsRatioNoConf <- function(med, condition) {
  numdocs <- length(docs) 
  #med <- read.table(file = "/home/kingfish/datafarm/aki/med1/med_Ibuprofen.txt", row.names = NULL)[,1]
  #med
  #length(med)
  #condition <- read.table(file = "/home/kingfish/datafarm/aki/condition", row.names = NULL)[,1]
  #condition
  #length(condition)
  #length(intersect(condition, med))
  # odds ratio
  a <- 0
  b <- 0
  c <- 0
  d <- 0
  as <- 0
  bs <- 0
  cs <- 0
  ds <- 0
  OR  <- 0
  ORs <- 0
  a <- length(intersect(med, condition))  # for NaN, add .5 for ad-hoc smoothing
  as <- a + .5
  b <- a - length(med) 
  b <- length(setdiff(med, intersect(med, condition)))
  bs <- b + .5
  c <- a - length(condition) 
  c <- length(setdiff(condition, intersect(med, condition)))
  cs <- c + .5
  d <- numdocs - (a + b + c) 
  d <- length(setdiff(docs, union(med, union(med, union(med, condition)))))
  ds <- d  + .5
  #OR <- (a*d)/(b*c)
  OR <- (a/b)/(c/d) # no smoothing
  ORs <- (as/bs)/(cs/ds) # smoothing
  print("odds ratio")
  print(paste('(a/b)/(c/d)=(a*d)/(b*c)=(', a, '*', d, ')/(', b, '*', c, ') = ', OR, sep=''))
  print(paste('(as/bs)/(cs/ds)=(as*ds)/(bs*cs)=(', as, '*', ds, ')/(', bs, '*', cs, ') = ', ORs, sep=''))
  print(OR)
  m <- matrix(c(a, b, c, d), nrow = 2)
  yes.no <- c("No", "Yes")
  exposure <- gl(k = 1, labels = yes.no, n = 2)
  #exposure
  #a
  #b
  #c
  #d
  n.tot <- c(d, b)
  n.ade <- c(c, a)
  #print(n.tot)
  #print(n.ade)
  data.frame(exposure, n.ade, n.tot)
  ade.tbl <- cbind(n.ade, n.tot - n.ade)
  #ade.tbl
  #n.tot
  #n.ade
  #ade.tbl
  glm.ade <- glm(ade.tbl ~ exposure, family=binomial(link = "logit"))
  zscore <- summary(glm.ade)$coe[,3]["exposureYes"]
  aic <- summary(glm.ade)$aic
  zscore
  medAdrCoef <- glm.ade$coefficients["exposureYes"]
  print("simple medAdrCoef")
  print(medAdrCoef)
  chi <- chisq.test(m)
  pvalue_ns <- chi$p.value
  chistat_ns <- chi$statistic
  ms <- matrix(c(as, bs, cs, ds), nrow = 2)
  chi_smoothing <- chisq.test(ms)
  pvalue_smoothing <- chi_smoothing$p.value
  chistat_smoothing <- chi_smoothing$statistic
  names(zscore) <- "zScore"
  names(OR) <- "OddsRatioNS"
  names(ORs) <- "OddsRatioSmth"
  names(a) <- "MedAdeIntersection"
  names(b) <- "MedNoAde"
  names(c) <- "AdeNoMed"
  names(d) <- "DocsMinusUnion"
  names(pvalue_ns) <- "PvalueNS"
  names(chistat_ns) <- "ChiSqStatNS" 
  names(pvalue_smoothing) <- "PvalueSmth"
  names(chistat_smoothing) <- "ChiSqStatSmth"
  names(medAdrCoef) <- "medAdrCoef"
  names(aic) <- "aic"
  #, a, b, c, d, pvalue_ns, chistat_ns, pvalue_smoothing, chistat_smoothing) <-  " "DocsMinusUnionNS", "PvalueNS", "ChiSqStatNS", "PvalueSmth", "ChiSqStatSmth")
  return(c(OR, ORs, a, b, c, d, pvalue_ns, chistat_ns, pvalue_smoothing, chistat_smoothing, medAdrCoef, zscore, aic)) # a_intersection, b_interMinusMed, c_interMinusAdr, d_numdocsMinusAdrMedInter))
  #return (c(ORs, chistat_ns))
}

readVec <- function(xyz) {
  if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
    unique(as.character(unlist(read.table(xyz, stringsAsFactors = F))))
  } else { c(99999999) }  
}

writeVec <- function(dat,fn) {
  write.table(unique(dat), file=fn, quote=FALSE, row.names=FALSE, eol="\n", col.names=FALSE, append=FALSE)
}


getFolderName <- function(adr, ci, d, path, cType) {
  fn <- paste("d3jss/", adr, "/", ci, "/", path, "/", drug, "/baseline_", cType, ".txt", sep="")
  print(fn)
  return(fn)
}

persistDataModels <- function(path, dm) {
  print("print dm")
  print(path)
  OR.fn <- getFolderName(dm$condition, dm$level, dm$drug, path, "OR")
  writeVec(dat = dm$OddsRatioNS, fn = OR.fn) #dm$oddsRatioNS, fn = OR.fn)
  ORs.fn <- getFolderName(dm$condition, dm$level, dm$drug, path, "ORs")
  writeVec(dat = dm$OddsRatioSmth, fn = ORs.fn)
  pvalue_ns.fn <- getFolderName(dm$condition, dm$level, dm$drug, path, "PvalueNS")
  writeVec(dat = dm$PvalueNS, fn = pvalue_ns.fn)
  chistat_ns.fn <- getFolderName(dm$condition, dm$level, dm$drug, path, "ChiSqStatNS")
  writeVec(dat = dm$ChiSqStatNS, fn = chistat_ns.fn)
  pvalue_smoothing.fn <- getFolderName(dm$condition, dm$level, dm$drug, path, "PvalueSmth")
  writeVec(dat = dm$PvalueSmth, fn = pvalue_smoothing.fn)
  chistat_smoothing.fn <- getFolderName(dm$condition, dm$level, dm$drug, path, "ChiSqStatSmth")
  writeVec(dat = dm$ChiSqStatSmth, fn = chistat_smoothing.fn)
  medAdrCoef.fn <- getFolderName(dm$condition, dm$level, dm$drug, path, "medAdrCoef")
  writeVec(dat = dm$medAdrCoef, fn = medAdrCoef.fn)
  zscore.fn <- getFolderName(dm$condition, dm$level, dm$drug, path,  "zScore")
  writeVec(dat = dm$zScore, fn = zscore.fn)
  aic.fn <- getFolderName(dm$condition, dm$level, dm$drug, path, "aic")
  writeVec(dat = dm$aic, fn = aic.fn)
}

paths <- c("CAUSES-INV==", "PREDISPOSES-INV==", "TREATS==COEXISTS_WITH-INV")

###################
dataModels <- c()

getDrugFromFolder <- function(n) {
  substr(n, start = 14, stop = nchar(n)-4) 
}

for (condition in conditions) {
  # level = LEV
  for (level in levels) {
    condPath <- generatePathForEntity(condition, level, "a")
    print(condPath)
    #adeData <- unique(read.table(condPath)[,1]) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SCORE
    adeData <- readVec(condPath)
    print(paste("datasize of cond: ", condPath, " >>>", length(adeData)))
    ########
    prelimDrugPath <- generatePathForEntity(condition, level, "d")
    print(prelimDrugPath)
    #    print(list.files(prelimDrugPath))
    #drugListFolder <- paste(condition, "_drug", level, ".txt", sep="")
    drugListPaths <- list.files(include.dirs = FALSE, path=prelimDrugPath)
    print(drugListPaths)
    #drugListPaths <- unique(read.table(drugListFolder)[,1])
    for (currDrugPath in drugListPaths) {
      ########
      ########### in real life, drug_$DRUGNAME.txt
      drugPath <- paste(prelimDrugPath, "/", currDrugPath, sep="")
      print(drugPath)
      drug <- getDrugFromFolder(drugPath)
      print("DRUGDRUGDRUG")
      print(drug)
      print(drug)
      if (file.exists(drugPath)) {
        print("FOUND A MED")
        medicationData <- c()
        command <- paste("cat ", drugPath, " | wc -l 2>&1", sep="")
        results <- system(command, wait=T, intern=T)
        print(paste("line count of results: ", results))
        if (results > 2) {
          #try(medicationData <- unique(read.table(drugPath)[,1])) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SCORE
          try(medicationData <- readVec(drugPath))
          print(paste("opening: ", drugPath))
          if (length(medicationData) > 0) {
            dataModel <- getOddsRatioNoConf(medicationData, adeData)
            stringsAsFactors = FALSE
            dataModel$condition <- condition
            dataModel$level <- as.numeric(level)
            dataModel$drug <- drug
            #dataModels <- rbind(dataModels, dataModel)
            for (p in paths) {
              #print(dataModel)
              persistDataModels(p, dataModel)
            }
          } 
        } else {
          medicationData <- list(c("1"))
          dataModel <- getOddsRatioNoConf(medicationData, adeData)
          dataModel$condition <- condition
          dataModel$level <- unlist(as.numeric(level))
          dataModel$drug <- drug
          #dataModels <- rbind(dataModels, dataModel)
          for (p in paths) {
            persistDataModels(p, dataModel)
            #print(drug)
          }
        }
      }
    }
  }
}


