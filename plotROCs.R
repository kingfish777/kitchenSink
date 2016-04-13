
library(RMySQL)
library(ggplot2)
library(data.table)
library(pROC)
#library(sendmailR)
library(ROCR)

setwd("/db/datafarm")

controlInds <- c("0", "1")

# pathways <- c("CAUSES-INV==", "PREDISPOSES-INV==") #
readVec <- function(xyz) {
  if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
    as.character(unlist(read.table(xyz, stringsAsFactors = F)))
  } else { "0" }  
}

readVecHeader <- function(xyz) {
  if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
    as.character(unlist(read.table(xyz, stringsAsFactors = F, header = TRUE)))
  } else { "0" }  
}

writeVec <- function(dat,fn) {
  write.table(dat, file=fn, quote=FALSE, row.names=FALSE, eol="\n", col.names=TRUE, append=FALSE, sep = ";")
}

writeVec_c <- compiler::cmpfun(writeVec)
readVec_c <- compiler::cmpfun(readVec)
readVecHeader_c <- compiler::cmpfun(readVecHeader)


cleanStrings <- function(xyz) {
  xyz <- gsub(",", replacement="COMMA", x=xyz)
  xyz <- gsub(",", replacement="COMMA", x=xyz)
  #xyz <- gsub("-", replacement="HYPHEN", x=xyz)
  #xyz <- gsub("-", replacement="HYPHEN", x=xyz)
}

generateROCs <- function(condition, path, eval_metric) {
  datamodel <- c()
  index <- 1
  #condition <- "ali"
  #path <- "CAUSES-INV=="
  # for (condition in conditions) {
    for (controlInd in controlInds) {
   #   for (path in pathways) {
        drugPaths <- paste("d3jss/", condition, "/", controlInd, "/", path, sep="")
        listPaths <- list.dirs(path = drugPaths)
        listPaths <- listPaths[2:length(listPaths)]
        for (dataPointPath in listPaths) {
          dataRow <- c()
          index <- index + 1
          print(index)
          onion <- unlist(strsplit(dataPointPath, "/"))
          lo <- length(onion)
          poison <- onion[5]
          coef.fn <- paste(dataPointPath, "/medConf.txt", sep="")
          print("COEFCOEFCOEFCOEFCOEFCOEF")
          medCoef <- readVec(coef.fn)
          print(medCoef)
          if (is.na(medCoef)) { medCoef <- 0 }
          glm.fn <- paste(dataPointPath, "/glm.txt", sep="")
          glm <- unlist(readVec(coef.fn)) # as.character(readVec(glm.fn))
          ixn.fn <- paste(dataPointPath, "/dataTotXN.txt", sep="")
          ixn <- readVecHeader(ixn.fn)
          med.fn <- paste(dataPointPath, "/dataTotMed.txt", sep="")
          medCnt <- as.numeric(unlist(readVecHeader(med.fn)))
          med0.fn <- paste(dataPointPath, "/medConf0.txt", sep="")
          medCoef0 <- as.numeric(unlist(readVec(med0.fn)))
          zScore.fn <- paste(dataPointPath, "/zDataZscore.txt", sep="")
          zScore <- as.double(readVecHeader(zScore.fn))
          zScore0.fn <- paste(dataPointPath, "/zDataZscore0.txt", sep="")
          zScore0 <- as.double(readVecHeader(zScore0.fn))
          converged.fn <- paste(dataPointPath, "/converged.txt", sep="")
          converged <- readVec(converged.fn)
          baselineOR.fn <- paste(dataPointPath, "/baseline_OR.txt", sep="")
          baselineOR <- readVec(baselineOR.fn)
          baselineORs.fn <- paste(dataPointPath, "/baseline_ORs.txt", sep="")
          baselineORs <- readVec(baselineORs.fn)
          baselineChiSqStatNS.fn <- paste(dataPointPath, "/baseline_ChiSqStatNS.txt", sep="")
          baselineChiSqStatNS <- readVec(baselineChiSqStatNS.fn)
          baselineChiSqStatSmth.fn <- paste(dataPointPath, "/baseline_ChiSqStatSmth.txt", sep="")
          baselineChiSqStatSmth <- readVec(baselineChiSqStatSmth.fn) 
          if (is.na(converged)) { converged <- 0 } 
          print(paste0(condition, controlInd, path, poison, medCoef, medCoef0, converged, ixn, medCnt, zScore, zScore0, baselineOR, baselineORs, baselineChiSqStatNS, baselineChiSqStatSmth, collapse = ": "))
          dataRow <- c(condition, controlInd, path, poison, medCoef, medCoef0, converged, ixn, medCnt, zScore, zScore0, baselineOR, baselineORs, baselineChiSqStatNS, baselineChiSqStatSmth)
          datamodel <- rbind(datamodel, dataRow)
    }
  }

  data.tbl <- data.table(datamodel)
  data <- data.frame(datamodel, stringsAsFactors=FALSE, check.names = FALSE, check.rows = FALSE, row.names = NULL)
  colnames(data.tbl) <- c("adr", "controlInd", "path", "med", "medCoef", "medCoef0", "converged", "ixn", "medCnt", "zScore", "zScore0", "baselineOR", "baselineORs", "baselineChiSqStatNS", "baselineChiSqStatSmth")
  colnames(data) <- c("adr", "controlInd", "path", "med", "medCoef", "medCoef0", "converged", "ixn", "medCnt", "zScore", "zScore0", "baselineOR", "baselineORs", "baselineChiSqStatNS", "baselineChiSqStatSmth")
  writeVec(dat = data.table(data), fn = paste0(adr, "_", path, ".csv"))
  head(data)
  data.tbl
  convergenceFlag <- "1"

  data
  cond_data <- data.table(subset(data, path %in% path))
  cond_data <- subset(cond_data, !is.na(med))
  #cond_data <- subset(cond_data, subset = med %in% c("Keto"))
  cond_data <- subset(cond_data, !is.na(medCnt))
  print(cond_data)
  dm.tbl <- data.table(cond_data)
  dm.df <- as.data.frame(dm.tbl)
  #dm.tbl
  #dm.df
  dm.tbl$rankOrder <- rank(dm.df[eval_metric], na.last=TRUE, ties.method="max") #/lengthIntersections  
  #head(dm.tbl)
  #dm.tbl
  glm.adeMedExp <- glm(as.numeric(dm.tbl$controlInd) ~ as.double(rankOrder), dm.tbl, family = binomial("logit"))
  attach(dm.tbl)
  roc <- plot.roc(data = dm.tbl, ci = TRUE, predictor = as.double(dm.tbl$rankOrder), response = as.numeric(dm.tbl$controlInd), as.numeric(dm.tbl$controlInd) ~ as.double(dm.tbl$rankOrder), col = "red", main = paste("ADR: ", unique(adr), "       discovery pathway: ", path,"  "))
  roc <- plot.roc(data = dm.tbl, ci = TRUE, predictor = as.double(dm.tbl$rankOrder), response = as.numeric(dm.tbl$controlInd), as.numeric(dm.tbl$controlInd) ~ as.double(dm.tbl$rankOrder), col = "red", main = paste("ADR: ", unique(adr), "       discovery pathway: ", path,"\n eval_metric = ", eval_metric, "\n AUC =", roc$auc, " "))
  #plot(roc)
  print(round(roc$auc,2))
  print(paste0("##############SUMMARY#########################", condition, path, eval_metric))
  print(paste0("####################", round(roc$auc, 4)))
  len <- nrow(dm.tbl)
  pc <- sum(as.numeric(dm.tbl$controlInd))
  nc <- len - pc
  print(paste0("##############SUMMARY#########total pairs    : ", len, "   positive controls: ", pc, "      negative controls: ", nc))
} 


aki <- c("aki") #aki", "ali", "gib", "mi") #, "gib") # "ali", "gib", "aki", "mi") #c("gib", "aki") #"gib", "aki")
ali <- c("ali")
gib <- c("gib")
mi <- c("mi")
#pathways <- c("TREATS==COEXISTS_WITH-INV", "CAUSES-INV==", "PREDISPOSES-INV==")
c <- "CAUSES-INV=="
p <- "PREDISPOSES-INV=="
tc <- "TREATS==COEXISTS_WITH-INV"
#evals
m0 <- "medCoef0"
m <- "medCoef"
z0 <- "zScore0"
z <- "zScore" 
ro <- "OR" # OddsRatioNS
ros <- "ORs" # OddsRatioSmth
baic <- "bestAIC"
cssnc <- "ChiSqStatNS"
cssncsm <- "ChiSqStatSmth"

generateROCs_c <- compiler::cmpfun(generateROCs)

args <- commandArgs(trailingOnly=TRUE)
adr <<- as.character(args[1])
pred <<- as.character(args[2])
metric <<- as.character(args[3])

#if (length(args) > 0){  generateROCs_c(adr, pred, metric) } else {  generateROCs_c(gib, c, m)   }
generateROCs_c(adr, pred, metric)
#generateROCs_c(gib, c, m)
#generateROCs_c(gib, tc, m)

####  sudo Rscript plotROCs.R "gib" "PREDISPOSES-INV==" "medCoef"
####  sudo Rscript plotROCs.R "gib" "TREATS==COEXISTS_WITH-INV" "medCoef0"

#Rscript plotROCs.R "aki" "CAUSES-INV==" "zScore0"
#Rscript plotROCs.R "aki" "CAUSES-INV==" "zScore"
#Rscript plotROCs.R "aki" "CAUSES-INV==" "medCoef0"
#Rscript plotROCs.R "aki" "CAUSES-INV==" "medCoef"
#Rscript plotROCs.R "aki" "PREDISPOSES-INV==" "zScore0" 
#Rscript plotROCs.R "aki" "PREDISPOSES-INV==" "zScore" 
#Rscript plotROCs.R "aki" "PREDISPOSES-INV==" "medCoef0" 
#Rscript plotROCs.R "aki" "PREDISPOSES-INV==" "medCoef" 
#Rscript plotROCs.R "aki" "TREATS==COEXISTS_WITH-INV==" "zScore0"
#Rscript plotROCs.R "aki" "TREATS==COEXISTS_WITH-INV==" "zScore"
#Rscript plotROCs.R "aki" "TREATS==COEXISTS_WITH-INV==" "medCoef0"
#Rscript plotROCs.R "aki" "TREATS==COEXISTS_WITH-INV==" "medCoef"

#Rscript plotROCs.R "ali" "CAUSES-INV==" "zScore0"
#Rscript plotROCs.R "ali" "CAUSES-INV==" "zScore"
#Rscript plotROCs.R "ali" "CAUSES-INV==" "medCoef0"
#Rscript plotROCs.R "ali" "CAUSES-INV==" "medCoef"
#Rscript plotROCs.R "ali" "PREDISPOSES-INV==" "zScore0" 
#Rscript plotROCs.R "ali" "PREDISPOSES-INV==" "zScore" 
#Rscript plotROCs.R "ali" "PREDISPOSES-INV==" "medCoef0" 
#Rscript plotROCs.R "ali" "PREDISPOSES-INV==" "medCoef" 
#Rscript plotROCs.R "ali" "TREATS==COEXISTS_WITH-INV==" "zScore0"
#Rscript plotROCs.R "ali" "TREATS==COEXISTS_WITH-INV==" "zScore"
#Rscript plotROCs.R "ali" "TREATS==COEXISTS_WITH-INV==" "medCoef0"
#Rscript plotROCs.R "ali" "TREATS==COEXISTS_WITH-INV==" "medCoef"

#Rscript plotROCs.R "gib" "CAUSES-INV==" "zScore0"
#Rscript plotROCs.R "gib" "CAUSES-INV==" "zScore"
#Rscript plotROCs.R "gib" "CAUSES-INV==" "medCoef0"
#Rscript plotROCs.R "gib" "CAUSES-INV==" "medCoef"
#Rscript plotROCs.R "gib" "PREDISPOSES-INV==" "zScore0" 
#Rscript plotROCs.R "gib" "PREDISPOSES-INV==" "zScore" 
#Rscript plotROCs.R "gib" "PREDISPOSES-INV==" "medCoef0" 
#Rscript plotROCs.R "gib" "PREDISPOSES-INV==" "medCoef" 
#Rscript plotROCs.R "gib" "TREATS==COEXISTS_WITH-INV==" "zScore0"
#Rscript plotROCs.R "gib" "TREATS==COEXISTS_WITH-INV==" "zScore"
#Rscript plotROCs.R "gib" "TREATS==COEXISTS_WITH-INV==" "medCoef0"
#Rscript plotROCs.R "gib" "TREATS==COEXISTS_WITH-INV==" "medCoef"

#Rscript plotROCs.R "mi" "CAUSES-INV==" "zScore0"
#Rscript plotROCs.R "mi" "CAUSES-INV==" "zScore0"
#Rscript plotROCs.R "mi" "CAUSES-INV==" "zScore0"
#Rscript plotROCs.R "mi" "CAUSES-INV==" "zScore0"
#Rscript plotROCs.R "mi" "PREDISPOSES-INV==" "medCoef0" 
#Rscript plotROCs.R "mi" "PREDISPOSES-INV==" "medCoef0" 
#Rscript plotROCs.R "mi" "PREDISPOSES-INV==" "medCoef0" 
#Rscript plotROCs.R "mi" "PREDISPOSES-INV==" "medCoef0" 
#Rscript plotROCs.R "mi" "TREATS==COEXISTS_WITH-INV==" "medCoef0"
#Rscript plotROCs.R "mi" "TREATS==COEXISTS_WITH-INV==" "medCoef0"
#Rscript plotROCs.R "mi" "TREATS==COEXISTS_WITH-INV==" "medCoef0"
#Rscript plotROCs.R "mi" "TREATS==COEXISTS_WITH-INV==" "medCoef0"




