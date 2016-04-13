
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

generateROCs <- function(fn, eval_metric) {
  
  print(fn)
  data.tbl <- read.table(file = fn, header=TRUE, sep=";")
  data <- data.frame(data.tbl)
  head(data)
  data.tbl
  convergenceFlag <- "1"
  cond_data <- data.table(subset(data, path %in% path))
  cond_data <- subset(cond_data, !is.na(med))
  cond_data <- subset(cond_data, !is.na(medCnt))
  condition <- cond_data$adr[[1]]
  print(cond_data)
  dm.tbl <- data.table(cond_data)
  dm.df <- as.data.frame(dm.tbl)
  #eval_metric = "medCoef0"
  dm.tbl$rankOrder <- rank(dm.df[eval_metric], na.last=TRUE, ties.method="max") #/lengthIntersections  
  glm.adeMedExp <- glm(as.numeric(dm.tbl$controlInd) ~ as.double(rankOrder), dm.tbl, family = binomial("logit"))
  attach(dm.tbl)
  roc <- plot.roc(data = dm.tbl, ci = TRUE, predictor = as.double(dm.tbl$rankOrder), response = as.numeric(dm.tbl$controlInd), as.numeric(dm.tbl$controlInd) ~ as.double(dm.tbl$rankOrder), col = "red")
  roc <- plot.roc(data = dm.tbl, ci = TRUE, predictor = as.double(dm.tbl$rankOrder), response = as.numeric(dm.tbl$controlInd), as.numeric(dm.tbl$controlInd) ~ as.double(dm.tbl$rankOrder), col = "red", main = paste("\nTisdale subset :: ADR: ", unique(adr)[[1]], " using ", path[[1]],"\n eval_metric = ", eval_metric[[1]], "\n AUC =", round(roc$auc,4), " "))
  print(round(roc$auc,2))
  print(paste0("##############SUMMARY#########################", condition, path, eval_metric))
  print(paste0("####################", round(roc$auc, 4)))
  len <- nrow(dm.tbl)
  pc <- sum(as.numeric(dm.tbl$controlInd))
  nc <- len - pc
  print(paste0("##############SUMMARY#########total pairs    : ", len, "   positive controls: ", pc, "      negative controls: ", nc))

}

generateROCs_c <- compiler::cmpfun(generateROCs)

args <- commandArgs(trailingOnly=TRUE)
filename <<- as.character(args[1])
evalmet <<- as.character(args[2])

#generateROCs("tisdale_gib_c.csv", "medCoef")
generateROCs(filename, evalmet)

#generate_ROCs_c(fn, eval_metric)

#if (length(args) > 0){  generateROCs_c(adr, pred, metric) } else {  generateROCs_c(gib, c, m)   }
#generateROCs_c(adr, pred, metric)
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


