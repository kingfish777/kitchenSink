
library(RMySQL)
library(ggplot2)
library(data.table)
library(pROC)
#library(sendmailR)
library(ROCR)

setwd("/db/datafarm")

controlInds <- c("0", "1")

aki_tisdale <- c("Capreomycin", "Cyclosporine", "Diflunisal", "Etodolac", "Fenoprofen", "Ibuprofen", "Ketoprofen", "Ketorolac", "Mefenamate", "Naproxen", "Piroxicam", "meloxicam", "oxaprozin")
ali_tisdale <- c("Lamivudine", "Stavudine", "Zalcitabine", "Zidovudine", "abacavir", "bosentan")
gib_tisdale <- c("Diflunisal", "Etodolac", "Fenoprofen", "Flurbiprofen", "Ibuprofen", "Indomethacin", "Ketoprofen", "Ketorolac", "Mefenamate", "Naproxen", "Piroxicam", "Potassium_Chloride", "Sulindac", "Tolmetin", "meloxicam", "nabumetone", "oxaprozin", "valdecoxib")
mi_tisdale <- c("Diflunisal", "Estradiol", "Estrogens,_Conjugated", "Etodolac", "Fenoprofen", "Flurbiprofen", "Indomethacin", "Ketoprofen", "Ketorolac", "Piroxicam", "Salsalate", "Sulindac", "Tolmetin", "estropipate", "nabumetone", "oxaprozin")

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

#generateROCs <- function(condition, path, eval_metric) {
  #tisdale <- c()
#genROCs(adr, discpat, evalmet)
genROCs <- function(condition, pathway, em) {
  if (pathway == 'p') { path = "PREDISPOSES-INV==" } else if (pathway == 'c') { path = "CAUSES-INV==" }
   else if (pathway == 'tc') { path = "TREATS==COEXISTS_WITH-INV" }
  if (em == "m") { eval_metric = "medCoef" } else if (em == "m0") { eval_metric = "medCoef0" }
   else if (em == "z") { eval_metric = "zScore" } else if (em == "z0") { eval_metric = "zScore0" }
  #condition <- "mi"
  if (condition == "aki") { tisdale <- aki_tisdale } else if (condition == "ali" ) { tisdale <- ali_tisdale } else if (condition == "gib") { tisdale <- gib_tisdale } else if (condition == "mi") { tisdale <- mi_tisdale }
  #path <- "TREATS==COEXISTS_WITH-INV"
  #path <- "CAUSES-INV=="
  #eval_metric <- "medCoef"
  datamodel <- c()
  index <- 1
  #condition <- "ali"
  #path <- "CAUSES-INV=="
  # for (condition in conditions) {
  for (controlInd in controlInds) {
    #   for (path in pathways) {
    drugPaths <- paste("d3jss/", condition, "/", controlInd, "/", path, sep="")
    listPaths <- list.dirs(path = drugPaths)
    #listPaths <- listPaths[2:length(listPaths)]
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
  writeVec(dat = data.table(data), fn = paste0("tisdale_", condition, "_", path, ".csv"))
  #data
  cond_data_all <- data.table(subset(data, path %in% path))
  cond_data_all <- subset(cond_data_all, !is.na(medCnt))
  cond_data_all <- subset(cond_data_all, !is.na(med))
  cond_data_all <- subset(cond_data_all, converged == 1)
  cond_data_all <- subset(cond_data_all, medCnt >= 5)
  cond_data_tis <- subset(cond_data_all, subset = med %in% tisdale)
  print("cond_data_tis")
  print(cond_data_tis)
  cond_data_neg <- subset(cond_data_all, subset = controlInd == "0")
  print("cond_data_neg")
  print(cond_data_neg)
  cond_data <- rbind(cond_data_neg, cond_data_tis)
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
  roc <- plot.roc(data = dm.tbl, ci = TRUE, predictor = as.double(dm.tbl$rankOrder), response = as.numeric(dm.tbl$controlInd), as.numeric(dm.tbl$controlInd) ~ as.double(dm.tbl$rankOrder), col = "red", main = paste("ADR: ", unique(condition), "       discovery pathway: ", path,"  "))
  roc <- plot.roc(data = dm.tbl, ci = TRUE, predictor = as.double(dm.tbl$rankOrder), response = as.numeric(dm.tbl$controlInd), as.numeric(dm.tbl$controlInd) ~ as.double(dm.tbl$rankOrder), col = "red", main = paste("ADR: ", unique(condition), "       discovery pathway: ", path,"\n eval_metric = ", eval_metric, "\n AUC =", round(roc$auc,4), " "))
  #plot(roc)
  print(round(roc$auc,4))
  print(paste0("##############SUMMARY#########################", condition, path, eval_metric))
  print(paste0("####################", round(roc$auc, 4)))
  len <- nrow(dm.tbl)
  pc <- sum(as.numeric(dm.tbl$controlInd))
  nc <- len - pc
  print(paste0("##############SUMMARY#########total pairs    : ", len, "   positive controls: ", pc, "      negative controls: ", nc))
  print(paste0("OUTPUT: tisdale_", condition, "_", path, ".csv"))
} 

genROCs_c <- compiler::cmpfun(genROCs)
args <- commandArgs(trailingOnly=TRUE)
adr <<- as.character(args[1])
discpat <<- as.character(args[2])
evalmet <<- as.character(args[3])

#generateROCs("tisdale_mi_c.csv", "medCoef")
genROCs_c(adr, discpat, evalmet)

#genROCs()
