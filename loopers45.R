
library(data.table)
#setwd("/home/smalec/NEWDATA2")
setwd("/db/datafarm") # /Users/smalec/data")
home <- getwd()
setwd(getwd())
docs <- seq(1:2197496)
alldocs <- seq(1:2197496)
newcol <- rep(0, 2197496)
 conditions <- c("gib") # , "gib", "ali")
 # conditions <- c("gib") #, "aki") #"aki", "ali", "gib", "mi")
levels <- c("1", "0") # , "1") # , "0") #, "0")
  predicatePaths <- c("PREDISPOSES-INV==", "CAUSES-INV==", "TREATS==COEXISTS_WITH-INV") #c("TREATS==COEXISTS_WITH-INV") #PREDISPOSES-INV==", "CAUSES-INV==")
 #predicatePaths <- c("TREATS==COEXISTS_WITH-INV") 
intersectionBehavior=0 #NA

# Rscript glmShiny.R $adr $cInd $med $pred $sq

args <- NULL

#
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  sq <- as.character(args[1])
} else { sq <- 15 }

readVec <- function(xyz) {
  if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
    as.character(unlist(read.table(xyz, stringsAsFactors = F)))
  } else { c() }  
}


generatePathForEntity <- function(ade, lev, pathType, medication, entityType) {
  if (entityType == "c") {
    resultPath <- paste(ade, "/confounding", lev, "/",  medication, "/", pathType, sep="")
  }
  if (entityType == "d") {
    resultPath <- paste(ade, "/", "med", lev, sep="")
  }
  if (entityType == "a") {
    resultPath <- paste(ade, "/", "condition", sep="")
  }
  return(resultPath)
}

runPVDATA <- function(conditions, levels, squelch) {
  for (condition in conditions) {
    print(condition)
    for (level in levels) {
      print(level)
      prelimDrugPath <- generatePathForEntity(condition, level, "",  "", "d")
      print(prelimDrugPath)
      #    print(list.files(prelimDrugPath))
      #drugList <- unique(read.table("gib_drug1.txt"))
      drugListFolder <- paste(condition, "_drug", level, ".txt", sep="")
      print(drugListFolder)
      drugList <- unique(read.table(drugListFolder, stringsAsFactors = F)[,1])
      #drugList <- sort(drugList, decreasing=T)
      for (path in predicatePaths) {
        
        for (drug in drugList) {
       # print(drug)
          #drug <- "Acyclovir"
         # system(paste("say -v Samantha -r 20000 'starting new drug: ", drug, "' ", sep=""))
	  numConfs <- 13
          #system("rm numConfs.txt")
 	  for (squelchFilter in 0:squelch) {
       #try(numConfs <- readVec("numConfs.txt"))
              #  if (squelchFilter == 0) {
                print(paste("Rscript city9.R", condition, level, drug, path, squelchFilter)) 
                  system(paste("Rscript city9.R", condition, level, drug, path, squelchFilter))
                 # numConfs <- readVec("numConfs.txt")
              #  } else if (squelchFilter > 0) {
               #     print(numConfs) 
           #       if (squelchFilter < readVec("numConfs.txt")) { 
         #           system(paste("say -v Samantha -r 20000 'squelchFilter: ", squelchFilter, "' ", sep=" "))
        #            system(paste("Rscript city233.R", condition, level, drug, path, squelchFilter))
#
         #         }  
        #      }
 	  }
        }
      }
    }
  }
}

runPVDATA(conditions, levels, sq)

