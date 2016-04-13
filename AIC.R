

#setwd("/home/kingfish/datafarm")
setwd("/db/datafarm")
conditions <- c("aki", "gib", "ali", "ami") #, "ali", "gib", "ami") #"aki") # mi", "gib", "ali", "aki") #, "gib") # "ali", "gib", "aki", "mi") #c("gib", "aki") #"gib", "aki")
controlInds <- c("0", "1")
pathways <- c("TREATS==COEXISTS_WITH-INV", "CAUSES-INV==", "PREDISPOSES-INV==")

readVec <- function(xyz) {
  if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
    unique(as.character(unlist(read.table(xyz, stringsAsFactors = F))))
  } else { c(99999999) }  
}

writeVec <- function(dat,fn) {
  write.table(unique(dat), file=fn, quote=FALSE, row.names=FALSE, eol="\n", col.names=FALSE, append=FALSE)
}

copyFiles <- function(x, i, suff) {
  x.curr <- paste(x, i, suff, sep="")
  x.new <- paste("best_", x, suff, sep="")
  command <- paste("cp ", x.curr, " ", x.new, sep="")
  print(command)
  system(command)
}



#dataAIC.ls <- list.files(pattern="^dataAIC")
#print(dataAIC.ls)

numAICs = 15


# construct file paths

for (condition in conditions) {
  for (controlInd in controlInds) {
    for (path in pathways) {
      drugPaths <- paste("d3jss/", condition, "/", controlInd, "/", path, sep="")
      #print(drugPaths)
      listPaths <- list.dirs(path = drugPaths)
      #print(listPaths)
      for (dataPointPath in listPaths) {
        lowVecVal = 999999999
        lowVecIndex = 0
        baseline.aic.fn <- paste(dataPointPath, "/baseline_aic.txt", sep="")
        baseline.aic <- 999999
        baseline.aic <- readVec(baseline.aic.fn)
        writeVec(dat = baseline.aic, fn = paste(dataPointPath, "/dataAIC0.txt", sep=""))
        writeVec(dat = baseline.aic, fn = paste(dataPointPath, "/bestDataAIC.txt", sep=""))
        lowVecVal = baseline.aic

        print(paste("be happy and merry with this baseline aic: ", baseline.aic))
        for (i in 0:numAICs) {
          print(dataPointPath)
          # for (aic.fn in dataAIC.ls) { 
          # dThreeDIRfn <- paste("d3jss/", ade, "/", pairRef, "/", pred, "/", med, sep="") 
          #aic.fn <- paste(dataPointPath, "/dataAIC", i, ".txt", sep="")
          #aicBEST.fn <- paste(dataPointPath, "/bestAIC", ".txt", sep="")
          #coef.fn.curr <- paste(dataPointPath, "/medConf", i, ".txt", sep="")
          #coef.fn <- paste(dataPointPath, "/medConf.txt", sep="")
          
          #datahtml.fn.curr <- paste(dataPointPath, "/data", i, ".html", sep="")
          #datahtml.fn <- paste(dataPointPath, "/data.html", sep="")
          #zDataChiSqStat.fn.curr <- paste(dataPointPath, "/zDataChiSqStat", i, ".txt", sep="")
          #zDataChiSqStat.fn <- paste(dataPointPath, "/zDataChiSqStat.txt", sep="")
          #zDataChiSqStat.fn.curr <- paste(dataPointPath, "/zDataChiSqStat", i, ".txt", sep="")
          #zDataChiSqStat.fn <- paste(dataPointPath, "/zDataChiSqStat.txt", sep="")
          #zDataZscore.fn.curr <- paste(dataPointPath, "/zDataZscore", i, ".txt", sep="")
          #zDataZscore.fn <- paste(dataPointPath, "/zDataZscore.txt", sep="")
          
          aic.fn.curr <- paste(dataPointPath, "/dataAIC", i, ".txt", sep="")
          print(paste("SILLY GOOSE: ", aic.fn.curr))
          if (file.exists(aic.fn.curr)) {
            aic.val.best.fn <- paste(dataPointPath, "/bestDataAIC.txt", sep="")
            aic.val.best <- readVec(xyz = aic.val.best.fn)
            aic.val.best.fn <- paste(dataPointPath, "/bestDataAIC", sep="")
            zDataChiSqStat.fn <- paste(dataPointPath, "/zDataChiSqStat", sep="")
            zDataChiSqPval.fn <- paste(dataPointPath, "/zDataChiSqPval", sep="")
            zDataZscore.fn <- paste(dataPointPath, "/zDataZscore", sep="")
            medConf.fn <- paste(dataPointPath, "/medConf", sep="")
            datahtml.fn <- paste(dataPointPath, "/data", sep="")
            aic.val <- 0
            aic.val.curr <- as.numeric(readVec(aic.fn.curr))
            print(paste("found a silly goose aic.val: ",  aic.val.curr, " aic.val.best: ", aic.val.best, " lowVecVal: ", lowVecVal))
           # print(aic.val)
            bestVecIndex.fn <- paste(dataPointPath, "/bestVecIndex.txt", sep="")
            best.aic.fn <- paste(dataPointPath, "/bestDataAIC.txt", sep="")
            
            # construct bestZdata.html, bestZreal.txt, etc. fn's
            if ((!is.na(aic.val.curr)) && (!is.null(aic.val.curr)) && (aic.val.curr <= aic.val.best)) {
              print("#############################")
              print("good vec")
              lowVecVal <- aic.val; lowVecIndex = i;
              print(aic.val)
              print("new lowVecVal")
              print(lowVecVal); 
              # move best data to new files
              print(i);
              print("##############################")
              writeVec(dat = i, fn = bestVecIndex.fn) 
              #copyFiles(aic.fn, i, ".txt")
              writeVec(dat = aic.val, fn = best.aic.fn)
              copyFiles(zDataChiSqStat.fn, i, ".txt")
              copyFiles(zDataChiSqPval.fn, i, ".txt")
              copyFiles(aic.fn.curr, i, ".txt")
              copyFiles(medConf.fn, i, ".txt")
              copyFiles(zDataZscore.fn, i, ".txt") 
            } else { print("##################"); 
                     print("bad vec"); 
                     print(aic.val); 
                     print(i);
                     print("##################") } 
          }
        }
      }
    }
  }
}



