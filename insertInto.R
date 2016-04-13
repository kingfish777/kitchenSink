

#library("pROC")
library("RMySQL")
#library("Deducer")
#setwd("/Users/smalec/data")

setwd("/db/datafarm")
conditions <- c("aki") #, "gib") # "ali", "gib", "aki", "mi") #c("gib", "aki") #"gib", "aki")
controlInds <- c("0", "1")
pathways <- c("TREATS==COEXISTS_WITH-INV", "CAUSES-INV==", "PREDISPOSES-INV==")
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

#readVecHeader("/db/datafarm/d3jss/aki/1/PREDISPOSES-INV==/moexipril/zDataZscore0.txt")

cleanStrings <- function(xyz) {
  xyz <- gsub(",", replacement="COMMA", x=xyz)
  xyz <- gsub(",", replacement="COMMA", x=xyz)
  #xyz <- gsub("-", replacement="HYPHEN", x=xyz)
  #xyz <- gsub("-", replacement="HYPHEN", x=xyz)
}

persistStats <- function(adr, ci, way, med, coef, glm, convergeInd, ixn, medCnt, zScore) {
  closeAllConnections()
  allCons <- dbListConnections(MySQL())
  for (con in allCons) {
    dbDisconnect(con)
  }
  glm <- cleanStrings(glm)
  m<-dbDriver("MySQL")
  con<-dbConnect(m, user="root", password='', host='139.52.155.151', dbname='shinypv')
  query = paste('INSERT INTO shinypv.ftabpvconfMARCH (adr, ci, way, med, coef, glm, converged, ixn, meds, zScore) VALUES ("', adr, '", "', ci, '", "', way, '", "', med, '", "', coef, '", "', glm, '", "', converged, '", "', ixn, '", "', medCnt, '", "', zScore, '");', sep='')
  print(query)
  res <- dbSendQuery(con, query)
  dbClearResult(res)
  dbDisconnect(con)
}


index <- 1
for (condition in conditions) {
  for (controlInd in controlInds) {
    for (path in pathways) {
      drugPaths <- paste("d3jss/", condition, "/", controlInd, "/", path, sep="")
      #print(drugPaths)
      listPaths <- list.dirs(path = drugPaths)
      #print(listPaths)
      for (dataPointPath in listPaths) {
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
        medCnt <- as.numeric(unlist(readVec(med.fn)))
        zScore.fn <- paste(dataPointPath, "/zDataZscore.txt", sep="")
        zScore <- as.numeric(readVecHeader(zScore.fn))
        #print(glm)
        #glm <- '""'
        converged.fn <- paste(dataPointPath, "/converged.txt", sep="")
        converged <- readVec(converged.fn)
        if (is.na(converged)) { converged <- 0 } 
        #print(paste0("condition controlInd path poison medCoef converged ixn medCnt zScore"))
        print(paste0(condition, controlInd, path, poison, medCoef, converged, ixn, medCnt, zScore, collapse = ": "))
        try(persistStats(condition, controlInd, path, poison, medCoef, glm, converged, ixn, medCnt, zScore))
        #print("GAINSBARRE")
        #print(poison)
      }
    }
  }
}
