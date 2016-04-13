
library(data.table)
require(plyr) 
library(dplyr)
home <- getwd()
setwd(getwd())
docs <- seq(1:2197496)
alldocs <- seq(1:2197496)
newcol <- rep(0, 2197496)
levels <- c("1", "0") 
intersectionBehavior=0 #NA

args <- NULL
dataModel <- c()
columnNames <- c()
confList <- c()
finalConfList <- c()
#clf <- c()
synonymList <- c()
# createDataModel(adr, cInd, med, pred, sq) 
args <- commandArgs(trailingOnly=TRUE)
adr <<- as.character(args[1])
cInd <<- as.character(args[2])
med <<- as.character(args[3])
pred <<- as.character(args[4])
sq <- as.character(args[5])

generatePathForEntity <- function(ade, lev, pathType, medication, entityType) {
  #if (entityType == "c") {
  #  resultPath <- paste(ade, "/confounding", lev, "/",  medication, "/", pathType, sep="")
  #}
  #if (entityType == "c") {
  #  resultPath <- paste(ade, "/confounders/confounding", lev, "/", medication, "/", pathType, sep="")
  #}
  if (entityType == "c") {
    resultPath <- paste(ade, "/confounding", lev, "/", medication, "/", pathType, sep="")
  }
  if (entityType == "d") {
    resultPath <- paste(ade, "/", "med", lev, sep="")
  }
  if (entityType == "a") {
    resultPath <- paste(ade, "/", "condition", sep="")
  }
  if (entityType == "h1") {
    resultPath <- paste0(ade, "/h1")
  }
  if (entityType == "h2") {
    resultPath <- paste0(ade ,"/h2")
  }
  return(resultPath)
}


createCTable <- function(cList) {
  #initialize variables
  no.yes <- c("No", "Yes")
  ade.df <- c()
  gLevs <- function(x, y) {
    print(paste0("x: ", x, " y:", y))
    d <- gl(2, x, y, no.yes)
    return(d)
  }
  genAde.df <- function(n) {
    sequence <- c(1)
    for (i in 1:n) {
      sequence <- c(sequence, 2^i) # was *
    }
    hyp.df <- c()
    print(sequence)
    if (n > 0) {
      hyp.df <- sapply(sequence, gLevs, (2^(2+n))/2) ###### screwy stuff here ---- was /2
      print("&&&&&&&&&&&&&&&&&&&&&&&&&&hyp.df")
      print(hyp.df)
      print("&&&&&&&&&&&&&&&&&&&&&&&&&&hyp.df")
    } else if (n == 0) {
      hyp.df <- gl(2, 1, 2, no.yes)
      hyp.df <- matrix(hyp.df, 2, 1)
      #names(hyp.df) <- "medication"
      #print(hyp.df)
    }
    return(hyp.df)
  }
  #generate yes/no/yes/no//yes/yes/no/no, etc.
  ade.df <- genAde.df(length(cList))
  #print(paste("88888888888888888888888888", cList))
  return(ade.df)
}



getData <- function(path, dataCategoryInd) {
  if (dataCategoryInd == 'medOrAde') {
    data <- tryCatch({read.table(file=path, header = F, quote = '')
    }, warning = function(w) {
      
    }, error <- function(e) { data <- '1'}, finally = { } )
  } else {
    data <- tryCatch({read.table(file=path, header = F, quote = '')
    }, warning = function(w) {
      
    }, error <- function(e) { data <- '1'; closeAllConnections(); }, finally = { } )
  }
  return(unlist(data))
}

filterConditions <- function(condition, confoundingName) {
  if (condition == 'aki') {
    if ((length(grep("nephr", confoundingName))) || (length(grep("renal", confoundingName))) || (length(grep("kidney", confoundingName)))) {
      result = FALSE
      print("!!!!!!!!!!!!!!!synonym alert!!!!!!!!!!!!!!")
      print(confoundingName)
      synonymList <- c(synonymList, confoundingName)
    } else { result = TRUE }
  } else if (condition == 'gib') {
    if ((length(grep("ulcer", confoundingName))) || (length(grep("gastro", confoundingName))) || (length(grep("bleed", confoundingName))) || (length(grep("melena", confoundingName)))) {
      result = FALSE
      print("!!!!!!!!!!!!!!!synonym alert!!!!!!!!!!!!!!")
      print(confoundingName)
      synonymList <- c(synonymList, confoundingName)
    } else { result = TRUE }
  } else if (condition == 'ali') {
    if ((length(grep("hepat", confoundingName))) || (length(grep("liver", confoundingName)))) {
      result = FALSE
      print("!!!!!!!!!!!!!!!synonym alert!!!!!!!!!!!!!!")
      print(confoundingName)
      synonymList <- c(synonymList, confoundingName)
    } else { result = TRUE }
  } else if (condition == 'mi') {
    if ((length(grep("card", confoundingName))) || (length(grep("heart", confoundingName)))) {
      result = FALSE # coron-ary
      print("!!!!!!!!!!!!! SYNONYM alert!!!!!!!!!")
      print(confoundingName)
      synonymList <- c(synonymList, confoundingName)
    } else { result = TRUE }  #... and so on
  }
  return(result)
}


cleanPunct <- function(xyz) {
  xyz <- gsub(",", replacement="COMMA", x=xyz)
  xyz <- gsub(",", replacement="COMMA", x=xyz)
  xyz <- gsub("-", replacement="HYPHEN", x=xyz)
  xyz <- gsub("-", replacement="HYPHEN", x=xyz)
}

reverseClean<- function(xyz) {
  if ((!(length(xyz)) > 0) && is.character(xyz)) {
    xyz <- gsub("COMMA", replacement=",", x=xyz)
    xyz <- gsub("COMMA", replacement=",", x=xyz)
    xyz <- gsub("HYPHEN", replacement="-", x=xyz)
    xyz <- gsub("HYPHEN", replacement="-", x=xyz)
  }
}

readVec <- function(xyz) {
  print(xyz)
  print(getwd())
  xyz <- paste0(getwd(),"/", xyz)
  print(xyz)
  dat <- c()
  if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
   try(dat <- unique(as.character(unlist(read.table(xyz)))))
   try(dat <- subset(dat, subset=(dat != "x")))
  } else { dat <- c() }
  return(dat)
}

makeUnique <- function(fName) {
  if (file.exists(fName) && (file.info(fName)$size > 0)) {
    datum <- unique(as.character(unlist(read.table(fName, stringsAsFactors = F))))
    write.table(unique(datum), file=fName, quote=FALSE, row.names=FALSE, eol="\n", col.names=FALSE, append=FALSE)
  } else { c() } 
}

writeVec <- function(dat,fn) {
  write.table(unique(dat), file=fn, quote=FALSE, row.names=FALSE, eol="\n", col.names=TRUE, append=FALSE)
}

writeVecAppend <- function(dat,fn) {
  write.table(unique(dat), file=fn, quote=FALSE, row.names=FALSE, eol="\n", col.names=FALSE, append=TRUE)
  makeUnique(fn)
}

getDataModelMetaData <- function (confList, conf, ade, medData, adeData, intersectionBehavior) {
  dataModelMetaData <- c()
  for (conf in confList) {
    confConcept <- substring(conf, 1, nchar(conf) -5)
    confDataPath <- paste("confoundingData/", confConcept, ".dat", sep="")
    columnNames <- c()
    if (file.exists(confDataPath)) {
      confData <- getData(confDataPath, "conf") ###################################
      #print(confData)
      metaData <- c()
      if (filterConditions(ade, confConcept)) { 
        print(confConcept)
        if (length(intersect(medData, intersect(confData, adeData))) > 0) { intersectionALL <- length(intersect(medData, intersect(confData, adeData))) 
        } else { intersectionALL = intersectionBehavior }
        
        #intersectionMedADE <- length(intersect(medData, adeData))
        if (length(intersect(medData, adeData)) > 0) { intersectionMedADE <- length(intersect(medData, adeData)) } else { intersectionMedADE <- intersectionBehavior }
        #intersectionConfADE <- length(intersect(confData, adeData))
        if (length(intersect(confData, adeData)) > 0) { intersectionConfADE <- length(intersect(confData, adeData)) } else { intersectionConfADE <- intersectionBehavior }
        #intersectionConfMed <- length(intersect(confData, medData))
        if (length(intersect(confData, medData)) > 0) { intersectionConfMed <- length(intersect(confData, medData)) } else { intersectionConfMed <- intersectionBehavior }
        
        metaData <- c(confConcept, intersectionALL, intersectionMedADE, intersectionConfADE, intersectionConfMed)
        dataModelMetaData <- rbind(dataModelMetaData, metaData)
        print(paste(intersectionALL, intersectionMedADE, intersectionConfADE, intersectionConfMed))
        print(paste("####################################"))
        closeAllConnections()
      }
    }
  }
  return(dataModelMetaData)
}

createDataModel <- function(ade, pairRef, med, pred, squelchFilter) {
  #SETUP
 # ade <- "gib"; pairRef <- "1"; med <- "valdecoxib"; pred <- "CAUSES-INV=="; squelchFilter = 5
  print("squelchFilter")
  print(squelchFilter)
  print(paste(ade, pairRef, med, pred))
  dThreeDIRfn <- paste("d3jss/", ade, "/", pairRef, "/", pred, "/", med, sep="")
  if ((squelchFilter == 0) && (file.exists(dThreeDIRfn))) {
    system(paste("rm ", dThreeDIRfn, "/*", sep=""))
  }
  system(paste("mkdir ", dThreeDIRfn, sep=""))
  columnNames <- c()
  dThreefn <- paste(dThreeDIRfn, "/data", sep="")
  high.fn <- paste(dThreefn, "HIGH", ".txt", sep="")
  real.fn <- paste(dThreefn, "REAL", ".txt", sep="")
  cand.fn <- paste(dThreefn, "NEW", squelchFilter, ".txt", sep="")
  totXn.fn <- paste(dThreefn, "TotXN.txt", sep="")
  totMed.fn <- paste(dThreefn, "TotMed.txt", sep="")
  newCand <- ""
  squelchFilter <- as.numeric(squelchFilter)
  oldInd <- squelchFilter - 1
  oldInd <- as.character(oldInd)
  oldnew.fn <- paste(dThreefn, "NEW", oldInd, ".txt", sep="")
  oldCand <- ""
  exclude.fn <- paste(dThreefn, "EXCLUSIONS", ".txt", sep="")
  dataModel <- c()
  confList <- NULL
  gems <- NULL
  closeAllConnections()
  ########## get confounding path
  confPath <- generatePathForEntity(ade, pairRef, pred, med, "c") # h2
  print(confPath)
  confList <- list.files(pattern=".conf", include.dirs = FALSE, path=confPath)
  print(confList)
  ########## get med and ade data
  print(paste0(ade, "/med", pairRef, "/med_", med, ".txt", sep=""))
  medData <- getData(paste(ade, "/med", pairRef, "/med_", med, ".txt", sep=""), "medOrAde")
  print(head(medData))
  print(length(medData))
  adeData <- getData(paste(ade, "/condition", sep=""), 'medOrAde')
  print(head(adeData))
  print(length(adeData))
  writeVec(length(intersect(adeData, medData)), totXn.fn)
  writeVec(length(medData), totMed.fn)
  dataModelMetaData <- getDataModelMetaData(confList, conf, ade, medData, adeData, intersectionBehavior) 
  print(dataModelMetaData)
  try(colnames(dataModelMetaData) <- c("confConceptName", "medAdrConf", "medAdr", "adrConf", "medConf"))
  try(rownames(dataModelMetaData) <- NULL)
  
  dmmd <- as.data.frame(dataModelMetaData, stringsAsFactors = F)
  for (index in 2:(ncol(dmmd))) {
    dmmd[,index] <- as.numeric(dmmd[,index]) 
  }
  print(dmmd)
  topMedAdrConf <- NULL
  topMedAdrConf <- subset(dmmd, medAdrConf >= 5)
  print(topMedAdrConf)
  
  topHighClub <- NULL
  if (ncol(topMedAdrConf) > 2) {
    for (index in 2:(ncol(topMedAdrConf))) {
      topMedAdrConf[,index] <- as.numeric(topMedAdrConf[,index]) -1
    }
    str(topMedAdrConf)
    topMedAdrConf[ , order("medAdrConf", decreasing = T) ]
    highClub <- topMedAdrConf[ order(-topMedAdrConf$medAdrConf, -topMedAdrConf$medAdr), ]
    print(highClub)
    topHighClub <- highClub[1:squelchFilter, ]
    topHighClub
  }
  
  confList <- c()
  roughGems <- c()
  try(roughGems <- as.character(highClub[, 1]))
  #try(roughGems <- roughGems)
  print("BUTTERGUPPIES")
  try(print(length(roughGems)))
  try(print(roughGems))
  if (squelchFilter <= length(roughGems)) {
    if (squelchFilter != 0) { 
      roughGems <- as.character(highClub[, 1])
      confList <- as.character(highClub[, 1]) 
      #clf <- confList
      write(length(roughGems), file="numConfs.txt")
      #confList
      #confList <- unlist(sapply(confList, cleanNames))
      writeVec(confList[1:squelchFilter], high.fn)
      exclusions <- unique(readVec(exclude.fn))
      print("**********exclusions")
      print(exclusions)
      real <- readVec(real.fn)
      print("**********real")
      print(real)
      high <- readVec(high.fn)
      print("**********high")
      print(high)
      penUlt <- squelchFilter -1
      oldCand <- confList[penUlt]
      newCand <- confList[squelchFilter]
      #oldCand <- readVec(oldnew.fn)
      print("**********oldCand")
      print(oldCand)
      if (squelchFilter == 1) {
        newCand <- high # identify new candidate ... no exclusions yet, obviously
        print(newCand) 
        confList <- c(newCand) # list
        writeVec(newCand, cand.fn) # store new candidate
      } else if (squelchFilter == 2) {         
        if (length(intersect(oldCand, real)) == 1) {
          
          penUlt <- squelchFilter -1
          oldCand <- confList[penUlt]
          #system(paste("say -v Samantha -r 20000 'hubba bubba ... entering into SEXYTIME intersection: ", length(intersect(oldCand, real)), " from model in squelchFilter", squelchFilter, " '", sep=""))
          print("^6^6^6^6^6^6^6^6^6^6^6^6^6")
          print(intersect(real, oldCand))
          # then no changes in confList
          #oldCand <- confList[penUlt]
          #confList <- c(newCand, oldCand)
          newCand <- confList[squelchFilter] #setdiff(high, oldCand) # calculate new candidate and store
          confList <- union(real, newCand)
          print("confList")
          print(confList)
          writeVec(newCand, cand.fn) # store new candidate
        } else if (length(intersect(oldCand, real)) == 0) { 
          # obviously first candidate ought to be excluded
          penUlt <- squelchFilter - 1
          oldCand <- confList[penUlt]
          #system(paste("say -v Samantha -r 20000 'excluding: ", oldCand, " from model in squelchFilter", squelchFilter, "' ", sep=""))
          
          newCand <- confList[squelchFilter]
          exclusions <- c(oldCand)
          writeVec(exclusions, exclude.fn) #write first exclusion
          confList <- setdiff(high, exclusions) # calculate confList
          writeVec(confList, high.fn) # write new high
          #new
          writeVec(newCand, cand.fn) # new candidate = confList in this special case
        } 
      } else if (squelchFilter > 2) {
        if (length(intersect(oldCand, real)) == 1) { # 11111111111111111
          if (length(exclusions) == 0) {
            #A
            #system(paste("say -v Samantha -r 20000 'YAK SE MASH, I BORAT!!!! intersection greater than 2 with no exclusions is where i am' "))
            # nothing to be done with confList
            # find new and save
            
            #newCand <-setdiff(high, real)
            print("#######")
            print(high)
            newCand <- roughGems[squelchFilter]
            confList <- union(real, newCand)
            #system(paste("say -v Samantha -r 20000 '", newCand, "'"))
            writeVec(confList, high.fn) # write new high
          } else if (length(exclusions) > 0) { # 1111111111111111111111111111
            #B
            #system(paste("say -v Samantha -r 20000 'I LIKE!!! intersection greater than 2 with exclusions is where i am' "))
            
            #identify the new candidiate
            #newCand <- setdiff(setdiff(high, real), exclusions)
            newCand <- confList[squelchFilter]
            
            #system(paste("say -v Samantha -r 20000 '", newCand, "'"))
            writeVec(newCand, cand.fn)
            #identify and store the new confList
            confList <- union(real, newCand)
            writeVec(confList, high.fn)
          }
        } else if (length(intersect(oldCand, real)) == 0) {
          if (length(exclusions) == 0) {
            #C
            # oldCand is exclusion, so just store it
            penUlt <- squelchFilter - 1
            oldCand <- roughGems[penUlt]
            exclusions <- unique(c(unique(exclusions), roughGems[penUlt])) # first exclusion
            #system(paste("say -v Samantha -r 20000 'OOOH LA LA! excluding: ", oldCand, " from model in squelchFilter", squelchFilter, "' ", sep=""))
            #system(paste("say -v Samantha -r 20000 'condition C ", oldCand, " dies a cruel death like a louche douche'"))
            writeVecAppend(exclusions, exclude.fn) # append exclusions
            # find new candidate and store
            newCand <- roughGems[squelchFilter] # n00b
            #system(paste("say -v Samantha -r 20000 '", newCand, "'"))
            writeVec(newCand, cand.fn) # new vector
            # generate new confList and store new high
            confList <- unlist(union(real, newCand))
            writeVec(confList, high.fn) 
          } else if (length(exclusions) > 0) {
            #D
            # find new exclusion, merge with old exclusions, and store
            
            penUlt <- squelchFilter - 1
            oldCand <- roughGems[penUlt]
            exclusions <- unique(c(unique(exclusions), oldCand)) # find exclusions
            #system(paste("say -v Samantha -r 20000 'Youre under arrest!", oldCand, " cuz youre the best. Serge Gainsbourg! excluding: ", oldCand, " from model in squelchFilter", squelchFilter, "' ", sep=""))
            #system(paste("say -v Samantha -r 20000 'condition D ", oldCand, " dies a cruel death like a louche douche'"))
            
            writeVecAppend(exclusions, exclude.fn) # append exclusions
            # find new candidate and store
            newCand <- roughGems[squelchFilter]
            #newCand <- setdiff(setdiff(high, real), exclusions) # n00b
            writeVec(newCand, cand.fn)
            # generate new confList and store new high
            confList <- unlist(union(real, newCand))
            writeVec(confList, high.fn) 
          }
        }
      }
      
      dataModel <- c()
      ######### get conf data
      ####

      print("finalconflist")
      print(finalConfList) ########### xconflist
      print("confList")
      print(confList)
      confList <- confList[1:length(confList)]
      print("end pottybrain")
      for (confConcept in confList) {  
        print("CONFFFFONCEPPTPTPTPTP")
        print(confList)
        print(confConcept)
        confDataPath <- paste("confoundingData/", confConcept, ".dat", sep="")
        confData <- getData(confDataPath, "conf") ###################################
        dataModel <- cbind(dataModel, newcol)
        columnNames <- c(columnNames, confConcept)
        #columnNames <- as.character(unlist(sapply(columnNames, cleanPunct)))
        colnames(dataModel) <- columnNames
        dataModel[ , confConcept][confData] <- 1
        print(paste("####################################"))
        print(head(dataModel))
        closeAllConnections()
      }
    }
    gems <- c()
    if (squelchFilter == 0) { columnNames <- c() } else { gems <- confList }
    print("pottybrain2")
    #print(confList)
    if (squelchFilter == 0) { finalConfList <- c("medication", "n.tot", "n.ade") } else if (squelchFilter > 0) { finalConfList <- c(confList, "medication", "n.tot", "n.ade") }
   # columnNames <- finalConfList
    print(columnNames)
    print(finalConfList)
    #if (squelchFilter == 0) { columnNames <- c() } else { gems <- confList[1:squelchFilter] }
    #print(paste0("gems: ", confList))
    
    ########## process med data
    dataModel <- cbind(dataModel, newcol)
    columnNames <- c(columnNames, "medication") # "med" pseudo name
    print(columnNames)
    colnames(dataModel) <- columnNames
    dataModel[ , "medication"][medData] <- 1 
    ########## process ade data
    dataModel <- cbind(dataModel, newcol)
    columnNames <- c(columnNames, "ade") # "ade" pseudo name
    colnames(dataModel) <- columnNames
    dataModel[ , "ade"][adeData] <- 1 
    ########## build hyp.df UNDER CONSTRUCTION
    print("building ade.df")
    print(paste("building empty contingency table for ", ade, pairRef, med, pred))
    if ((length(confList) == 0) || is.na(confList)) { ade.df <- createCTable(c()) } else { ade.df <- createCTable(confList) } 
    writeVec(dat = confList, fn = real.fn)
    print(confList)
    print(colnames(ade.df))
    print(ade.df)
    ########## create contingency table
    print("building dfm")
    print(paste("computing contingency table for ", ade, pairRef, med, pred))
    dfm <- as.matrix(ftable(as.data.frame(dataModel)))
    ########## merging tables
    print(dfm)
    print("merging tables: ade.df and dfm -empty contingency table with frequency given y/n")
    ncolAde.df <- ncol(ade.df)
    nrowAde.df <- nrow(ade.df)
    print(ncolAde.df)
    print(nrowAde.df)
    ade.df <- cbind(ade.df, rep(0, nrowAde.df))
    ade.df <- cbind(ade.df, rep(0, nrowAde.df))
    print(colnames(ade.df))
   # colnames(ade.df) <- c(colnames(ade.df[,1:ncolAde.df]), "n.tot", "n.ade")
    rowStr <- c("")
    print(paste0("nrowAde.df: ", nrowAde.df))
    print(paste0("ncolAde.df: ", ncolAde.df))
    for (nr in 1:nrowAde.df) {
      #nr = 1
      rowStr <- c()
      rowString <- ""
      rowFreq <- c()
      for (nc in 1:ncolAde.df) {
        elem <- ade.df[nr,nc]
        print(paste0(nr, " ", nc))
        #print(elem)
        if ((elem == "No") && (nc < ncolAde.df)) { rowStr <- c(rowStr, "0_"); print(rowStr) } 
         else if ((elem == "No") && (nc == ncolAde.df)) { rowStr <- c(rowStr, "0"); print(rowStr) }  
         else if ((elem == "Yes") && (nc < ncolAde.df)) { rowStr <- c(rowStr, "1_"); print(rowStr) }  
         else if ((elem == "Yes") && (nc == ncolAde.df)) { rowStr <- c(rowStr, "1") }; print(rowStr) }
        rowString <- paste0(rowStr, collapse="")
        print(rowString)
        rowFreq <- subset(dfm, rownames(dfm) == rowString)
        ade.df[nr,ncolAde.df+1] <- rowFreq[1]
        ade.df[nr,ncolAde.df+2] <- rowFreq[2]
    }
    colnames(ade.df) <- finalConfList
    print(ade.df)
    #data.frame(ade.df)
    
    writeVec(dat <- data.table(ade.df), fn <- paste0("ade.df", squelchFilter,".txt"))
    #print(head(ade.df))
    addf <- as.data.frame(ade.df, stringsAsFactors = FALSE)
   # addf <- data.frame(ade.df)
    data.frame(addf)
    attach(data.frame(ade.df, stringsAsFactors = FALSE))
    print("n.ade / n.tot")
    n.ade <- as.numeric(addf[,ncol(addf)])
    n.tot <- as.numeric(addf[,ncol(addf)-1])
   print(n.ade)
   print(n.tot)
    #n.tot <- as.numeric(ade.df[,ncol(ade.df)-1])
   print("monkeybrains")
    diff.tot_ade <- abs(as.numeric(n.tot) - as.numeric(n.ade))
    ade.tbl <- cbind(n.ade, diff.tot_ade)
    
    #ade.tbl <- cbind(n.ade, n.tot-n.ade)
    #ade.tbl
    print(addf)   
    #addf[, "n.tot"] <- as.numeric(addf[, "n.tot"])
    #addf[, "n.ade"] <- as.numeric(addf[, "n.ade"])

    #columnNames <- colnames(addf)
    columnNames <- finalConfList #c(finalConfList[1:length(finalConfList)], "medication", "n.tot", "n.ade")
    print(paste0("columnNames: ", columnNames, " from finalConfList"))
    print("printaddf")
    print(addf)
    # XXXXXXX
    columnNames <- as.character(unlist(sapply(columnNames[1:length(columnNames)], cleanPunct)))
    colnames(addf) <- columnNames
    print("922922922922")
    print(addf)
    addf_clean <- addf
    colnames(addf_clean) <- columnNames
    colnames(addf) <- columnNames
   
    print(addf)
    print("piss pot")
    #addf_clean <- na.omit(addf)
    dThreeDIRfn <- paste("d3jss/", ade, "/", pairRef, "/", pred, "/", med, sep="") 
    attach(addf_clean)
   attach(addf)
    #print(addf_clean)
   print("lengths: ")
   #print(length(aspirin))
  # print(str(aspirin))
   #print(length(medication))
   #print(str(medication))
  print("ade.tbl: ")
   print(ade.tbl)
   print(str(ade.tbl))
    writeVec(dat = data.table(addf_clean), fn = paste0("addf_clean", squelchFilter, ".txt"))
    getChiSq <- function(adr.tbl, ade, pairRef, pred, med, squelchFilter) { 
      chisq <- chisq.test(x = adr.tbl)
      chisqstat <- chisq$statistic
      chisqpval <- chisq$p.value
      chisqdf <- chisq$parameter
      dThreeDIRfn <- paste("d3jss/", ade, "/", pairRef, "/", pred, "/", med, sep="")
      chisqstat.fn <- paste(dThreeDIRfn, "/zDataChiSqStat.txt", sep="")
      chisqstatSQ.fn <-  paste(dThreeDIRfn, "/zDataChiSqStat", squelchFilter, ".txt", sep="")
      chisqPval.fn <- paste(dThreeDIRfn, "/zDataChiSqPval.txt", sep="")
      chisqPvalSQ.fn <-  paste(dThreeDIRfn, "/zDataChiSqPval", squelchFilter, ".txt", sep="")
      chisqDF.fn <- paste(dThreeDIRfn, "/zDataChiSqDF.txt", sep="")
      chisqDFSQ.fn <-  paste(dThreeDIRfn, "/zDataChiSqDF", squelchFilter, ".txt", sep="")
      print(chisqstat.fn)
      writeVec(fn=chisqstat.fn,dat=chisqstat)
      writeVec(fn=chisqstatSQ.fn, dat=chisqstat)
      writeVec(fn=chisqPval.fn,dat=chisqpval)
      writeVec(fn=chisqPvalSQ.fn, dat=chisqpval)
      writeVec(fn=chisqDF.fn,dat=chisqdf)
      writeVec(fn=chisqDFSQ.fn, dat=chisqdf)
    }
    
    try(getChiSq(ade.tbl, ade, pairRef, pred, med, squelchFilter))
    
    # print(head(ade.tbl))
    #sum(n.ade, n.tot)
    #sum(ade.tbl)
    rownames(ade.tbl) <- NULL
    rownames(ade.df) <- NULL
    prop.ade <- abs(as.numeric(n.ade))/as.numeric(n.tot)/sum(as.numeric(n.ade)/as.numeric(n.tot))  #### funny stuff
    names(prop.ade) <- NULL
    prop.ade
    sum(prop.ade)
    glm.ade <- NULL
  print("HERE NOW NOWNOWNOW")
    print(finalConfList)
    getStartVals <- function(n) { 
      startVals <- ""; 
      startVals <- paste("c(", paste(rep(paste("0"), n), collapse=","), ")", sep="")
      return(startVals) }
    columnNames <- as.character(unlist(sapply(columnNames, cleanPunct)))
    #glm.str <- paste("glm.ade <- glm(ade.tbl ~ ", paste(columnNames[1:(length(columnNames)-2)], collapse="+"), ", na.action(na.pass(c(NA, 0))), start = ", getStartVals(length(confList)+2), ", family = binomial(link='logit'), control=glm.control(trace=TRUE))", sep="")
    print("columnNames Names Names")
    print(columnNames)
    print("#################")
    print(columnNames[1:(length(columnNames)-2)])
    
    glm.str <- paste("glm.ade <- glm(ade.tbl ~ ", paste(columnNames[1:(length(columnNames)-2)], collapse="+"), ", family = binomial(link='logit'))", sep="")
    
    tryCatch({
      system("rm glm.R")
      print(glm.str)
      write(glm.str, file="glm.R", append=FALSE, sep="")
      source("glm.R", local=TRUE, echo = TRUE, verbose = TRUE)
      print((paste("computing contingency table for ", ade, pairRef, med, pred)))
      print(paste(ade, pairRef, med, pred))
      dataNames <- names(glm.ade$coefficients[2:length(glm.ade$coefficients)])
      print(paste0("DATANAMES", dataNames))
      cleanNames <- function(n) { result <- substr(n, start = 1, stop = nchar(n)-3); return(result) }
      cleanDataNames <- sapply(dataNames, cleanNames)
      print(paste0("CLEAN DATA NAMES: ", cleanDataNames))
      print("SUMMMMMMARY!!!")
      print(summary(glm.ade))
      cat(as.character(as.matrix(summary(glm.ade))),file=paste(dThreeDIRfn, "/glmSummary", squelchFilter, ".txt", sep=""),sep="\n",append=FALSE)
      cat(as.character(as.data.frame.model.matrix(addf_clean)),file=paste(dThreeDIRfn, "/ade-tbl", squelchFilter, ".txt", sep=""),sep="\n",append=FALSE)
      cat(as.character(glm.ade$coefficients["medicationYes"]),file=paste(dThreeDIRfn, "/medConf", squelchFilter, ".txt", sep=""),sep="\n",append=FALSE)
      cat(as.character(glm.ade$coefficients["medicationYes"]),file=paste(dThreeDIRfn, "/medConf.txt", sep=""),sep="\n",append=FALSE)
      medAdrzscore <- summary(glm.ade)$coefficients[,3]["medicationYes"]
      print("coefficients")
      print(glm.ade$coefficients["medicationYes"])
      if (glm.ade$converged) { convergeInd <- 1 } else { convergeInd <- 0 }
      #NACOEF
      #  glm.ade$coefficients["medicationYes"] <- NA
      if (is.na(glm.ade$coefficients["medicationYes"])) {
        medAdrCoef <- 0
        medAdrCoef <- glm.ade$coefficients["medicationYes"]
        
        print("simple medAdrCoef")
        print(medAdrCoef)
        glm.ade.rawname <- names(medAdrCoef)
        glm.ade.rawname.index <- nchar(glm.ade.rawname)
        glm.ade.name <- substr(glm.ade.rawname, start = 1, stop = glm.ade.rawname.index -3)
        glm.ade.name 
        # adr, cInd, med, pred, sq
        dThreeDIRfn <- paste("d3jss/", adr, "/", cInd, "/", pred, "/", med, sep="")
        # add to exclusions
        # 1.) check if exclusions exists, 2.) read exclusions 3.) concat to previous exclusions or create new exclusions file
        exclude.fn <- paste(dThreeDIRfn, "/dataEXCLUSIONS", ".txt", sep="")
        excludeSQ.fn <- paste(dThreeDIRfn, "/dataEXCLUSIONS", squelchFilter, ".txt", sep="")
        print("NACOEFFEXCLUDE NA COEFF EXCLUDE")
        print(exclude.fn)
        if (!file.exists(exclude.fn)) {
          exclusions <- c(glm.ade.name)
          exclusions <- exclusions[ ! exclusions %in% "NA" ]
          writeVec(fn = exclude.fn, dat = exclusions)
          writeVec(fn = excludeSQ.fn, dat = exclusions)
        } else if (file.exists(exclude.fn)) {
          exclusions <- unique(c(readVec(exclude.fn), glm.ade.name))
          #system(paste("rm", exclude.fn))
          # exclusions <- exclusions[!is.na(exclusions)]
          exclusions <- exclusions[ ! exclusions %in% "NA" ] 
          writeVecAppend(dat = exclusions, fn = exclude.fn)
          writeVec(fn = excludeSQ.fn, dat = exclusions)
        }
        # add to na coeff list
        #1.) check if na coeff list exists, 2.) read nacoeff, 3.) concat to previous nacoeffs or create new
        nas.fn <- paste(dThreeDIRfn, "/dataNAS", ".txt", sep="")
        if (!file.exists(nas.fn)) {
          nas <- c(glm.ade.name)
          nas <- nas[ ! nas %in% "NA" ]
          writeVec(fn = nas.fn, dat = nas)
          print("NEW NACOEFf")
          print(nas)
          print("NEW NACOEFFf")
        } else if (file.exists(nas.fn)) {
          nas <- unique(c(readVec(nas.fn)), glm.ade.name)
          nas <- nas[!is.na(nas)]
          nas <- nas[ ! nas %in% "NA" ]
          print(nas) 
          writeVecAppend(dat = nas, fn = nas.fn)
        }
      } else { 
        cat(as.character(convergeInd),file=paste(dThreeDIRfn, "/converged", squelchFilter, ".txt", sep=""),sep="\n",append=FALSE)
        cat(as.character(convergeInd),file=paste(dThreeDIRfn, "/converged.txt", sep=""),sep="\n",append=FALSE)
        cat(as.character(glm.str),file=paste(dThreeDIRfn, "/glm", squelchFilter, ".txt", sep=""),sep="\n",append=FALSE)
        cat(as.character(glm.str),file=paste(dThreeDIRfn, "/glm.txt", sep=""),sep="\n",append=FALSE)
        ### include ftable in future iterations
        
        data <- glm.ade$coefficients[2:length(glm.ade$coefficients)]
        datTab <- cbind(cleanDataNames, data)
        
        zData.raw <- summary(glm.ade)$coe[,3]
        zData <- zData.raw[2:length(zData.raw)]
        zDatTab <- cbind(cleanDataNames, zData)
        
        dThreeDIRfn <- paste("d3jss/", ade, "/", pairRef, "/", pred, "/", med, sep="") 
        
        aic.fn <- paste(dThreeDIRfn, "/dataAIC", ".txt", sep="")
        aicSQ.fn <- paste(dThreeDIRfn, "/dataAIC", squelchFilter, ".txt", sep="")
        
        aic <- summary(glm.ade)$aic
        writeVec(dat=aic, fn=aic.fn)
        writeVec(dat=aic, fn=aicSQ.fn)
        
        zscore.fn <- paste(dThreeDIRfn, '/zDataZscore.txt', sep='')
        zscoreSQ.fn <- paste(dThreeDIRfn, '/zDataZscore', squelchFilter, '.txt', sep='')
        writeVec(fn=zscore.fn, dat=medAdrzscore)
        writeVec(fn=zscoreSQ.fn, dat=medAdrzscore)
        dtab.fn <- paste(dThreeDIRfn, '/zDataTABLE', '.txt', sep='')
        dtabSQ.fn <- paste(dThreeDIRfn, '/zDataTABLE', squelchFilter, '.txt', sep='')
        writeVec(fn=dtab.fn, dat =datTab)
        writeVec(fn=dtabSQ.fn, dat = datTab)
        
        zdtab.fn <- paste(dThreeDIRfn, '/zScoreDataTABLE', '.txt', sep='')
        zdtabSQ.fn <- paste(dThreeDIRfn, '/zScoreDataTABLE', squelchFilter, '.txt', sep='')
        writeVec(fn=zdtab.fn, dat =zDatTab)
        writeVec(fn=zdtabSQ.fn, dat = zDatTab)

        print(datTab)
        sepStr=paste(',', '\n', '\t', sep='')
        getPVNODES <- function(x) {
          jsonStart <- paste('pv.nodes = [', '\n', '\t', '{ "name" : "ADR" },', sep='')
          jsonEnd <- paste('\n', ' ]; ', '\n', sep='')
          getPVNODE <- function(y) {
            nodeStr <- paste(' { "name" : "', y, '" }', sep='')
            #print(gNodes)
            return(nodeStr)
          }
          gNodes <- paste(jsonStart, paste(unlist(sapply(x, getPVNODE)), collapse=sepStr, sep=''), jsonEnd, sep='')
          #print(gNodes)
        }
        pvNodes <- getPVNODES(zDatTab[,1])
        ##print(pvNodes)
        
        getPVLinks <- function(x) {
          jsonStart <- paste("pv.links = [ ", "\n", "\t", sep="")
          jsonEnd <- paste("\n", " ]; ", "\n", sep="");
          getPVLink <- function(y, z) {
            linkStr <- paste('{ "a" : ["', y, '", "', y, '"], "z" : ["ADR", "', round(as.numeric(z), 3), '"] }', sep='')
            #print(nodeStr)
            return(linkStr)
          }
          #gNodes <- paste(jsonStart, paste(unlist(sapply(x[,1], getPVNODE, x[,2])), collapse=","), jsonEnd, sep="")
          gLinks <- c()
          for (index in 1:nrow(x)) {
            gLinks <- c(gLinks, paste(unlist(sapply(x[,1][index], getPVLink, x[,2][index])), collapse=''))
          }
          #gNodes <- paste(unlist(sapply(x[,1][index], getPVNODE, x[,2][index])), collapse=",")
          gLinks <- paste(jsonStart, paste(unlist(gLinks), collapse=sepStr, sep=""), jsonEnd, sep=" ")
          return(gLinks)
        }
        pvLinks <- getPVLinks(zDatTab)
        ##print(pvLinks)
        write(file="nodes", pvNodes)
        write(file="links", pvLinks)
        system('cat top.html nodes links bottom.html > index.html')
        system("cat index.html")
        system(paste("cp index.html ", dThreefn, ".html", sep=""))
        system(paste("cp index.html ", dThreefn, squelchFilter, ".html", sep=""))
        print(paste0("cleanDataNames: ", gems))
        if (squelchFilter == 0) { } else { writeVec(gems, real.fn) }
        if (squelchFilter == 0) { } else { writeVec(roughGems, high.fn) }
        
        print(paste(ade, pairRef, med, pred))
        #############print(anova(glm.ade, test="Chisq"))
        #system("say -v Samantha -r 200 'zip'") # #Good work!'")
      }
    }, warning = function(w) {
      print(w)
    }, error <- function(e) {
      print(e)
      write.table("bad", file="bad.txt")
      #system(paste("say -v Samantha -r 20000 'Do you do the douche'"))
      #system(paste("say -v Samantha -r 20000 '", squelchFilter, " was the louche'"))
      print(paste("FAILURE FAILURE FAILURE FAILURE FAILURE FAILURE"))
      print(paste("FAILURE!!!!", ade, pairRef, med, pred, squelchFilter)) 
      print("###############");print("###############");print("###############");print("###############");print("###############");
      try(detach(name="addf_clean"))
      closeAllConnections(); #next; 
    }, finally = {
      write.table(x=1, quote=F, row.names=F, col.names=F, file="good.txt")
      system(paste("diff ", dThreefn, "HIGH.txt ", dThreefn, "REAL.txt > highreal.txt", sep=""))
    } )
    
    search()
    try(detach(name="addf_clean"))
    search()
    write.table(x=1, quote=F, row.names=F, col.names=F, file="good.txt")
    system(paste("diff ", dThreefn, "HIGH.txt ", dThreefn, "REAL.txt > highreal.txt", sep=""))
  }
}

createDataModel_c <- compiler::cmpfun(createDataModel)

print(" adr cInd med pred sql ") 
print(paste0(adr, cInd, med, pred, sq))
#createDataModel_c(adr, cInd, med, pred, sq)

#print(length(args)) #Citalopram 7
 # if (length(args) == 0){ for (i in 1:15) { createDataModel_c("gib", "1", "Citalopram", "PREDISPOSES-INV==", i) } } else { createDataModel_c(adr, cInd, med, pred, sq)  }
#if (length(args) == 0){  createDataModel_c("gib", "1", "Citalopram", "PREDISPOSES-INV==", 1) }  else { createDataModel_c(adr, cInd, med, pred, sq)  }
createDataModel_c(adr, cInd, med, pred, sq) 


#runPVDATA(2)
###createDataModel <- function(ade, pairRef, med, pred, squelchFilter)
# gib 0 Adenosine PREDISPOSES-INV==
#createDataModel("gib", "0", "Adenosine", "PREDISPOSES-INV==", 5)
#createDataModel_c("gib", "1", "Citalopram", "CAUSES-INV==", 0)
##createDataModel("gib", "1", "valdecoxib", "PREDISPOSES-INV==", 90)
#createDataModel("gib", "0", "rosiglitazone", "PREDISPOSES-INV==", 90)
#createDataModel("gib", "1", "Ibuprofen", "PREDISPOSES-INV==", 90)
#createDataModel("gib", "1", "Clindamycin", "PREDISPOSES-INV==", 4)
#gib 1 Clindamycin PREDISPOSES-INV==

###############################################################runPVDATA()
#list.files(path="aki/med1")
