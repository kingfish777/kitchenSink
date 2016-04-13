####################################################################
# (C) University of Texas
# Health Science Center at Houston
# School of Biomedical Informatics
#
# @author: Scott Malec
# @last_updated: 06/04/2015
# DISLAIMER: {insert boilerplate prose}
# BSD license: {insert boilerplate prose}
#####################################################################
#
# USAGE: Rscript RxNormMedTermExpansion.R
# assuming setup of the UMLS database in MySQL and a list
# of input files with drug names
#
# USE CASE: INPUT: Ibuprofen ("trivial name" in biomed literature),
# DESIRED OUPUT: Advil, Motrin (tradenames)
#
# in /dataRaw folder, aki1, aki0, ali1, ali0 ...
# in aki1 file, list of 'trivial' names of meds
#
# requires that you will have a local installation of SemMedDB
#
#  ... works best if executed locally from root
#####################################################################
library("XML")
library("tm")
library("RCurl")
library("RJSONIO")
library("RMySQL")
library("rjson")
library("stringr")
is_simple_error <- function(x) inherits(x, "simpleError")
#homeFolderName <- "/home/smalec/pv/pv_conditions"
homeFolderName <- "/db/mayo"
getCUI <- function(drugName) {
  #drugName <- "Ibuprofen"
  m<-dbDriver("MySQL")
  con<-dbConnect(m, user='root', password='', host='139.52.155.223', dbname='umls');
  query = paste("select c.CUI from umls.CONCEPT AS c where UPPER(c.PREFERRED_NAME) = UPPER('", drugName, "') AND c.TYPE='META';", sep="")
  print(query)
  res <- dbSendQuery(con, query)
  data <- dbFetch(res)
  print(data)
  #print(i)
  #print(unlist(data[["TVTERM"]]))
  #print(sqlTab$condition[[i]])
  #newdf <- data.frame(adeNameGen<-sqlTab$condition[[i]], drugGen<-drug_orig, adeNameSpec<-basicAdrNames[[i]], adeCUISpec<-sqlTab$basicAdrCUIs[[i]], drugNameSpec<-unlist(data[["PREFERRED_NAME"]]), drugCUISpec<-unlist(data[["CUI"]], drugtvterm<-unlist(data[["TVTERM"]])))
  #pvDF <- rbind(pvDF, newdf)
  #fName <- paste("/home/smalec/pv/PV_conditions/dataCUI/", as.character(pairNames[drugIndex]), "/", drug, ".drug", sep="")
  #write.table(data, file=fName, quote = FALSE, sep = ' ', row.names = FALSE, col.names = FALSE)
  dbClearResult(res)
  dbDisconnect(con)
  return(data)
}
getLukeyCUI <- function(drugCUI) {
  drugName <- "Ibuprofen"
  drugCUI <- "C0020740"
  m<-dbDriver("MySQL")
  con<-dbConnect(m, user='root', password='', host='139.52.155.223', dbname='umls');
  # do some more stuff here
  #query = paste("select CONCAT(c.CUI, c.PREFERRED_NAME) AS lukeyCUI from umls.CONCEPT AS c where UPPER(c.PREFERRED_NAME) = UPPER('", drugName, "') AND c.TYPE='META';", sep="")
  query = paste("select LOWER(CONCAT('umls:', c.CUI, '_', c.PREFERRED_NAME)) AS lukeyCUI from umls.CONCEPT AS c where UPPER(c.CUI) = UPPER('", drugCUI, "') AND c.TYPE='META';", sep="")
  #print(query)
  res <- dbSendQuery(con, query)
  data <- dbFetch(res)
  print(data)
  #data <- "pref, name, 2"
  indexSpace <- regexpr(" ", data)[1]
  indexSpace
  if (indexSpace > 0) {
    data <- substring(data, 0, indexSpace-1)
    print(data)
  }
  indexComma <- regexpr(",", data)[1]
  indexComma
  if (indexComma > 0) {
    data <- substring(data, 0, indexComma-1)
    data
  }
  indexHyphen <- regexpr("-", data)[1]
  indexHyphen
  if (indexHyphen > 0) {
    data <- substring(data, 0, indexHyphen-1)
    data
  }
  indexForwardSlash <- regexpr("/", data)[1]
  indexForwardSlash
  if (indexForwardSlash > 0) {
    data <- substring(data, 0, indexForwardSlash-1)
    data
  }
  dbClearResult(res)
  dbDisconnect(con)
  return(data)
}
#getLukeyCUI("C0020740")
#drugs <- c( "/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_1.txt", "/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_0.txt",
#"/home/kingfish/pv/pv_conditions/dataRaw/drugs_ali_1.txt", "/home/kingfish/pv/pv_conditions/dataRaw/drugs_ali_0.txt",
#"/home/kingfish/pv/pv_conditions/dataRaw/drugs_mi_1.txt", "/home/kingfish/pv/pv_conditions/dataRaw/drugs_mi_0.txt",
#"/home/kingfish/pv/pv_conditions/dataRaw/drugs_gib_0.txt", "/home/kingfish/pv/pv_conditions/dataRaw/drugs_gib_1.txt")
setwd("/home/smalec/pv/pv_conditions/dataRaw/")
drugs <- c( "aki1", "aki0", "ali1", "ali0", "gib1", "gib0", "mi1", "mi0" )
#pairNames <- c("aki0", "aki1")
# drugs <- c( "/home/smalec/pv/pv_conditions/dataRaw/aki1", "/home/smalec/pv/pv_conditions/dataRaw/aki0",
# "/home/kingfish/pv/pv_conditions/dataRaw/ali1", "/home/kingfish/pv/pv_conditions/dataRaw/ali0",
# "/home/kingfish/pv/pv_conditions/dataRaw/mi1", "/home/kingfish/pv/pv_conditions/dataRaw/mi0",
# "/home/kingfish/pv/pv_conditions/dataRaw/gib1", "/home/kingfish/pv/pv_conditions/dataRaw/gib0")
#drugs <- c( "/home/smalec/pv/pv_conditions/dataRaw/Ibuprofen" )
#################
getExpandedForms <- function (eidUrl, pn, drug) {
  medVariants <- c(stringsAsFactors=FALSE)
  json_file <- paste(homeFolderName, "json_file.txt", sep="")
  download.file(eidUrl, json_file)
  #document <- fromJSON(file = json_file)
  #eid <- fromJSON(getURL(eidUrl))
  eid <- fromJSON(file = json_file)
  print("#########obtained json object")
  print(eid$idGroup["rxnormId"])
  options(show.error.messages = FALSE)
  if (!is.null(eid$idGroup["rxnormId"]))
    try({entityRxNormId <- eid$idGroup["rxnormId"]
    options(show.error.messages = FALSE)
    print(entityRxNormId)
    print("########")
    ##### http://rxnav.nlm.nih.gov/REST/rxcui/174742/related?rela=tradename_of+has_precise_ingredient
    url <- paste("http://rxnav.nlm.nih.gov/REST/rxcui/", entityRxNormId, "/related?tty=SCD+DF+DFG+GPCK+SBD+SBDF+SBDC+SBDG+BN+IN+MIN+PIN", sep="") #scd SBDC BN
    #url <- paste("http://rxnav.nlm.nih.gov/REST/rxcui/", entityRxNormId, "/allrelated/", sep="") # formerly /allrelated/ in RxNorm SOAP API
    print(url)
    print("########")
    tale <- xmlTreeParse(getURL(url), useInternal = T)
    #print(tale)
    #print(length(tale))
    #options(show.error.messages = FALSE)
    #print(xmlValue(getNodeSet(tale, "//rxnormdata//allRelatedGroup//conceptGroup//conceptProperties//umlscui")[[1]]))
    entities <- c()
    for (i in 1:100000) {
      options(show.error.messages = FALSE)
      try({
        newEntity <- xmlValue(getNodeSet(tale, "//rxnormdata//relatedGroup//conceptGroup//conceptProperties//umlscui")[[i]])
        #print(xmlValue(getNodeSet(tale, "//rxnormdata//relatedGroup//conceptGroup//conceptProperties//name")[[i]]))
        #print(newEntity)
        entities <- c(entities, newEntity)
        medVariants <- c(medVariants, entities)
        print(entities)
        print(medVariants)
        #entities
        #fName <- paste("/home/kingfish/pv/pv_conditions/dataExpanded/", pn, "/", drug, sep="")
        #stopCharNum <- nchar(fName) - 5
        #fileName <- substr(fName, 0, stopCharNum)
        #print("writing to file")
        #print(fileName)
        #write.table(entities, file=fileName, quote = FALSE, sep = ' ', row.names = FALSE, col.names = FALSE, append = TRUE)
        #close(fName)
      })
      if (class(newEntity) == "try-error") next;
    }
    })
  if (class(entityRxNormId) == "try-error") next;
  return(unique(medVariants))
}
#olmesartan_medoxomil
variants <- getExpandedForms("http://rxnav.nlm.nih.gov/REST/rxcui.json?idtype=UMLSCUI&id=C0386393","/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_1.txt","/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_1.txt")
#niacin
variants <- getExpandedForms("http://rxnav.nlm.nih.gov/REST/rxcui.json?idtype=UMLSCUI&id=C0027996","/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_1.txt","/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_1.txt")
#darunavir
variants <- getExpandedForms("http://rxnav.nlm.nih.gov/REST/rxcui.json?idtype=UMLSCUI&id=C1435444","/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_1.txt","/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_1.txt")

#sitagliptin
variants <- getExpandedForms("http://rxnav.nlm.nih.gov/REST/rxcui.json?idtype=UMLSCUI&id=C1565750","/home/smalec/pv/PV_conditions/dataRaw/drugs_ali_0.txt","/home/smalec/pv/PV_conditions/dataRaw/drugs_ali_0.txt")


#variants <- getExpandedForms(eidUrl, as.character(pn), as.character(df))

#drugList[df][dL][i] <-
print("##############DUKEY LUKEYS########################")
dope <- unlist(unique(variants))

fileName <- paste("/db/mayo/", "drug_sitagliptin", ".txt", sep="")
print(fileName)
write.table(dope, file=fileName, quote = FALSE, sep = '', row.names = FALSE, col.names = FALSE, append=FALSE)



#####################
#####################
drugList <- list()
dope <- list(stringsAsFactors=FALSE)
#drugs <- c( "/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_0.txt", "/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_1.txt") #,
stringsAsFactors = FALSE
for (df in drugs) {
  dl <- read.csv2(file=df, header=F)
  #print(dl)
  drugList[df] <- dl
  #drugList["/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_1.txt"][[1]]
  for (dL in drugList[df]) {
    #print(dL)
    print(df)
    drugVars <- c(stringsAsFactors=FALSE)
    for (i in 1:length(dL)) {
      pn <- df
      print(as.character(dL[i]))
      #drugList[df][as.character(dL[i])]
      cui <- getCUI(as.character(dL[i]))
      cuiNew <- str_replace_all(as.character(cui), fixed(" "), "")
      eidUrl <- paste("http://rxnav.nlm.nih.gov/REST/rxcui.json?idtype=UMLSCUI&id=", cuiNew, sep="")
      #print(eidUrl)
      #########print(paste("getExpandedForms(", eidUrl, ",", pn, ",", df, ")", sep=""))
      drugVars <- c()
      variants <- getExpandedForms(eidUrl, as.character(pn), as.character(df))
      print("#### OOMPA LOOMPA #############################################")
      print(variants)
      print("###############################################################")
      for (v in variants) {
        varDrug <- getLukeyCUI(v)
        drugVars <- c(drugVars, varDrug)
      }
      #drugList[df][dL][i] <-
      print("##############DUKEY LUKEYS########################")
      dope <- unlist(unique(drugVars))
      #length(dope)
      ###########print(dope) #collect the cookies - safe in target folders - clean up file names
      #dope[[1]] <- NULL
      #dope
      # at this point, save to disk instead of messing with data structures
      #print(drugList[df][dL][i])
      #fileName = paste("drugLukeys/", df, dL)
      fileName <- paste("/home/smalec/pv/pv_conditions/dataCooked", "/", df, "/drug_", dL[i], ".txt", sep="")
      print(fileName)
      write.table(dope, file=fileName, quote = FALSE, sep = '', row.names = FALSE, col.names = FALSE, append=FALSE)
      system(paste("cat", fileName, sep=" "))
      print("##################################################")
      print("##################################################")
      #length(dope)
      #length(drugVars)
      #mode(dope)
      #mode(drugVars)
    }
  }
}

#getExpandedForms("http://rxnav.nlm.nih.gov/REST/rxcui.json?idtype=UMLSCUI&id=C0248719","/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_1.txt","/home/smalec/pv/PV_conditions/dataRaw/drugs_aki_1.txt")
#getCUI("Ibuprofen")
#cui <- getCUI("Ibuprofen")
#cui



