

setwd("/db/datafarm/d3jss")
homeDir <- getwd()

conditions <-  c("gib") #"mi", "ali", "gib", "aki") #, "aki") # ali", "gib", "aki", "mi")
levs <- c("0") #, "1") #c("1", "0")
#pathways <- c("CAUSES-INV==", "PREDISPOSES-INV==", "TREATS==COEXISTS_WITH-INV")
#pathways <- c("TREATS==COEXISTS_WITH-INV", "PREDISPOSES-INV==")
pathways <- c("TREATS==COEXISTS_WITH-INV") ##"PREDISPOSES-INV==", "CAUSES-INV==", "TREATS==COEXISTS_WITH-INV")
cTypes <- c("REAL") #, "EXCLUSIONS")


generateWordCloud <- function(conditions, pathways) {
  
  for (c in conditions) {
    for (l in levs) {
      for (p in pathways) {
        sp1 ="/"
        path <- paste(c,sp1,l,sp1,p,sp1,sep="")
        print(path)
        drugs <- list.files(path = path, include.dirs = FALSE)
        print(drugs)
        #drugs <- c("Acarbose", "Adenosine", "benzonatate", "ferrous_gluconate", "Nevirapine", "orlistat", "Paromomycin", "Phentermine", "pioglitazone", "rosiglitazone", "Scopolamine", "Stavudine", "Tetrahydrocannabinol", "Urea")
        drugs <- c("abacavir", "benzonatate", "Dicyclomine", "fluticasone", "Ketoconazole", "Lactulose", "Loratadine", "metaxalone", "Methocarbamol", "Nitrofurantoin", "pioglitazone", "rosiglitazone", "salmeterol", "Temazepam")
        sp2="__"
        dumpPathREAL <- paste(c,sp2,l,sp2,p,sp2,"REAL", sep="")
        system(paste("rm -r ", dumpPathREAL, sep=""))
        system(paste("mkdir ", dumpPathREAL, sep=""))
        print(dumpPathREAL)
        dumpPathEXCLUSIONS <- paste(c,sp2,l,sp2,p,sp2,"EXCLUSIONS", sep="")
        print(dumpPathEXCLUSIONS)
        system(paste("rm -r ", dumpPathEXCLUSIONS, sep=""))
        system(paste("mkdir ", dumpPathEXCLUSIONS, sep=""))
        system(paste("ls ", dumpPathEXCLUSIONS, sep=""))
        drugs <- 
        for (d in drugs) {
          cReal <- paste("cp ", path, d, sp1, "dataREAL.txt ", dumpPathREAL, sp1, d, ".txt", sep="")
          print(cReal)
          system(cReal)
          # cEXC <- paste("cp ", path, d, sp1, "dataEXCLUSIONS.txt ", dumpPathEXCLUSIONS, sp1, d, ".txt", sep="")
          #system(cEXC)
          #print(cEXC)
        }
      }
    }
  }
  system(paste("ls gib__1__CAUSES-iNV==__EXCLUSIONS", sep=""))
  
  # clean
  
  readVec <- function(xyz) {
    if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
      unique(as.character(unlist(read.table(xyz, stringsAsFactors = F, skipNul = TRUE, header = FALSE))))
    } else { c() }  
  }
  
  writeVec <- function(dat,fn) {
    write.table(unique(dat), file=fn, quote=FALSE, row.names=FALSE, eol="\n", col.names=FALSE, append=FALSE)
  }
  
  getwd()
  for (c in conditions) {
    for (l in levs) {
      for (p in pathways) {
        for (t in cTypes) {
          sp0 <- "/"
          sp3 <- "|"
          sp4 <- "__"
          targetDir <- paste(homeDir,sp0,c,sp4,l,sp4,p,sp4,t,sp0, sep="")
          print(targetDir)
          drugs <- list.files(path=targetDir, recursive=TRUE)
          print(d)
          for (d in drugs) {
            tFile <- paste(homeDir,sp0,c,sp4,l,sp4,p,sp4,t,sp0,d, sep="")
            print(tFile)
            system(paste("cat ", tFile, sep=""))
            try(dat <- readVec(tFile))
            writeVec(dat = dat, fn = tFile)
          }
        }
      }
    }
  }
  
  
  
  library("tm")
  library("SnowballC")
  library("wordcloud")
  library("RColorBrewer")
  
  #getWordcloud <- function(conditions, cTypes, levs, pathways) {
  
  # create wordclouds
  for (c in conditions) {
    for (l in levs) {
      for (p in pathways) {
        for (t in cTypes) {
          sp0 <- "/"
          sp3 <- "|"
          sp4 <- "__"
          targetDir <- paste(homeDir,sp0,c,sp4,l,sp4,p,sp4,t,sp0, sep="")
          print(targetDir)
          setwd(targetDir)      
          text <- system.file("texts", "txt", package="tm")
          # read in corpus
          cname <- file.path(".")
          corpus <- Corpus(DirSource(cname))
          #print(corpus[[3]])
          #corpus <- tm_map(FUN = removeNumbers, corpus)
          corpus <- tm_map(FUN = removePunctuation, corpus)
          #corpus <- tm_map(FUN = stemDocument, corpus)
          corpus <- tm_map(corpus, removeWords, stopwords(kind = "en"))
          #corpus <- tm_map(corpus, removeWords, stopwords(c("x")))
          
          #corpus[[1]]
          ## <<PlainTextDocument>>
          ## Metadata:  7
          ## Content:  chars: 558
          inspect(corpus)
          ## <<VCorpus>>
          ## Metadata:  corpus specific: 0, document level (indexed): 0
          ## Content:  documents: 32
          dtm <- TermDocumentMatrix(corpus)
          dtm
          m <- as.matrix(dtm) 
          m
          v <- sort(rowSums(m),decreasing=TRUE)
          v
          d <- data.frame(word = names(v),freq=v)
          
          #set.seed(1234)
          #setwd(homeDir)
          fn <- paste(homeDir,sp0,c,sp4,l,sp4,p,sp4,t,".jpg", sep="")
          print(fn)
          #jpeg(fn, units = "px", width = 1900, height = 1900)
          wordcloud(use.r.layout = TRUE, words = d$word, freq = d$freq, min.freq = 1,
                    random.order=FALSE, rot.per=0.35, #, 
                    colors=brewer.pal(8, name = "BrBG"))
          #dev.off()
        }
      }
    }
  }
  
}

#conditions <-  c("gib") #"mi", "ali", "gib", "aki") #, "aki") # ali", "gib", "aki", "mi")
#levs <- c("0", "1") #c("1", "0")
#pathways <- c("CAUSES-INV==", "PREDISPOSES-INV==", "TREATS==COEXISTS_WITH-INV")
#pathways <- c("TREATS==COEXISTS_WITH-INV", "PREDISPOSES-INV==")
#pathways <- c("TREATS==COEXISTS_WITH-INV") ##"PREDISPOSES-INV==", "CAUSES-INV==", "TREATS==COEXISTS_WITH-INV")
#cTypes <- c("REAL") #, "EXCLUSIONS")

#generateROCs_c <- compiler::cmpfun(generateROCs)

args <- commandArgs(trailingOnly=TRUE)
adr <<- as.character(args[1])
pred <<- as.character(args[2])
metric <<- as.character(args[3])

#if (length(args) > 0){  generateROCs_c(adr, pred, metric) } else {  generateROCs_c(gib, tc, m)   }

generateWordCloud_c <- compiler::cmpfun(generateWordCloud)
#if (length(args) > 0){  generateROCs_c(adr, pred, metric) } else {  generateROCs_c(gib, c, m)   }

generateWordCloud_c(conditions = "aki", pathways = "CAUSES-INV==")
#generateWordCloud_c(conditions = "aki", pathways = "PREDISPOSES-INV==")
generateWordCloud_c(conditions = "aki", pathways = "TREATS==COEXISTS_WITH-INV")

generateWordCloud_c(conditions = "ali", pathways = "CAUSES-INV==")
generateWordCloud_c(conditions = "ali", pathways = "PREDISPOSES-INV==")
generateWordCloud_c(conditions = "ali", pathways = "TREATS==COEXISTS_WITH-INV")

generateWordCloud_c(conditions = "gib", pathways = "CAUSES-INV==")
generateWordCloud_c(conditions = "gib", pathways = "PREDISPOSES-INV==")
generateWordCloud(conditions = "gib", pathways = "TREATS==COEXISTS_WITH-INV")

generateWordCloud_c(conditions = "mi", pathways = "CAUSES-INV==")
generateWordCloud_c(conditions = "mi", pathways = "PREDISPOSES-INV==")
generateWordCloud_c(conditions = "mi", pathways = "TREATS==COEXISTS_WITH-INV")


