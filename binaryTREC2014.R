
library(ggplot2)
library(data.table)
#setwd("/data3/cricket/results/TREC2014/trecresults/") # prev: TREC2014
setwd("/home/smalec/pm/results/")
setwd("/db/ball/cricket/results/TREC2014/trecresults")
tCycleMin = 1
tCycleMax = 3

dimInterval = 512
dimStartInterval = 1
dimEndInterval = 96

queryType = list("SUM", "LUCENE")
queryWeighting = list("idf", "logentropy")
queryLevel = list("summaries", "descriptions")

cricketDF <- data.frame(row.names=c("tCycle", "dims", "qType", "qWeighting", "qLevel", "cricket.map", "cricket.gm_ap", "cricket.Rprec", "cricket.bpref"), stringsAsFactors = FALSE)
for (tC in tCycleMin:tCycleMax) {
  for (dI in dimStartInterval:dimEndInterval) {
    for (qT in queryType) {
      for (qW in queryWeighting) {
        for (qL in queryLevel) {
          dimensions = dI*dimInterval #, qL before sep on next line
          fName = paste(paste("tCycle","_",as.character(tC),"_","binaryDim","_",as.character(dimensions),"_",qT,"_",qW,"_", qL, sep=""), ".txt", sep="")
          print(fName)
          if (file.exists(fName)) {
            #try(cricketData <- scan(as.character(fName)))
            #print(fName)
            x <- scan(fName, what="numeric",skip=1, quiet=TRUE, sep="\n")  
            map <- as.numeric(substring(x[4], 21, 26))
            gm_ap <- as.numeric(substring(x[5], 21, 26))
            Rprec <- as.numeric(substring(x[6], 21, 26))
            bpref <- as.numeric(substring(x[7], 21, 26))
            newDF <- data.frame(tCycle=tC, dims=as.character(dimensions), qType=as.character(qT), qWeighting=as.character(qW), qLevel=as.character(qL), cricket.map=map, cricket.gm_ap=gm_ap, cricket.Rprec=Rprec, cricket.bpref=bpref)
            cricketDF <- rbind(cricketDF, newDF)
          }
        }        
      }
    } 
  }
}

print(cricketDF)


######## USED TO GENERATE REPORT GRAPHIC

#cricketDFss <- subset(subset(subset(subset(cricketDF, qLevel == "needs"), tCycle == "2"), qWeighting %in% c("logentropy")))
cricketDF <- subset(subset(subset(cricketDF, qLevel %in% c("descriptions", "summaries")), tCycle == "1"), qWeighting %in% c("logentropy", "idf"))
cricketDF$qType <- as.character(cricketDF$qType)
cricketDF$qType <- replace(cricketDF$qType, cricketDF$qType=="SUM", "LUCENE")
realTREC2014 <- ggplot(cricketDF, aes(x = dims, y = cricket.map, shape = qWeighting, size = qType, color = qType))
realTREC2014 <- realTREC2014 + geom_point(size=5)
realTREC2014 <- realTREC2014 + xlab("dimensions in binary vector space") + ylab("mean average precision (map)") + ggtitle("TREC2014 results: Lucene vs RRI binary vectors") + ylim(0, .18)
#realTREC2014 <- realTREC2014 + facet_wrap(tCycle ~ qLevel) # ~ qWeighting)
realTREC2014



#cricketDFss <- subset(subset(subset(subset(cricketDF, qLevel == "needs"), tCycle == "2"), qWeighting %in% c("logentropy")))
cricketDF <- subset(subset(subset(cricketDF, qLevel %in% c("descriptions", "summaries")), tCycle %in% c("1", "2")), qWeighting %in% c("logentropy"))
cricketDF$qType <- as.character(cricketDF$qType)
cricketDF$qType <- replace(cricketDF$qType, cricketDF$qType=="SUM",  "LUCENE")
realTREC2014 <- ggplot(cricketDF, aes(x = dims, y = cricket.Rprec, shape = qWeighting, size = qType, color = qType))
realTREC2014 <- realTREC2014 + geom_point(size=5)
realTREC2014 <- realTREC2014 + xlab("dimensions in binary vector space") + ylab("mean average precision (map)") + ggtitle("TREC2014 results: Lucene vs RRI binary vectors") + ylim(0, .18)
#realTREC2014 <- realTREC2014 + facet_wrap(tCycle ~ qLevel) # ~ qWeighting)
realTREC2014


ggsave("binaryTREC2014.png", plot = realTREC2014, width = 4, height = 4)

##########

binaryOHSU <- ggplot(cricketDF, aes(x = dims, y = cricket.map, shape = qWeighting, size = qType, color = qType))
binaryOHSU <- binaryOHSU + geom_point(size=list(5, 10))
binaryOHSU <- binaryOHSU + xlab("dimensions") + ylab("mean adjusted precision") + ggtitle("TREC2014 binary")
binaryOHSU <- binaryOHSU + facet_wrap(tCycle ~ qLevel) # ~ qWeighting)
binaryOHSU




binaryOHSU <- ggplot(cricketDF, aes(x = dims, y = cricket.map, shape = qLevel, size = qType, color = qWeighting))
binaryOHSU <- binaryOHSU + geom_point(size=list(300, 500, 1000))
binaryOHSU <- binaryOHSU + xlab("dimensions") + ylab("mean adjusted precision") + ggtitle("TREC2014 binary")
binaryOHSU <- binaryOHSU + facet_wrap(tCycle ~ qLevel)
binaryOHSU


binaryOHSU <- ggplot(cricketDF, aes(x = dims, y = cricket.map, shape = qLevel, size = qType, color = qWeighting))
binaryOHSU <- binaryOHSU + geom_point(200, 500, 900) #size=list(200, 500, 800))
binaryOHSU <- binaryOHSU + xlab("dimensions") + ylab("mean adjusted precision") + ggtitle("TREC2014 binary")
binaryOHSU <- binaryOHSU + facet_wrap(tCycle ~ qLevel) # ~ qWeighting)
binaryOHSU












convert.factors.to.strings.in.dataframe <- function(dataframe)
{
  class.data <- sapply(dataframe, class)
  factor.vars <- class.data[clsas.data == "factor"]
  for (colname in names(factor.vars))
  {
    dataframe[,colname] <- as.character(dataframe[,colname])
  }
  return(dataframe)
}

fileList  = list.files(include.dirs = FALSE)
unlist(fileList)


DF <- data.frame(stringsAsFactors = FALSE)


library(rgl)
# "cycle", "dims", "qType", "qWeighting" 
#print(x[4]) # map x[5] gm_ap   x[6] R-prec  x[7]    bpref 

qtCheck <- cricketDF$qType
qtCol <- c() 
for (i in 1:dim(cricketDF)[1] ) {
  if (qtCheck[i] == 'LUCENE') {
    cricketDF <- cbind(cricketDF, qtCol[i] <- "red")
  } else {
    cricketDF <- cbind(cricketDF, qtCol[i] <- "blue")
  }
  #i = i + 1
}
qtCol

qwCheck <- cricketDF$qWeighting
qwShape <- c() 
for (j in  1:dim(cricketDF)[1]) { 
  if (qwCheck[j] == as.character('logentropy')) { 
    cricketDF <- cbind(cricketDF, qwShape[j] <- "s")
  } else { cricketDF <- cbind(cricketDF, qwShape[j] <- "i")}
}
qwShape

qlCheck <- cricketDF$qLevel
qlFacet <- c()
for (m in 1:dim(cricketDF)[1]){
 if (qlCheck[m] == as.character('robust')) {
   cricketDF <- cbind(cricketDF, qlFacet[m] <- "robust")
 } else  { cricketDF <- cbind(cricketDF, qlFacet[m] <- "skeleton")  }
}

plot3d(x=unlist(cricketDF[1]), y=unlist(cricketDF[2]), z=unlist(cricketDF[5]), xlab="cycle", ylab="dims", zlab="map", col=qtCol, size=0.5, type="s", main="mean average precision", sub="as a function of training cycles and dimensions") 

, type='t')

# "cycle", "dims", "qType", "qWeighting" 
#print(x[4]) # map x[5] gm_ap   x[6] R-prec  x[7]    bpref 
#plot3d(x=unlist(cricketDF[1]), y=unlist(cricketDF[2]), z=unlist(cricketDF[7]), xlab="cycle", ylab="dims", zlab="R-prec")
#plot3d(x=unlist(cricketDF[1]), y=unlist(cricketDF[2]), z=unlist(cricketDF[8]), xlab="cycle", ylab="dims", zlab="bpref")
# N L M funding
# retrieval based  --- not "relevance", but "utility"
# -------> CALEB   |    THORSTEN JOACHIM - SVM-lite
#HAMTMC
# ... queries + things that people have clicked on + downloads
#### how do you use server logs for query result evaluation (use rank order)
# you make up for the sparseness and quality of the data --> MEDLINE (multiple search engines)
# do taxonomies 
# I 2 B 2 
#################### temporal relation extraction
## not generalized to any data set 
# learn on the literature
#
# statistics:          variables
# machine learning:    features
#
#  CANDIDATE + PROJECT ------ deadline is March
# ----------MELISSA
#        NLM training --------- SCMB
#
#

# "cycle", "dims", "qType", "qWeighting" 
#print(x[4]) # map x[5] gm_ap   x[6] R-prec  x[7]    bpref

trainingCycle = cricketDF$cycle
queryWeighting = cricketDF$qWeighting
queryType = cricketDF$qType

# training cycle-map relationship
qplot(x = unlist(cricketDF[1]), y = unlist(cricketDF[5]), xlab = "training cycles", ylab = "map", main = "training cycles-map relationship", color = factor(queryWeighting), shape=factor(queryType), size = factor(trainingCycle))
# dimensions-map relationship
qplot(x = unlist(cricketDF[2]), y = unlist(cricketDF[5]), xlab = "dimensions", ylab = "map", main = "binary dimensions-map relationship", color = factor(queryWeighting), shape=factor(queryType), size=factor(trainingCycle))


# dimensions-gm_ap relationship
p <- qplot(x = unlist(cricketDF[2]), y = unlist(cricketDF[6]), xlab = "dimensions", ylab = "gm_ap", main = "binary dimensions-gm_ap relationship", color = factor(queryWeighting), shape=factor(queryType), size=factor(trainingCycle)) + facet_grid(unlist(cricketDF[5]) ~ .)
# dimensions-R-prec relationship
qplot(x = unlist(cricketDF[2]), y = unlist(cricketDF[7]), xlab = "dimensions", ylab = "R-prec", main = "binary dimensions-R-prec relationship", color = factor(queryWeighting), shape=factor(queryType), size=factor(trainingCycle))
# dimensions-bpref relationship
qplot(x = unlist(cricketDF[2]), y = unlist(cricketDF[8]), xlab = "dimensions", ylab = "bpref", main = "binary dimensions-bpref relationship", color = factor(queryWeighting), shape=factor(queryType), size=factor(trainingCycle))

######## another factor: binary/real
###### time to construct vector spaces
###### performance time on queries
######    ........ more data from trec_eval output? 
#
#
#
#
dimension <- unlist(cricketDF[2][[1]])
dimension
map <- unlist(cricketDF[7][[1]])
map
tCycle <- unlist(cricketDF[1][[1]])
tCycle
library(ggplot2)
sp <- ggplot(cricketDF, aes(x=tCycle, y = map)) + geom_point(shape=1)
sp + facet_grid(. ~ map)

binaryOHSU <- ggplot(cricketDF, aes(x = dims, y = cricket.map, shape = qLevel, size = qType, color = qWeighting))
binaryOHSU <- binaryOHSU + geom_point(size=list(1,3,5))
binaryOHSU <- binaryOHSU + xlab("dimensions") + ylab("mean adjusted precision") + ggtitle("TREC2014 binary")
binaryOHSU <- binaryOHSU + facet_wrap(tCycle ~ qLevel)
binaryOHSU

