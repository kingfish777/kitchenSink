
#deLukeyCUIize
setwd("/db/mayo")
conditionList <- c("aki", "ali", "gib", "mi")

readVec <- function(xyz) {
  if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
    unique(as.character(unlist(read.table(xyz, stringsAsFactors = F, skipNul = TRUE, header = FALSE))))
  } else { c() }  
}

writeVec <- function(dat,fn) {
  write.table(unique(dat), file=fn, quote=FALSE, row.names=FALSE, eol="\n", col.names=FALSE, append=FALSE)
}

cleanLukeyCUI <- function(x) {
  paste("C", substr(x, start = 7, stop = 13), sep="")
}

for (condition in conditionList) {
  
  fileList <- unlist(list.files(path=condition, pattern="txt", include.dirs = FALSE, recursive=TRUE))
  print("NEW CONDITION")
  print(condition)
  print("##############################")
  for (entity in fileList) {
    print("NEW ENTITY")
    print(entity)
    print("##############################")
    fileName <- paste(condition, "/", entity, sep="")
    data <- readVec(xyz = fileName)
    cleanData <- c()
    for (datum in data) {
      freshDatum <- ""
      print("old")
      print(datum)
      freshDatum <- cleanLukeyCUI(datum)
      print("new")
      print(freshDatum)
      cleanData <- c(cleanData, freshDatum)
    }
    writeVec(dat = cleanData, fn = fileName)
    print(paste("finish with", fileName))
  }
}



