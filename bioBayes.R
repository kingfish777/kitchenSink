
library(deal)
library(bnlearn)
library(stringr)


#df <- table(read.csv2("sjstenData_clean_obs_v2.csv", header = TRUE, sep = ";", stringsAsFactors=FALSE, na.strings = "x")))
#df <- as.data.frame(read.csv("sjstenData_clean_obs_v3.csv", header = TRUE, sep = ";", stringsAsFactors=FALSE))
df <- as.data.frame(read.csv("sj.csv", header = FALSE, sep = "\t", stringsAsFactors=FALSE))
head(df)
df
count(df)
head(df)
#colnames(df) <- c("source", "target", "controlInd", "relOrder")
#colnames(df) <- NULL
df[,1] <- as.factor(df[,1])
df[,2] <- as.factor(df[,2])
df[,3] <- as.factor(df[,3])
df[,4] <- as.factor(df[,4])
df
df <- droplevels(df)
#df[,5] <- as.factor(df[,5])
head(head(df))
#complete.cases(df)
df <- na.omit(df)
x <- df[complete.cases(df), ]


# combinatorial explosion
# divide and conquer
count(substr(df[,2], start = 1, stop = 3) %in% 'adr')
count(substr(df[,2], start = 1, stop = 4) %in% 'drug')
count(substr(df[,2], start = 1, stop = 7) %in% 'protein')
count(substr(df[,2], start = 1, stop = 7) %in% 'pathway')
count(substr(df[,2], start = 1, stop = 4) %in% 'gene')
count(substr(df[,2], start = 1, stop = 4) %in% c('drug', 'drug')) # && count(substr(df[,2], start = 1, stop = 4) %in% 'gene'))

# subset to rows with genes in the first column 
gene <- as.data.frame(df[grep("gene", df[,1]),])
# subset to rows with drugs in the second column
drug <- as.data.frame(df[grep("drug", df[,2]),])
# subset data frame to second rcolumns
adr <- as.data.frame(df[grep("adr", df[,2]),])

adm <- merge.data.frame(x = adr, y = drug, by = 1)
# try others

adm.nw <- network(adm)
adm.prior <- jointprior(adm.nw)
adm.nw.v2 <- learn(adm.nw, adm, adm.prior)$nw
adm.search <- autosearch(adm.nw.v2, adm, adm.prior, trace = TRUE)
adm.heuristic <- heuristic(adm.search$nw, adm, adm.prior, restart = 2, degree = 10, trace, TRUE, trylist = adm.search$trylist)

adm.heuristic <- heuristic(adm.search$nw, adm, adm.prior, trace = TRUE, trylist = adm.search$trylist)

best.adm <- adm.heuristic$nw




gene <- na.omit(df)
gene.nw <- network(gene)
gene.prior <- network(gene.nw)
#writeVec(dat=df, fn="dat.csv")
#raw case list to aggregated case list
#df <- as.data.frame(ftable(df))
#df
df.nw <- network(df)

df.prior <- jointprior(df.nw)

df.nw.v2 <- learn(df.nw, df, df.prior)$nw
# perform structural search
df.search <- autosearch(df.nw.v2, df, df.prior, trace=TRUE)
df.heuristic <- heuristic(df.search$nw, df, df.prior, restart = 2, degree = 10, trace = TRUE, trylist = df.search$trylist)
best.df <- df.heuristic$nw

buildBN <- function(x) {
  x <- gene
  x <- na.omit(x)
  #x <- as.data.frame(x)
  x.nw <- network(x)
  x.prior <- jointprior(x.nw)
  x.nw.v2 <- learn(x.nw, x, x.prior)$nw
  x.search <- autosearch(x.nw.v2, x, x.prior, trace=TRUE)
  x.heuristic <- heuristic(x.search$nw, x, x.prior, restart = 2, degree = 10, trace = TRUE, trylist = x.search$trylist)
  best.x <- x.heuristic$nw
}

buildBN(gene)


#df.tab2 <- xtabs(Freq~., data=df.tab)
#aggregated case list to table
xtabs(Freq~., data=df)
pv <- df.v2

pv.nw <- network(pv)

pv.prior <- jointprior(pv.nw)
pv.nw.v2 <- learn(pv.nw, pv, pv.prior)$nw
# perform structural search
pv.search <- autosearch(pv.nw.v2, pv, pv.prior, trace=TRUE)
pv.heuristic <- heuristic(pv.search$nw, pv, pv.prior, restart = 2, degree = 10, trace = TRUE, trylist = pv.search$trylist)
best.pv <- pv.heuristic$nw


factor.nc <- ncol(pv)-3
print(factor.nc)
confoundingAndMedNames <- names(pv[,1:factor.nc])
print(confoundingAndMedNames)
for (pv.df.i in 1:factor.nc) {
  print(pv.df.i)
  pv[,pv.df.i] <- as.factor(pv[,pv.df.i])
  #print(confoundingAndMedNames[pv.df.i])
  #fieldName <-  substr(x = confoundingAndMedNames[pv.df.i], start = 7, stop = nchar(confoundingAndMedNames[pv.df.i]))
  #print(fieldName)
  #names(pv[,pv.df.i]) <- fieldName 
}
head(pv)




pv.nw <- network(pv)

pv.prior <- jointprior(pv.nw)
pv.nw.v2 <- learn(pv.nw, pv, pv.prior)$nw
# perform structural search
pv.search <- autosearch(pv.nw.v2, pv, pv.prior, trace=TRUE)
pv.heuristic <- heuristic(pv.search$nw, pv, pv.prior, restart = 2, degree = 10, trace = TRUE, trylist = pv.search$trylist)
best.pv <- pv.heuristic$nw
