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
adm.heuristic <- heuristic(adm.search$nw, adm, adm.prior, restart = 2, degree = 8, trace, TRUE, trylist = adm.search$trylist)
best.adm <- adm.heuristic$nw
