setwd("/home/kingfish/statnetTutorial")
#  /Users/smalec/statnetTutorial
 # install.packages("RCurl")
 # nstall.packages("ergm")
library(RCurl); library(ergm)

### STEP 1 
 # system("ls *.tsv")
#First, read in the sociomatrix
 # ga.mat<-getURL("https://docs.google.com/spreadsheet/pub?key=0Ai--oOZQWBHSdDE3Ynp2cThMamg1b0VhbEs0al9zV0E&single=true&gid=0&output=txt",
    #           ssl.verifypeer = FALSE)
 # str(ga.mat)
 # print(ga.mat)
 # ga.mat<-getURL("grey_adjacency.tsv",
 #              ssl.verifypeer = FALSE)
 # xmlTreeParse("grey_adjacency.tsv", useInternalNodes = TRUE)
 # str(ga.mat)

ga.mat <- as.matrix(read.csv2("greys_adj_mat.tsv",quote = "\"", header = TRUE, sep = "\t", row.names=1))

 # ga.mat<-as.matrix(read.table(textConnection(ga.mat), sep="\t", 
 #                            header=T, row.names=1, quote="\", stringsAsFactors = F))
ga.mat # SUCCESS!!!!

#ga.adj<-as.matrix(read.table(('grey_adjacency.tsv'), sep="\t",
#                             header=T, row.names=1, quote="\""))


library(statnet)

# Now to make the data into a network called grey.net



ga.atts <- as.matrix(as.table(read.csv("nodes.tsv", header = TRUE, sep = "\t", row.names=1, stringsAsFactors = FALSE)))
#ga.atts <- as.matrix(read.csv("nodes.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE))
#system("ls *.tsv")
ga.atts


grey.net.a <- network( ga.mat, vertex.attr=ga.atts, vertex.attrnames=rownames(ga.atts), directed=FALSE, hyper=FALSE, loops=FALSE, multiple=FALSE, bipartite=FALSE )
plot(grey.net.a)
#ga.attr <- ga.attrs
#length(vertex.attr)
#length(vertex.attrnames)
#Second, read in the network attributes
#ga.atts<-getURL("https://docs.google.com/spreadsheet/pub?key=0Ai--oOZQWBHSdDE3Ynp2cThMamg1b0VhbEs0al9zV0E&single=true&gid=1&output=txt",
 #               ssl.verifypeer = FALSE)
#ga.atts<-read.table(textConnection(ga.atts), sep="\t", header=T, quote="\"", stringsAsFactors=F, strip.white=T, as.is=T)
#str(ga.mat)
#str(ga.atts)
#Third, create a network object using the sociomatrix and its corresponding attributes
ga.net <- network(ga.mat, vertex.attr=ga.atts, vertex.attrnames=rownames(ga.atts), directed=F, hyper=F, loops=F, multiple=F, bipartite=F)
#ga.net <- network(ga.mat, vertex.attr=ga.atts)
### BREAKDOWN =====> Error in vertex.attrnames[[i]] : subscript out of bounds

plot(ga.net, vertex.col=c("blue","pink")[1+(get.vertex.attribute(ga.net, "sex")=="F")], label=get.vertex.attribute(ga.net, "name"), label.cex=.75)

ga.base<-ergm(ga.net~edges+nodematch("sex")) #Estimate the model
summary(ga.base) #Summarize the model

plot(simulate(ga.base),
     vertex.col=c("blue","pink")[ 1+(get.vertex.attribute(ga.net, "sex")=="F")])

ga.base.gof<-gof(ga.base)
summary(ga.base.gof) #Summarize the goodness of fit
par(mfrow=c(3,1)); plot(ga.base.gof) #Plot three windows. It's OK to ignore the warning.

ga.base.d1<-ergm(ga.net~edges+nodematch("sex")+degree(1))
summary(ga.base.d1)
plot(simulate(ga.base.d1),
     vertex.col=c("blue","pink")[1+(get.vertex.attribute(ga.net, "sex")=="F")])
ga.base.d1.gof<-gof(ga.base.d1); summary(ga.base.d1.gof)
par(mfrow=c(3,1)); plot(ga.base.d1.gof)

mcmc.diagnostics(ga.base.d1)

ga.base.d1<-ergm(ga.net~edges+nodematch("sex")+degree(1),
                 control=control.ergm(MCMC.burnin=50000, MCMC.interval=5000))

ga.base.d1<-ergm(ga.net~edges+nodematch("sex")+degree(1),
                 control=control.ergm(MCMC.burnin=50000,
                                      MCMC.interval=5000, MCMC.samplesize=50000))

summary(ga.base.d1) #The figures are mostly the same
mcmc.diagnostics(ga.base.d1)

ga.base.d1.age<-ergm(ga.net~edges+nodematch("sex")+degree(1)+absdiff("birthyear"),
                     control=control.ergm(MCMC.burnin=50000, MCMC.interval=5000))
summary(ga.base.d1.age)
mcmc.diagnostics(ga.base.d1.age)
summary(gof(ga.base.d1.age))


ga.base.d1.age.race<-ergm(ga.net~edges+nodematch("sex")+degree(1)
                          +absdiff("birthyear")+nodematch("race"),
                          control=control.ergm(MCMC.burnin=50000,MCMC.interval=5000))
summary(ga.base.d1.age.race)
mcmc.diagnostics(ga.base.d1.age.race)
summary(gof(ga.base.d1.age.race))
#We do find a positive and significant racial assortativity effect, suggesting that tie formation is more likely between characters of the same race.  In this show, sexual relationships are more likely to be intraracial than interracial compared to chance expectations.  The MCMC diagnostics look good and the goodness of fit measures are in the ballpark.  Though the model looks good by most indicators, the BIC backslides from 297.51 in the previous model on age to 298.21 in this model.  If we’re interested in a simple explanation characterizing network structure, the previous model on age provides the better answer.

#But is it possible that the white characters are homophilous and the black characters less so?  Or vice-versa?  We can refer to this process as differential homophily.  We can model it by adding the parameter “diff=T” to the racial nodematch term.

ga.base.d1.age.racediff<-ergm(ga.net ~ edges + nodematch("sex") + degree(1)
                              +absdiff("birthyear") + nodematch("race", diff=T),
                              control=control.ergm(MCMC.burnin=50000, MCMC.interval=5000))
summary(ga.base.d1.age.racediff)


ga.base.d1.age.racediff<-ergm(ga.net~edges+nodematch("sex")+degree(1)
                              +absdiff("birthyear")+nodematch("race", diff=T, keep=c(1,3)),
                              control=control.ergm(MCMC.burnin=50000, MCMC.interval=5000))
summary(ga.base.d1.age.racediff)
mcmc.diagnostics(ga.base.d1.age.racediff)
summary(gof(ga.base.d1.age.racediff))


plot(simulate(ga.base.d1.age.racediff),
     vertex.col=c("blue","pink")[1+(get.vertex.attribute(ga.net, "sex")=="F")])

#Perhaps men and women have different tendencies to form sexual partnerships?
ga.base.d1.age.sex<-ergm(ga.net~edges+nodematch("sex")+degree(1)
                         +absdiff("birthyear")+nodefactor("sex"),
                         control=control.ergm(MCMC.burnin=100000, MCMC.interval=5000))
#Maybe less "traditional" than previously concluded?
#Perhaps you're interested in the assortative mixture among roles in the hospital?
#Resident-resident sexual contacts (28) are the reference group,
#and unobserved pairings between positions (3-5, 9, 12-15, 17-21, 24, 27) are omitted.
ga.base.d1.age.rolemix<-ergm(ga.net~edges+nodematch("sex")+degree(1)+absdiff("birthyear")
                             +nodemix("position", base=c(3:5, 9, 12:15, 17:21, 24, 27, 28)),
                             control=control.ergm(MCMC.burnin=50000, MCMC.interval=5000))
list.vertex.attributes(ga.net) #List the different vertex attributes in the network object
get.vertex.attribute(ga.net, "sign") #Provides each actor's astrological sign.
?ergm.terms #This command will provide a list of other terms that ergm() can model
?control.ergm #This command will present different options on the estimation methods.


