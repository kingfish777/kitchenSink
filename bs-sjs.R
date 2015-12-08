library(deal)
library(stringr)

#fN <- "/home/kingfish/datafarm/d3jss/gib/1/TREATS==COEXISTS_WITH-INV/Citalopram/dataEXCLUSIONS.txt"

#datum <- readVec("/home/kingfish/datafarm/d3jss/gib/1/TREATS==COEXISTS_WITH-INV/Citalopram/dataEXCLUSIONS.txt")
#makeUnique(fN)
#system("cat /home/kingfish/datafarm/d3jss/gib/1/TREATS==COEXISTS_WITH-INV/Citalopram/dataEXCLUSIONS.txt")

readVecSep <- function(xyz, sp) {
  if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
    unique(as.character(unlist(read.table(xyz, stringsAsFactors = F, sep = sp, skipNul = TRUE, header = FALSE, allowEscapes = TRUE))))
  } else { c() }  
}



readVec <- function(xyz) {
  if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
    unique(as.character(unlist(read.table(xyz, stringsAsFactors = F, skipNul = TRUE, header = FALSE))))
  } else { c() }  
}

readVecT <- function(xyz) {
  if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
    unique(as.character(unlist(read.table(xyz, stringsAsFactors = F, skipNul = TRUE, header = TRUE))))
  } else { c() }  
}

readVecSepT <- function(xyz, sp) {
  if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
    unique(as.character(unlist(read.table(xyz, stringsAsFactors = F, skipNul = TRUE, header = TRUE))))
  } else { c() }  
}

writeVec <- function(dat,fn) {
  write.table(unique(dat), file=fn, quote=FALSE, row.names=FALSE, eol="\n", col.names=FALSE, append=FALSE)
}

makeUnique <- function(fName) {
  data <- readVec(fName)
  writeVec(fn = fName, dat = data)
}


#writeVec(dat = dat, fn = "/home/kingfish/datafarm/d3jss/gib/1/TREATS==COEXISTS_WITH-INV/Citalopram/dataEXCLUSIONS.txt")

setwd("/home/kingfish/UTH/getPathwaysStandalone/")

adrs <- c("sjsten")

#adrs <- c("gib", "ali", "mi", "sjsten", "aki")

# read in ADR 

# write to doc
#adr <- "aki"
#command <- paste("sudo ./qdis.sh ", adr, sep="")
#print(command)
#system(command)
#adrProteins.fn <- paste(adr, "Proteins.txt", sep="")
#print(adrProteins.fn)
#command <- paste("sudo ./sedDisease.sh ", adr, " > ",  adrProteins.fn, sep="")
#print(command)
#system("cat gibProteins.txt")
#adrProteins <- readVecSep(adrProteins.fn, sp=";")
#print(adrProteins)
#writeVec(dat = adrProteins, fn = adrProteins.fn)
#system(paste("sudo ./clean.sh ", adr, " > ", adr, "Proteins.txt", sep=""))

#system("cat gibProteins.txt")
#readVecSep(adrProteins.fn, sp=";")

getIntersection <- function(fn1, fn2, fn3) { #}, fn4) {
  system(paste("java -jar ADRPathways.jar ", fn1, " ", fn2, " ", fn3, " datafarm datafarm", sep=""))#  2>&1 ", fn4, sep="")) # 2>&1 | tee -a ", fn4, sep=""))
}


controlInds <- c("1", "0")

for (adr in adrs) {
  # write each doc to file
  try(system(paste("sudo rm proteinPathways.txt")), silent = TRUE)
  try(system(paste("sudo rm drugsnps.txt")), silent = TRUE)
  try(system(paste("sudo rm drugProteinPathways.txt")), silent = TRUE)
  try(system(paste("sudo rm drugUniProtID.txt")), silent = TRUE)
  try(system(paste("sudo rm drugTargets.txt")), silent = TRUE)
  #getIntersection("drugProteinPathways.txt", adrProteins.fn, "output.txt")
  adrProteins.fn <- paste(adr, "Proteins.txt", sep="")
  print(adrProteins.fn)
  for (ci in controlInds) {
    # read in drug list
    #print(ci)
    drugs <- readVec(paste(adr, "_drug", ci, ".txt", sep=""))
    print(drugs)
    for (drug in drugs) {
      print(drug)
      system(paste("mkdir ", adr, "/", ci, "/", drug, sep=""))
      print("bunga bunga")
      #gdid #hits runADR
      command <- paste("sudo ./gdid.sh ", drug, " > drugUniProtID.txt", sep="")
      print(command)
      system(command)
      system("cat drugUniProtID.txt ")
      uniProtID.drug <- readVec("drugUniProtID.txt")
      print(uniProtID.drug)
      
      
      tryCatch({command <- paste("./gdpathways.sh ", uniProtID.drug, " > drugProteinPathways.txt ", sep="")
                print(command)
                system(command)
                
                if (file.exists("drugProteinPathways.txt")) {
                  print(file.info("drugProteinPathways.txt"))
                }
                if (file.exists("drugProteinPathways.txt") && file.info("drugProteinPathways.txt")$size > 0) {
                  print("drugProteinPathways.txt exists")
                  print(paste("checking for regulatory protein Pathways for ADR: ", adr, sep=""))
                  try(drugPathways <- read.table("drugProteinPathways.txt", sep="\t", header = FALSE)[,3:5])
                  uniProtIDs.proteinPathways <- drugPathways[2]
                  print(uniProtIDs.proteinPathways)
                  try(writeVec(dat = uniProtIDs.proteinPathways, fn = "drugProteinPathways.txt"))
                  drugPathways <- readVec("drugProteinPathways.txt")
                  output.fn <- paste(adr, "/", ci, "/", drug, "/proteinPathway.txt", sep="")
                  #outputDat.fn <- paste(adr, "/", ci, "/", drug, "/proteinPathway.dat", sep="")
                  lookup.fn <- paste(adr, "/", ci, "/", drug, "/drugProteinPathways_master.txt", sep="")
                  writeVec(dat = drugPathways, fn = lookup.fn)
                  print("here i am in protein pathway land")
                  print(drug)
                  print("MONKNEY TWADDDLE")
                  print(" bad data? ")
                  system(paste("grep -i 'uniprot' drugProteinPathways.txt "))
                  getIntersection("drugProteinPathways.txt", adrProteins.fn, output.fn)
                  getIntersection("drugProteinPathways.txt", adrProteins.fn, "proteinPathways.txt") #, outputDat.fn)
                  system(paste("cp ", "proteinPathways.txt ", adr, "/", ci, "/", drug, "/.", sep=""))
                  
                } else { }
                
                
      }, warning = function(w) {
        #warning-handler-code
      }, error = function(e) {
        #error-handler-code
      }, finally = {
        #cleanup-code
      })
      
      
      
      tryCatch({print("#######################################################")
                print(paste("checking for target proteins, pathways, genes in common with ", adr, sep=""))
                print("#######################################################")
                
                command <- paste("./gdtarget.sh ", uniProtID.drug, " > drugTargets.txt ", sep="")
                print(command)
                system(command)
                if (file.exists("drugTargets.txt")) {
                  print(file.info("drugTargets.txt"))
                  print(read.table("drugTargets.txt"))
                  
                  
                  print("WWOOOOOHOOOO")
                }
                
                if (file.exists("drugTargets.txt") && file.info("drugTargets.txt")$size > 2) {
                  #try(drugTargets <- read.table("drugTargets.txt", sep="\t", header = FALSE)[,3:5])
                  #uniProtIDs.targets <- drugTargets[2]
                  #print(uniProtIDs.targets)
                  #try(writeVec(dat = uniProtIDs.targets, fn = "drugTargets.txt"))
                  
                  drugTargets <- readVec("drugTargets.txt") #read.table("drugTargets.txt", header = FALSE, stringsAsFactors = FALSE)
                  #geneTargets <- drugTargets
                  print(drugTargets)
                  
                  #try(writeVec(dat = geneTargets, fn = "drugTargets.txt"))
                  #readVec("drugTargets.txt")
                  for (gene in drugTargets) {
                    adrCui.fn <- paste(adr, "Cui.txt", sep="")
                    adrcuis <- readVec(adrCui.fn)
                    for (ac in adrcuis) { 
                      print(ac)
                      print(gene)
                      output.fn <- paste(adr, "/", ci, "/", drug, "/targetgene_", gene, ".txt", sep="")
                      outputDat.fn <- paste(adr, "/", ci, "/", drug, "/targetgene_", gene, ".dat", sep="")
                      lookup.fn <- paste(adr, "/", ci, "/", drug, "/targetgene_dgi_", gene, ".txt", sep="")
                      #writeVec(dat = drugTargets, fn = lookup.fn)
                      #grep PTGS2 befree_gene_disease_associations.txt | grep Myocardial
                      geneDrug.fn <- paste(adr, "/", ci, "/", drug, "/targetgene_", gene, ".txt", sep="")
                      system(paste("grep '", str_trim(gene, side = c("both")), "' befree_gene_disease_associations.txt | grep '", str_trim(ac, side = c("both")), "' > ", geneDrug.fn, sep=""))
                      print(geneDrug.fn)
                      print("gene drug target data")
                      print("--------------")
                      system(paste("cat ", geneDrug.fn, sep=""))
                      print("--------------")
                    }
                  }
                } else { }
                
      }, warning = function(w) {
        #warning-handler-code
      }, error = function(e) {
        #error-handler-code
      }, finally = {
        #cleanup-code
      })
      
      tryCatch({
        
        print(paste("checking for SNPs associated with ADVERSE DRUG REACTIONs for ADR: ", adr, sep=""))
        command <- paste("./gdsnp.sh ", uniProtID.drug, " > drugsnps.txt ", sep="")
        print(command)
        system(command)
        if (file.exists("drugsnps.txt")) {
          print(file.info("drugsnps.txt")) 
        }
        if (file.exists("drugsnps.txt") && file.info("drugsnps.txt")$size > 0) {
          snpADRs <- read.table("drugsnps.txt", sep="\t", header = FALSE)[,3:5]
          uniProtIDs.snpADRs <- snpADRs[2]
          print(uniProtIDs.snpADRs)
          length(uniProtIDs.snpADRs)
          writeVec(dat = uniProtIDs.snpADRs, fn = "drugsnps.txt")
          output.fn <- paste(adr, "/", ci, "/", drug, "/proteinSNPS.txt", sep="")
          #outputDat.fn <- paste(adr, "/", ci, "/", drug, "/proteinSNPS.dat", sep="")
          lookup.fn <- paste(adr, "/", ci, "/", drug, "/drugSNPS.txt", sep="") ### context data for surviving ixn SNPs
          writeVec(dat = snpADRs, fn = lookup.fn)
          print("here i am")
          print("MONKNEY TWADDDLE")
          print(" bad data? ")
          system(paste("grep -i 'uniprot' drugsnps.txt "))
          getIntersection("drugsnps.txt", adrProteins.fn, output.fn) #, outputDat.fn)
        } else { }
        #system(paste("sudo rm ", adrProteins.fn, sep=""))
        
        
      }, warning = function(w) {
        #warning-handler-code
      }, error = function(e) {
        #error-handler-code
      }, finally = {
        #cleanup-code
      })
      
      
      
      try(system(paste("sudo rm proteinPathways.txt")), silent = TRUE)
      try(system(paste("sudo rm drugsnps.txt")), silent = TRUE)
      try(system(paste("sudo rm drugProteinPathways.txt")), silent = TRUE)
      try(system(paste("sudo rm drugUniProtID.txt")), silent = TRUE)
      try(system(paste("sudo rm drugTargets.txt")), silent = TRUE)
    }
  }
}
