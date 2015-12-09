

library(reactome.db)
library(UniProt.ws)
library(stringr)

reactome_dbconn()

limit = 20

setwd("/home/kingfish/UTH/getPathwaysStandalone/")

getJeans <- function(adr, ci, med, description, pid) {
  # description: SNPS or PTPW
    # PTPW: get genes from pathways in ixnProteinPathway.txt in folder
    # SNPS: get genes from pathways in proteinSNPS.txt in folder
  
  genefilename.fn <- paste(adr, "/", ci, "/", med, "/", description, "genes_", pathID, ".txt", sep="")
  
  queryReactomeForGenesWithPathID <- function(pathID, reslimit) {
    queryString = paste("SELECT gene_id FROM pathway2gene WHERE DB_ID = ", pathID, " LIMIT ", reslimit, sep="")
    genesFromPathwayIDs <- dbGetQuery(reactome_dbconn(), queryString)
    return(genesFromPathwayIDs)
  }

  # get pathway id proteinpathway and SNP output (field 1)

  keys <- as.character(queryReactomeForGenesWithPathID(pid, limit))

  up <- UniProt.ws(taxId=9606)
  #up
  species(up)
  taxId(up) <- 9606

  if (interactive()) {
    egs = keys(up, keytype = "ENTREZ_GENE")
  }

  columns <- c("UNIREF100", "HGNC") #SEQUENCE
  kt <- "ENTREZ_GENE"
  genes.UniProt <- select(up, keys, columns, kt)
  #genes.UniProt[2]
  

  cleanUniProtGeneNames <- function(upgName) {
    substr(x = upgName, start = 11, stop = 100)
  }

  upgNames <- sapply(genes.UniProt[2], cleanUniProtGeneNames)
  #description = "FUNGI"
  for (gene in upgNames[2:length(upgNames)]) {
    print(gene)
    system(paste("./sedUP2HGNCsymb.sh ", gene, " >> ", genefilename.fn, sep="")) # to correct place
  }

  print(system(paste("cat", genefilename.fn)))

}

readVec <- function(xyz) {
  if (file.exists(xyz) && (file.info(xyz)$size > 0)) {
    #unique(as.character(unlist(read.table(xyz, stringsAsFactors = F, skipNul = TRUE, header = TRUE, sep = "\t"))))
    unique(as.character(unlist(read.table(xyz, stringsAsFactors = F, header = TRUE, sep = "\t"))))
    
  } else { c() }  
}

writeVec <- function(dat,fn) {
  write.table(unique(dat), file=fn, quote=FALSE, row.names=FALSE, eol="\n", col.names=FALSE, append=FALSE)
}

adrs <- c("sjsten")

controlInds <- c("1", "0")

for (adr in adrs) {
  # write each doc to file
  #adrProteins.fn <- paste(adr, "Proteins.txt", sep="")
  #print(adrProteins.fn)
  for (ci in controlInds) {
    # read in drug list
    #print(ci)
    drugs <- readVec(paste(adr, "_drug", ci, ".txt", sep=""))
    print(drugs)
    for (drug in drugs) {
      folderPath <- paste(adr, "/", ci, "/", drug, sep="")
      print(folderPath)
      
      tryCatch({print(" compiling list of target genes from drug bank for med ")
      try(tg <- as.character(list.files(path = folderPath, pattern = "targetgene.")))
      print(tg)
      rtg.fn <- paste(folderPath, "/drugGeneTargets.txt", sep="")
      
      if (length(tg) > 0) {
        targetGenes <- as.character(list.files(path = folderPath, pattern = "targetgene."))
        print(targetGenes)
        cleanTG <- function(s) {
          s <- substring(text = s, first = 12, last = str_locate(string = s, pattern = "[.]")-1)
          #return(s)    
        }
        targetGenes.new <- as.character(sapply(targetGenes, cleanTG))
        targetGenes.new
        writeVec
        writeVec(fn = rtg.fn, dat = targetGenes.new)
      }
        }, warning = function(w) {
        #warning-handler-code
      }, error = function(e) {
        #error-handler-code
      }, finally = {
        #cleanup-code
      })
      
      tryCatch({print("performing drug-ADR protein pathway lookup for related genes")
      # PTPW: get genes from pathways in ixnProteinPathways.txt in folder
      ptpw.fn <- paste(folderPath, "/ixnProteinPathways.txt", sep="")
      PTPWproteins.fn <- paste(folderPath, "/PTPWproteins.txt", sep="")
      #print(file.info(ptpw.fn))
      if (file.exists(ptpw.fn) && file.info(ptpw.fn)$size > 0) {
        pidsPP.fn <- paste(folderPath, "/pidsPP.txt", sep="")
        system(paste("./sedGetProteins.sh ", ptpw.fn, " >  ", PTPWproteins.fn, sep=""))
        #system(paste("./sedGetProteins.sh ", ptpw.fn, " >  ", proteins.fn, sep=""))
        system(paste("sudo ./sedGetPathways.sh ", ptpw.fn, " > ", pidsPP.fn,  sep=""))
        pids <- readVec(pidsPP.fn)
        print(pids)
        for (pid in pids) {
          print(pid)
          #getJeans(adr, ci, med, description, pid)
        }
      }
      }, warning = function(w) {
        #warning-handler-code
      }, error = function(e) {
        #error-handler-code
      }, finally = {
        #cleanup-code
      })
      
      
      tryCatch({print("performing drug-ADR SNP protein pathway lookup for related genes")
      # SNPS: get genes from pathways in proteinSNPS.txt in folder
      snps.fn <- paste(folderPath, "/proteinSNPS.txt", sep="")
      SNPSproteins.fn <- paste(folderPath, "/SNPSproteins.txt", sep="")
      if (file.exists(snps.fn) && file.info(snps.fn)$size > 0) {
        pidsSNP.fn <- paste(folderPath, "/pidsSNPS.txt", sep="")
        system(paste("./sedGetProteins.sh ", snps.fn, " >  ", SNPSproteins.fn, sep=""))
        system(paste("./sedGetPathways.sh ", snps.fn, " >  ", pidsSNP.fn,  sep=""))
        pids <- readVec(pidsSNP.fn)
        for (pid in pids) {
          print(pid)
          #getJeans(adr, ci, med, description, pid)
        }
      }
      }, warning = function(w) {
        #warning-handler-code
      }, error = function(e) {
        #error-handler-code
      }, finally = {
        #cleanup-code
      })
    }
  }
}


