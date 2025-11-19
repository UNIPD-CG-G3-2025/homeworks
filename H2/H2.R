############################################################
# Homework 2 - RNAseq data analysis
# Course: Computational Genomics 2025/2026
# Date: 19/11/2025
#
# Description:
# This script implements three functions for RNAseq data analysis:
# 1. ...
# 2. ...
# 3. ...
#
# Group members:
# Lucas Cappelletti, 
# Soukaina Elboulrhaiti, 
# Maria Camilla Rigo,
# Stella Rinaldi,
# Filippo Tiberio
############################################################


# FIRST FUNCTION
MvAplot <- function(exprData, pdffilename, pcolor="black", lcolor="red") {
  reference <- exprData[1]
  reference_name <- colnames(exprData)[1]
  samples <- (length(exprData))
  
  pdf(file = pdffilename)
  for(i in 2:samples){
    current <- exprData[i]
    
    M<-(log2(reference)-log2(current))[[1]]
    A<-((log2(reference)+log2(current))/2)[[1]]
    
    max_max <- max(M[!is.infinite(M)], na.rm=T)
    max_min <- max(abs(M[!is.infinite(M)]), na.rm=T)
    
    ylim <- max(c(max_max, max_min), na.rm=T)
    
    plot(x = A, y = M,
         pch = ".", col = pcolor,
         ylim=c(-ylim, ylim),
         main=paste("MvA plot of sample 1 (", reference_name, ") vs sample ",i, " (", colnames(exprData)[i], ")"),
         xlab="M",
         ylab="A",
    )
    abline(h = 0, col = lcolor)
  }
  
  dev.off()
}

# SECOND FUNCTION
TMMnorm <- function(exprData, annot, Mtrim = 0.02, Atrim = c(0,8)) {
  seq_depths <- colSums(exprData)
  exprData <- exprData / seq_depths * (10^6) # Va fatto anche sul primo???
  
  lengths <- annot[rownames(exprData), "Length"] #Estraggo le lunghezze, impostando a NA quelli missing
  median_length <- median(annot$Length, na.rm = TRUE) #Estraggo la mediana
  lengths[is.na(lengths)] <- median_length #La imposto come valore per i missing values
  
  
  normalized <- log2(exprData) # mi porto in log per uniformare il calcolo
  reference <- normalized[1]
  samples <- (length(exprData))
  
  sfs <- rep(0, samples)
  
  for(i in 2:samples){
    current <- normalized[i]
    
    M <- (reference - current)[[1]]
    A <- ((reference + current)/2)[[1]]
    
    cond <- (A > Atrim[1] & A < Atrim[2])
    M_filtered <- M[cond]
    M_sorted <- sort(M_filtered)
    
    sfs[i] <- mean(M_sorted, trim = Mtrim, na.rm = T)
    
    normalized[,i] <- normalized[,i] + sfs[i]
  }
  
  # Qui va iniziato lo scaling con la lunghezza e fatto x 10^3
  normalized <- (2^normalized) * lengths / (10^3)  #In che scala??
  
  return (list(2^sfs, normalized)) #ritorno riportando in lineare
}

DATA <- read.table("raw_count.txt", sep="\t", row.names=1, header=TRUE)
annotations <- read.table("gene_annot.txt", sep="\t", row.names=2, header=TRUE, quote = "\"")

MvAplot(DATA, "MvA_Plot.pdf")
normalized <- TMMnorm(DATA, annotations)
MvAplot(normalized[[2]], "TEST_normalized.pdf")

