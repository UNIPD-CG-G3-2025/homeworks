############################################################
# Homework 2 - RNAseq data analysis
# Course: Computational Genomics 2025/2026
# Date: 22/11/2025
#
# Description:
# This script implements three functions for RNAseq data analysis:
# 1. MvA plot
# 2. TMMnorm
# 3. DEbyEdgeR
#
# Group members:
# Lucas Cappelletti, 
# Soukaina Elboulrhaiti, 
# Maria Camilla Rigo,
# Stella Rinaldi,
# Filippo Tiberio
############################################################

# First Function
MvAplot <- function(exprData, pdffilename, pcolor="black", lcolor="red") {
  reference <- exprData[1]
  reference_name <- colnames(exprData[1])
  samples <- (length(exprData))
  
  pdf(file = pdffilename)
  for(i in 2:samples){
    current <- exprData[,i]
    
    M<-(log2(reference)-log2(current))[[1]]
    A<-((log2(reference)+log2(current))/2)[[1]]
    
    max_max <- max(M[!is.infinite(M)], na.rm=T)
    max_min <- max(abs(M[!is.infinite(M)]), na.rm=T)
    
    ylim <- max(c(max_max, max_min), na.rm=T)
    
    plot(x = A, y = M,
         pch = ".", col = pcolor,
         ylim=c(-ylim, ylim),
         main=paste("MvA plot of sample 1 (", reference_name, ") vs sample ",i, " (", colnames(exprData)[i], ")"),
         xlab="A",
         ylab="M",
    )
    abline(h = 0, col = lcolor)
  }
  
  dev.off()
}

# Second Function
TMMnorm <- function(exprData, annot, Mtrim=0.02, Atrim = c(0,8)) {
  #scaling by the sequencing depth and multiplying by 10^6
  
  SD<-colSums(exprData[,1:ncol(exprData)], na.rm = T)  # sequencing depth
  data<-sweep(exprData[,1:ncol(exprData)],2,SD, FUN='/')*10^6 
  
  
  normData<-data   
  
  #SCALING FACTORS
  ni<-ncol(data)
  # initialization of SF vector
  SF<-rep(1,times=ni)

  for (i in 2:ni){
    # computing A and M with an offset
    offset=0.0001
    M<-log2(data[,1]+offset)-log2(data[,i]+offset)
    A<-(log2(data[,1]+offset)+log2(data[,i]+offset))/2
    # computing the trimmed mean
    indA<-A>Atrim[1]&A<Atrim[2]
    SF[i]<-mean(M[indA], trim=Mtrim, na.rm=T)
    SF[i]<-2^(SF[i])
    normData[,i]<-normData[,i]*SF[i]
  }
  
  # Scaling for the length 
  annot <- annot[match(rownames(exprData), rownames(annot)), ] #check row names 
  median_length <- median(annot$Length, na.rm = TRUE)
  len<-annot$Length
  len[is.na(len)] <- median_length
  normData <- sweep((normData),1,len, FUN="/")*(10^3)
  
  # returning the list
  
  l<-list(normData= normData,SF = SF)
  return(l)
}

# Third Function
DEbyEdgeR <- function(rawdata, groups, alpha = 0.05) {
  # N.B.: REQUIRES THE "edgeR" LIBRARY
  library(edgeR)
  
  # Finding the indices of columns that contain 'CTRL' (case-insensitive)
  indices_controls <- grep("CTRL", groups, ignore.case = TRUE)
  
  # Creating the group factor
  label <- factor(ifelse(seq_along(groups) %in% indices_controls,
                         "Control", "Case"))
  
  
  # Creating the DGEList object and computing normalization factors
  y <- DGEList(counts = rawdata)
  y <- calcNormFactors(y)
  
  # Building the design matrix
  design <- model.matrix(~ label) 
  rownames(design) <- colnames(rawdata) 
  
  # Fit values of phi and fitting the model, taking the results
  y <- estimateGLMCommonDisp(y, design, verbose=TRUE) 
  y <- estimateGLMTrendedDisp(y, design) 
  y <- estimateGLMTagwiseDisp(y, design) 
  fit <- glmFit(y, design)
  RES <- glmLRT(fit)
  out <- topTags(RES, n = "Inf")$table
  
  # G: total number of tests run
  G <- nrow(rawdata)
  
  # Computing the QValues
  out$QValue <- out$PValue*G/rank(out$PValue, ties.method = 'max')
  
  # Choosing the differentially expressed genes evaluating if their PValue < alpha
  DEGs <- out[out$PValue < alpha, c("PValue", "QValue", "logFC")]
  
  # nDEGs: selected genes
  # Computing False Positives, True Positives, True Negatives, False Negatives and FDR
  
  nDEGs <- sum(out$PValue < alpha)
  G0 <- 0.8*G
  FPs <- min(nDEGs, G0*alpha)
  TPs <- max(0, nDEGs - FPs)
  TNs <- min(G-nDEGs, G0 - FPs)
  FNs <- G - FPs - TPs - TNs
  FDR <- FPs/nDEGs
  
  # Returning results
  vec <- c("Selected_genes" = nDEGs,
           "FP" = FPs,
           "FN" = FNs,
           "FDR" = FDR)
  
  mat <- as.matrix(DEGs)
  results <- list(vec, mat)
  return(results)
}
