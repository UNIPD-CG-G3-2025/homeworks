############################################################
# Homework 1 - SNP Data Analysis
# Course: Computational Genomics 2025/2026
# Date: 20/10/2025
#
# Description:
# This script implements three functions for SNP analysis:
# 1. compute_MAF(SNPdata): computes Minor Allele Frequencies (MAF)
# 2. HWE_test(SNPdata, pdf_file_name): performs Hardy-Weinberg Equilibrium test
# 3. SNP_association_test(filepath, indCTRL, MAFth=0.01, HWEalpha=0.01):
#    filters SNPs by MAF and HWE, performs association tests, and
#    computes Benjamini–Hochberg corrected q-values.
#
# Group members:

############################################################

compute_MAF <- function(SNPdata) {
  nAA_c <- rowSums(SNPdata == 0, na.rm = TRUE)
  nAa_c <- rowSums(SNPdata == 1, na.rm = TRUE)
  naa_c <- rowSums(SNPdata == 2, na.rm = TRUE)
  
  N_c <- nAA_c + nAa_c + naa_c
  
  q <- (2*naa_c + nAa_c ) / pmax(2*N_c, 1)
  
  return(q)
}

HWE_test <- function(SNPdata, pdf_file_name, q) {
  
  # Count the alleles in patients
  nAA_p <- rowSums(SNPdata == 0, na.rm = TRUE)
  nAa_p <- rowSums(SNPdata == 1, na.rm = TRUE)
  naa_p <- rowSums(SNPdata == 2, na.rm = TRUE)
  
  # Count the total alleles in patients
  N_p <- nAA_p + nAa_p + naa_p
  
  # Calculate p
  p <- 1 - q
  
  
  # Calculate Chi-Square
  chisq <- ((nAA_p - (N_p*p^2))^2) / (N_p*p^2) 
            + ((nAa_p - (N_p*p*q*2))^2) / (N_p*p*q*2) 
            + ((naa_p - (N_p*q^2))^2) / (N_p*q^2)
  
  
  # And relative p-value
  p_val <- pchisq(q = chisq, df = 1)
  names(p_val) <- rownames(SNPdata)
  
  
  #Find and print the terminal the name(s) of the SNP(s) having the lowest p-value
  minp <- min(p_val, na.rm = TRUE)
  snp_min <- names(p_val)[p_val == minp]
  cat("Smallest p-value SNP: ", snp_min, " -> with p-value = ", sprintf("%.3e", minp), "\n")
  
  
  # Boxplot of the computed values
  pdf(file = pdf_file_name)
  boxplot(p_val, main = "HWE - PValues Boxplot", ylab = 'p-value', xlab = "SNPs")
  dev.off()
  
  # Returning p-vals
  return(p_val)
}


SNP_association_test <- function(filepath, indCTRL, MAFth=0.01, HWEalpha=0.01) {
  #BODY OF THE FUNCTION
  
  SNPdata = as.matrix(read.table(filepath, header = TRUE, sep = '\t', row.names = 1))
  
  ctrl <- SNPdata[, indCTRL, drop = FALSE]
  ps <- SNPdata[, -indCTRL, drop = FALSE]
  
  q <- compute_MAF(ctrl)
  
  
  HWE_test(ps, pdf_file_name = 'test_boxplot.pdf', q)
  
  return(0)
}

SNP_association_test("SNPdata.txt", 1201:2000)





