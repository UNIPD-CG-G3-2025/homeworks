############################################################
# Homework 1 - SNP Data Analysis
# Course: Computational Genomics 2025/2026
# Date: 26/10/2025
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
# Filippo Tiberio, 
# Maria Camilla Rigo,
# Soukaina Elboulrhaiti, 
# Lucas Cappelletti, 
# Stella Rinaldi

############################################################

# LOADING THE DATA BANANA

#setwd("C:/Users/stell/OneDrive/Documenti/magistrale/CG")
#bc <- read.table("SNPdata.dat", stringsAsFactors = T)

#current_working_dir <- getwd()
#setwd(current_working_dir)
#bc <- read.csv("SNPdata.txt", sep = '\t')

# FIRST FUNCTION

compute_MAF <- function(SNPdata) {
  
  # Counting how many individuals for each phenotype, given:
  # AA -> 0, Aa -> 1, aa -> 2; removing NAs
  nAA_c <- rowSums(SNPdata == 0, na.rm = TRUE)
  nAa_c <- rowSums(SNPdata == 1, na.rm = TRUE)
  naa_c <- rowSums(SNPdata == 2, na.rm = TRUE)
  
  # Computing the total number of individuals in the sample
  N_c <- nAA_c + nAa_c + naa_c
  
  # Computing the Minor Allele Frequency avoiding the 0/0 division
  q <- (2*naa_c + nAa_c ) / pmax(2*N_c, 1)
  
  # Printing on terminal the number of SNPs (nrows) and subjects (ncols)
  nSNPs <- dim(SNPdata)[1]
  nSubj <- dim(SNPdata)[2]
  
  print(paste('The number of SNPs is:', nSNPs))
  print(paste('The number of Subjects is:', nSubj))
  
  return(q)
}



# PROVA CALCOLO BANANA
#q<-compute_MAF(bc[,])




# SECOND FUNCTION

HWE_test <- function(SNPdata, pdf_file_name) {
  # Calculate q and p
  q <- compute_MAF(SNPdata)
  p <- 1 - q
  
  # Count the alleles in patients
  nAA_p <- rowSums(SNPdata == 0, na.rm = TRUE)
  nAa_p <- rowSums(SNPdata == 1, na.rm = TRUE)
  naa_p <- rowSums(SNPdata == 2, na.rm = TRUE)
  
  # Count the total alleles in patients
  N_p <- nAA_p + nAa_p + naa_p
  
  # Calculate Chi-Square
  exp_AA <- N_p * p^2
  exp_Aa <- N_p * 2 * p * q
  exp_aa <- N_p * q^2
  
  chisq <- ((nAA_p - exp_AA)^2 / exp_AA) +
    ((nAa_p - exp_Aa)^2 /exp_Aa) + 
    ((naa_p - exp_aa)^2 / exp_aa)
  
  # Computing the relative p-value
  p_val <-pchisq(q = chisq, df = 1, lower.tail = F)
  names(p_val) <- rownames(SNPdata)
  
  
  #Find and print on terminal the name(s) of the SNP(s) having the lowest p-value
  minp <- min(p_val, na.rm = TRUE)
  snp_min <- names(p_val)[p_val == minp]
  cat("Smallest p-value SNP: [", snp_min, "] -> with p-value = [", sprintf("%.3e", minp), "]\n")
  
  
  # Boxplot of the computed values
  #pdf(file = pdf_file_name)
  boxplot(p_val, main = "HWE - PValues Boxplot", ylab = 'p-value', xlab = "SNPs")
  #dev.off()
  
  # Returning p-vals
  return(p_val)
}

#HWE<-HWE_test(bc[,])


# THIRD FUNCTION

SNP_association_test <- function(filepath, indCTRL, MAFth=0.01, HWEalpha=0.01) {
  MIN_THR_COUNT <- 5
  
  SNPdata = as.matrix(read.table(filepath, 
                                 header = TRUE, 
                                 sep = '\t', 
                                 row.names = 1))
  
  # Data preprocessing - removing SNPs and patients that have more than 10% NAs
  na<-is.na(SNPdata)
  # Counting the number of NAs in rows - checking if they are less than 10%
  na_sub<-colSums(na)/nrow(na)<=0.1
  # Counting the number of NAs in columns - checking if they are less than 10%
  na_snp<-rowSums(na)/ncol(na)<=0.1
  # Filtering step
  SNPdata_filt<-SNPdata[na_snp,na_sub]
  
  ctrl <- SNPdata_filt[, indCTRL, drop = FALSE]
  #nAA_ctrl <- rowSums(ctrl == 0, na.rm = TRUE)
  #nAa_ctrl <- rowSums(ctrl == 1, na.rm = TRUE)
  #naa_ctrl <- rowSums(ctrl == 2, na.rm = TRUE)
  
  #cond_ctrl <- nAA_ctrl < MIN_THR_COUNT | nAa_ctrl < MIN_THR_COUNT | naa_ctrl < MIN_THR_COUNT
  #SNPdata_filt <- SNPdata_filt[cond_ctrl==F, , drop = FALSE]
  
  
  ps <- SNPdata_filt[, -indCTRL, drop = FALSE]
  #nAA_ps <- rowSums(ctrl == 0, na.rm = TRUE)
  #nAa_ps <- rowSums(ctrl == 1, na.rm = TRUE)
  #naa_ps <- rowSums(ctrl == 2, na.rm = TRUE)
  
  #cond_ps <- nAA_ps < MIN_THR_COUNT | nAa_ps < MIN_THR_COUNT | naa_ps < MIN_THR_COUNT
  #SNPdata_filt <- SNPdata_filt[cond_ps==F, , drop = FALSE]
  
  # Computing the Minor Allele Frequency and testing for HWE on the controls
  MAF <- compute_MAF(ctrl)
  OK_maf <- MAF < MAFth
  
  HWE <- HWE_test(ctrl, pdf_file_name = 'test_boxplot.pdf')
  
  # Filtering for SNPs that have a MAF smaller than the threshold
  # and are not in HWE in the control population
  
  OK_hwe <- HWE < HWEalpha
  
  OK <- OK_maf | OK_hwe
  
  SNPdata_filt<-SNPdata_filt[OK==F, , drop = FALSE]
  
  
  # Contingency table for controls
  AA_c <- rowSums(SNPdata_filt[,indCTRL] == 0, na.rm = TRUE)
  Aa_c <- rowSums(SNPdata_filt[,indCTRL] == 1, na.rm = TRUE)
  aa_c <- rowSums(SNPdata_filt[,indCTRL] == 2, na.rm = TRUE)
  
  N_c <- AA_c + Aa_c + aa_c
  
  SNPs_c <- data.frame(AA_ctrl = AA_c,
                       Aa_ctrl = Aa_c,
                       aa_ctrl = aa_c,
                       row.names = rownames(SNPdata_filt))
  
  # Contingency table for patients
  AA_p <- rowSums(SNPdata_filt[,-indCTRL] == 0, na.rm = TRUE)
  Aa_p <- rowSums(SNPdata_filt[,-indCTRL] == 1, na.rm = TRUE)
  aa_p <- rowSums(SNPdata_filt[,-indCTRL] == 2, na.rm = TRUE)
  
  N_p <- AA_p + Aa_p + aa_p
  
  SNPs_p <- data.frame(AA_case = AA_p,
                       Aa_case = Aa_p,
                       aa_case = aa_p,
                       row.names = rownames(SNPdata_filt))
  N_tot<-N_p+N_c
  
  # Expected occurances:
  # totality of patients for AA, Aa and aa cases separately
  EX<-SNPs_p+SNPs_c
  # expectation for patients
  EX_p<-((EX/N_tot)*rowSums(SNPs_p))
  # expectation for controls
  EX_c<-((EX/N_tot)*rowSums(SNPs_c))

  cond_c <- rowSums(EX_c < 5) == 0
  cond_p <- rowSums(EX_p < 5) == 0
  cond_all <- cond_p & cond_c
  
  # Computing the Chi^2, handling NAs
  c1<-(((SNPs_p-EX_p)^2)/EX_p)
  c1[is.na(c1)]<-0
  c2<-(((SNPs_c-EX_c)^2)/EX_c)
  c2[is.na(c2)]<-0
  Chi<-rowSums(c1)+rowSums(c2)
  pval <-pchisq(q = Chi, df = 2, lower.tail = F)
  
  # ragionare BANANA
  #pval<- pval[cond_all==T, , drop = FALSE]
  mat<-cbind(SNPs_c,SNPs_p,pval)
  mat<- mat[cond_all==T, , drop = FALSE] #Filtering rows with count lower than THR
  
  mat$qval<-mat$pval*(length(mat$pval))/rank(mat$pval)
  
  #colnames(mat)<-c("AA_ctrl","Aa_ctrl","aa_ctrl","AA_case","Aa_case","aa_case","pval","qval")
  
  return(mat)
}

# togliere le prove di calcolo alla fine BANANA
SNP_ass_test <- SNP_association_test("SNPdata.txt", 1201:2000)




