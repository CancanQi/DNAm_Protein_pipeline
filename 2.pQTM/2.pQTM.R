### run pQTM using MatrixEQTL
### reference: https://github.com/andreyshabalin/MatrixEQTL

library(MatrixEQTL)

setwd("...")
useModel = modelLINEAR; ### set model to use
SNP_file_name="./data/input_matrixeqtl/snp_for_eqtl.txt";  ### here we input the data of CpG sites, same format as the SNP file, see demo of this file
snps_location_file_name = "./data/input_matrixeqtl/snpsloc.txt";  ### location information of the CpG sites

expression_file_name="./data/input_matrixeqtl/gene_for_eqtl.txt"; ### here we input the data of proteins, same format as the gene expression file, see demo 
gene_location_file_name = "./data/input_matrixeqtl/geneloc.txt";  ### location information of the gene that encode the protein

covariates_file_name="./data/input_matrixeqtl/cov_for_eqtl.txt"  ### covariates used for pQTM analysis, usually, age, sex, batch effect

output_file_name_cis = "./results/Res_all_cis.txt";  ### location and file name for cis-pQTM results
output_file_name_tra = "./results/Res_all_trans.txt";  ### location and file name for trans-pQTM results

pvOutputThreshold_cis = 1  ### p value threshold for output, set as 1 if you want all the results
pvOutputThreshold_tra = 1e-6

cisDist = 1000000 ### distance to define which is "cis" region, usually 1mb

errorCovariance = numeric();

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;
snps$fileSkipColumns = 1;
snps$fileSliceSize = 2000;
snps$LoadFile(SNP_file_name);

gene = SlicedData$new();
gene$fileDelimiter = "\t";
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;
gene$fileSkipColumns = 1;
gene$fileSliceSize = 2000;
gene$LoadFile(expression_file_name);

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$LoadFile(covariates_file_name);

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis=output_file_name_cis,
  pvOutputThreshold.cis=pvOutputThreshold_cis,
  snpspos=snpspos,
  genepos=genepos,
  cisDist=cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

## ================================================
## extract HLA snp results re-calculate FDR
## if you need to pay attention to 

hla<-read.table("data/input_matrixeqtl/snpsloc_hla.txt",header = T)

res.hla<-res.cis[res.cis$snps %in% hla$snp,]
res.hla$FDR_new<-p.adjust(res.hla$pvalue,method = "fdr")
dim(res.hla) # 616

res.non.hla<-res.cis[!(res.cis$snps %in% hla$snp),]
res.non.hla$FDR_new<-p.adjust(res.non.hla$pvalue,method = "fdr")
dim(res.non.hla) # 1634

res.hla.sig<-res.hla[which(res.hla$FDR_new<0.05),]
res.non.hla.sig<-res.non.hla[which(res.non.hla$FDR_new<0.05),]

write.csv2(res.hla.sig,file = "Res_cis_sig_hla.csv")
write.csv2(res.non.hla.sig,file = "Res_cis_sig_non_hla.csv")

