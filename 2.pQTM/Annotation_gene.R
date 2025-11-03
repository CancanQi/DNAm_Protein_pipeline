### use biomart to get the location information of gene that encodes the protein

library(biomaRt)
setwd("~/Documents/Projects/airway_eqtl/PIAMA/")
pr.matrix<-read.table("...")

gene<-rownames(pr.matrix) ### if the rownames are names of the protein
mart <- useEnsembl("ENSEMBL_MART_ENSEMBL",GRCh = 38) ### need to use GRCh 38, consistent with CpG annotation
ensemble<-useDataset(dataset="hsapiens_gene_ensembl",mart = mart)

### get position information of the gene that encodes the protein
anno_pr<-getBM(attributes = c('ensembl_gene_id','hgnc_symbol', 'chromosome_name',
                               'start_position', 'end_position', 'band', "description"),
                filters = 'hgnc_symbol', 
                values = gene, 
                mart = ensemble)

### change the data format required for matrixeqtl
gene.loc<-anno_rna[,c("hgnc_symbol","chromosome_name","start_position","end_position")]
colnames(gene.loc)<-c("geneid","chr","s1","s2")
write.table(gene.loc,file = "./data/input_matrixeqtl/geneloc.txt",sep = "\t",quote = F,row.names = F)
