#working Jan 30 2020 - try not to touch

#setwd("/home/jm33a/domain_evolution/1KT_searches/CLV1/") 
#run <- "CLV1" 
.libPaths("/share/pkg/R/PACKAGES/3.6.1/tidyverse/1.3.0") #this is needed here to prevent R from seeking tidyverse in personal directory
library(tidyverse)


run <- commandArgs(TRUE) #the name of the run should be passed from shell script
#run<-"CLV1"


#debugging
#both_domains_genes <- matrix(data = NA, nrow = 2, ncol = 2)
#write.csv(both_domains_genes, file = paste(run, "test.csv", sep = "."))
#quit(save = "no")


#import the table outputted by pfam domain scans
#raw_pfam_df <- data.frame(read_tblout(paste(run,"pfamout.tsv", sep = "."))) #old version, uses rhmmer read_tblout which I can't get to run on cluster
raw_pfam_df <- read.table(paste(run,"pfamout.tsv", sep = "."), fill=TRUE)[,c(1,3)] #import, just the gene name and domain call. empty spaces are annoying so use fill=true to allow, the extra rows are not a problem since they don't have Lrr or kinase and will be screened out later
  colnames(raw_pfam_df) <- c("domain_name", "query_name")
  raw_pfam_df$query_name <- as.character(raw_pfam_df$query_name) #convert to character strings so they won't mess with the later steps

#summarize all domain finds by gene
pfam_df <- raw_pfam_df %>% 
  select(domain_name, query_name) %>%
  group_by(query_name)  %>% 
  dplyr::summarise("domains_found" = paste(domain_name, collapse=", "))  #summarize() is shared by several packages. dplyer::summarize() is necessary to make this work

#make new table of genes with both domains found and print to a file
both_domains_genes <- dplyr::filter(pfam_df, grepl('Pkinase', domains_found) & grepl('LRR', domains_found))[,1] 
write.table(both_domains_genes, row.names = F, col.names = F, quote = F, file = paste(run,"hits_with_both_domains_IDs.txt", sep = ".")) 

quit(save = "no")
