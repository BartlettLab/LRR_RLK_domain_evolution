install.packages("zoo")
library(zoo)

install.packages("cowplot")
library(cowplot)
BiocManager::install("Biostrings")
library(Biostrings)

BiocManager::install("msa")
library(msa)

library(tidyverse)

setwd("/home/jm33a/domain_evolution/1KT_searches/CLV1/tree_from_hits/follow_up_tree/")
setwd("/home/jm33a/domain_evolution/1KT_searches/HAE/site_rate_evo/")


setwd("/home/jm33a/mTERF/shot1_and_sister_branch_1KT_search/round2_from_node_3926/round3_with_genome_seqs/round4_round3_without_weirdos/infer_rate/")


#hae residues in contact with ligand  
contact_residues <- data.frame(c(409,407,385,361,364,364,383,339,337,337,315,313,290,290,268,268,313,290,266,266,264,264,240,242,218,196,196,196,196,172,148))
colnames(contact_residues) <- "Site"
df <- inner_join(contact_residues, tab4, by = "Site")
contact_residues <- df[df$metric == "Rate",c(1,3)]
typeof(contact_residues$evolution_rate)
contact_residues$evolution_rate <- as.numeric(contact_residues$evolution_rate)
contact_residues$Site <- as.numeric(contact_residues$Site)


tab1 = read.table('align.fasta.rate',header=TRUE)[, c("Site", "Rate")] #read in the site rates for the alignment

#remove alignment columns that do not have sequence in arabidopsisreadAAMultipleAlignment
aa1 <-  readAAMultipleAlignment(filepath = "input_align.fasta", format = "fasta") #bring in the entire MSA
aa2 <- as(aa1, "BStringSet") #convert to msa object


####This can be used for several or just one gene#
for (target_gene in c(" AT1G75820.1_CLV1", " AT1G73080.1_PEPR1")){
  #    target_gene <-  "AT4G28490.1" #HAESA
  #    target_gene <-  "AT3G60400.1_SHOT1" #SHOT1
 
  target_seq <- as.character(aa2[target_gene]) #get the sequence of gene in alignment, including dashes
   tab1$target_seq <- unlist(strsplit(target_seq, "")) #print these to a new column
   tab2 <- tab1[tab1$target_seq != "-",] #tab2 is the site rates for each position from gene but missing all the gaps
  joined_table <- full_join(tab1, tab2, by = "Site")
  tab1[, target_gene] <- joined_table$Rate.y

   names(tab2)[names(tab2) == "Site"] <- "orig_site"
   tab2 <- tab2 %>% mutate(Site = row_number())
   
sliding_window_size <- 50
   
      tab3 <- tab2 %>% #tab 3 adds sliding window mean
     select(Site, Rate) %>%
     mutate(sliding_window = rollmean(Rate, k = sliding_window_size, fill = NA))#,
       sliding_window_2 = rollmean(Rate, k = 81, fill = NA),
       sliding_window_3 = rollmean(Rate, k = 61, fill = NA),
       sliding_window_4 = rollmean(Rate, k = 41, fill = NA),
       sliding_window_5 = rollmean(Rate, k = 21, fill = NA),
       sliding_window_6 = rollmean(Rate, k = 11, fill = NA),
       sliding_window_7 = rollmean(Rate, k = 7, fill = NA))
 
  # collection_table$target_gene <- tab3$sliding_window_100 
}


tab4 <- tab3 %>%
  #gather(metric, evolution_rate, c(Rate, sliding_window_1, sliding_window_2, sliding_window_3, sliding_window_4, sliding_window_5, sliding_window_6)) 
  gather(metric, evolution_rate, c(Rate, sliding_window))


lrr_range <- 90:593
tm_range <- 622:641
rlk_range <- 683:968
  
ggplot() +
    #geom_rect(data=contact_residues, mapping=aes(xmin= Site -2, xmax=Site + 2, ymin=.5, ymax=1), fill = "orange", alpha = 1) +
   geom_line(data = tab4, aes(Site, evolution_rate, color = metric, alpha=metric)) + 
     xlim(100,400) +
  scale_color_manual(values=c('gray','forestgreen')) + 
  scale_alpha_manual(values = c(.4,1),guide=F) +
  theme_bw() +
  #geom_point(data = contact_residues, aes(Site, evolution_rate), color = "red") + 
  #geom_segment(aes(x=min(lrr_range),xend=max(lrr_range),y=mean(na.omit(tab2[lrr_range,'Rate'])),yend=mean(na.omit(tab2[lrr_range,'Rate']))),linetype = "dotted") + #LRR average
    annotate("text", hjust = 1, x=min(lrr_range)-10, y=mean(na.omit(tab2[lrr_range,'Rate'])), label = "LRR mean", size = 2) +
  geom_segment(aes(x=min(tm_range),xend=max(tm_range),y=mean(na.omit(tab2[tm_range,'Rate'])),yend=mean(na.omit(tab2[tm_range,'Rate']))),linetype = "dotted") + #TM average
   annotate("text", hjust = 1, x=min(tm_range)-10, y=mean(na.omit(tab2[tm_range,'Rate'])), label = "TM mean", size = 2) +
    geom_segment(aes(x=min(rlk_range),xend=max(rlk_range),y=mean(na.omit(tab2[rlk_range,'Rate'])),yend=mean(na.omit(tab2[rlk_range,'Rate']))),linetype = "dotted") + #RLK average
      annotate("text", hjust = 0, x=max(rlk_range), y=mean(na.omit(tab2[rlk_range,'Rate'])), label = "RLK mean", size = 2) +
  theme(legend.position = c(.85, .85), legend.background = element_rect(fill="white"), legend.box.background = element_rect(colour = "black", size = 2))
  #this geom_ribbon fills in area between line and average. But still needs some work...
  # geom_ribbon(data = tab4, mapping = aes(x = Site, y = evolution_rate, ymin = mean(tab2$Rate), ymax = evolution_rate, fill = metric), linetype = "dotted") 



ggplot() + geom_point(data = contact_residues, aes(Site, evolution_rate)) 





#hae_seq <- as.character(aa2$AT4G28490.1) #print sequence of HAE in alignment, including dashes
#hae_seq <- as.character(aa2$AT1G28440.1) #print sequence of HSL1 in alignment, including dashes
#hae_seq <- as.character(aa2$Solyc03g006300.1.1) #print sequence of SlHAE in alignment, including dashes
hae_seq <- as.character(aa2$LOC_Os01g13800.1) #print sequence of OsHAE in alignment, including dashes

