#Relative abundance - co-assembly
library(readxl)
library(tidyverse)
library(ggplot2)
library(viridis)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(tibble)
library(plyr)
library(dplyr)
library(hrbrthemes)
library(ggpubr)
library(patchwork)
library(cowplot)
library(ggsignif)
library(rstatix)
library(ggsci)
library(circlize)
library(ComplexHeatmap)
library(pheatmap)
library(gridExtra)

#import data
sid_trans <- read_excel("~/Desktop/BINGO2/MAGs/Fegenie/FeGenie-geneSummary-clusters.xlsx", 
                                               sheet = "Sid_trans_clean")

duplicated(sid_trans$MAGs)

####Stacked bar plot - glacial vs non-glacial ASV level#####
X1_bingo_sequence_table_SV_incubation <- read_excel("~/Desktop/BINGO1/ASV_bingo-result/1.bingo-sequence_table_SV_incubation.xlsx", 
                                                       sheet = "Sheet4")
bingo1_asv <- X1_bingo_sequence_table_SV_incubation

bingo1_asv1<- aggregate( . ~ Family, bingo1_asv, sum)
bingo1_asv2 <- melt(setDT(bingo1_asv1), id.vars = c("Family"), variable.name = "Replicates")
bingo1_asv3<-bingo1_asv2[!(bingo1_asv2$Family=="NA"),]

colourCount = length(unique(bingo1_asv3$Family))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

ggplot(bingo1_asv3, aes(fill=Family, y=value, x=Replicates)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = getPalette(colourCount)) +
labs(y = "Relative Abundance (Family > 1%)") 


display.brewer.all()

ggplot(bingo1_asv3) +
  aes(x = Replicates, y = value, fill = Family) +
  geom_bar(position="fill", stat="identity") +
  geom_col() +
  scale_fill_hue(direction = 1) +
  labs(y = "Relative Abundance (Family > 1%)") 

# install.packages("Polychrome")
library(Polychrome)

# build-in color palette
Glasbey = glasbey.colors(32)
swatch(Glasbey)



#########MAGs siderophore synthesis#####
Cov_MAG <- read_excel("anvio_bins_dastool/bins_across_samples/mean_coverage.xlsx", 
                           sheet = "Data_sid_trans")

Cov_MAG$Group <- as.factor(Cov_MAG$Group)
levels(Cov_MAG$Group)

###### Stacked bar plot
ggplot(Cov_MAG, aes(fill=MAG, y=Coverage, x=factor(Group, level=c('LacAMP', 'LacADN')),)) + 
  geom_bar(position="stack", stat="identity") + 
  geom_col() +
  xlab('Group')
dev.off()

ggplot(Cov_MAG) +
  aes(fill=MAG, y=Coverage, x=factor(Group, level=c('LacAMP', 'LacADN')),) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_brewer(palette="Dark2") +
  labs(y = "Mean coverage of MAGs") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 13L,
                                face = "bold"),
    axis.title.x = element_blank())

### Grouped bar plot
ggplot(Cov_MAG, aes(fill=MAG, y=Coverage, x=factor(Group, level=c('LacAMP', 'LacADN')),)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab('Group')
dev.off()

#########MAGs siderophore transport #####
trans_MAG <- read_excel("antismash/Siderophore_transport.xlsx", 
                            sheet = "Sheet2")
trans_MAG$Group <- as.factor(trans_MAG$Group)
levels(trans_MAG$Group)


###### Stacked bar plot
#tiff("MAG_sid_trans.tiff", units="in", width=5, height=5, res=300)
ggplot(trans_MAG) +
  aes(fill=Family, y=Coverage, x=factor(Group, level=c('LacAMP', 'LacADN')),) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_brewer(palette="Dark2") +
  labs(y = "Mean coverage of MAGs") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 13L,
                                face = "bold"),
    axis.title.x = element_blank()) +
  #geom_text(aes(label=MAG_quality),color="white",size=3,position=position_stack(vjust=0.5)) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
   ggtitle("MAGs containing siderophore transport genes") + 
  theme(plot.title = element_text(hjust = 0.5)) 

### MAGs containing only siderophore transport genes
trans_MAG_1 <- read_excel("antismash/Siderophore_transport.xlsx", 
                        sheet = "Clean")
trans_MAG_1$Group <- as.factor(trans_MAG_1$Group)
levels(trans_MAG_1$Group)
ggplot(trans_MAG_1) +
  aes(fill=Family, y=Coverage, x=factor(Group, level=c('LacAMP', 'LacADN')),) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_brewer(palette="Dark2") +
  labs(y = "Mean coverage of MAGs") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 13L,
                                face = "bold"),
    axis.title.x = element_blank()) +
  #geom_text(aes(label=MAG_quality),color="white",size=3,position=position_stack(vjust=0.5)) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
  ggtitle("MAGs containing only siderophore transport genes") + 
  theme(plot.title = element_text(hjust = 0.5)) 


#### Siderophore synthesis and siderophore transport in Nitrocolaceae and Psuedomonadaecaea MAG ####
Df_all <- read_excel("mean_coverage.xlsx", 
                     sheet = "Data_sid_trans")
Df_all$MAG <- as.factor(Df_all$MAG)
levels(Df_all$MAG)
Df_all$number <- 1

sum_df_all <- ddply(Df_all,.(MAG, HMM_siderophore_combination), #find the mean of all data, list all separate groups
                       summarize, count=sum(number)) 
ggplot(sum_df_all) +
  aes(x = MAG, y = count, fill = HMM_siderophore_combination) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_brewer(palette="Paired") +
  labs(y = "Gene count") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 13L,
                                face = "bold"),
    axis.title.x = element_blank())

ggplot(sum_df_all, aes(x=MAG, y=HMM_siderophore_combination, size = count)) +
  geom_point(alpha=0.7)

sum_df_all %>%
  arrange(desc(count)) %>%
  #mutate(MAG = factor(MAG_1, MAG_2)) %>%
  ggplot(aes(x=MAG, y=HMM_siderophore_combination, size=count, color=HMM_siderophore_combination)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(.1, 24), name="Gene count")


########## bidirectional plot - MAGs ####
all_mags <- read_excel("anvio_bins_dastool/bins_across_samples/mean_coverage.xlsx", 
                          sheet = "All")
all_mags$Group <- as.factor(all_mags$Group)
factor(all_mags$Group)
all_mags$Group <- factor(all_mags$Group, levels=c('LacAMP', 'LacADN'))
ggplot(all_mags,aes(x=MAG, y = ifelse(Group == "LacADN", Coverage, -Coverage) , fill=Group))+ 
  geom_bar(stat="identity", position="identity")+
  labs(y= "Mean Coverage",
       x = "MAGs")+
  scale_y_continuous(limits = c(-max(all_mags$Coverage), max(all_mags$Coverage))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15))+
  coord_flip()


#### only siderophore transport 
trans_MAG$Group <- factor(trans_MAG$Group, levels=c('LacAMP', 'LacADN'))
ggplot(trans_MAG,aes(x=Family, y = ifelse(Group == "LacADN", Coverage, -Coverage) , fill=Group))+ 
  geom_bar(stat="identity", position="identity") +
  labs(y= "Mean Coverage",
       x = "MAGs")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) + 
  scale_y_continuous(limits = c(-max(trans_MAG$Coverage), max(trans_MAG$Coverage))) +
  coord_flip()


##### Siderophores - New plots ####
##  both siderophore synthesis and transporters
both_mags <- read_excel("antismash/Replicates.xlsx", 
                            sheet = "Both ")
both_mags$Coverage <- as.numeric(both_mags$Coverage)

both_mags$Group <- as.factor(both_mags$Group)
both_mags$Group <- factor(both_mags$Group, levels=c('LacAMP', 'LacADN'))
summary(both_mags)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
  
both_mags_1 <- data_summary(both_mags, varname="Coverage", 
                      groupnames=c("Group", "Family"))
# Convert Group to a factor variable
both_mags_1$Group=as.factor(both_mags_1$Group)
head(both_mags_1)

### simple bar plot 
ggplot(both_mags_1) +
  aes(fill=Family, y=Coverage, x=factor(Group, level=c('LacAMP', 'LacADN')),) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_brewer(palette="Paired") +
  labs(y = "Mean coverage of MAGs") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 13L,
                                face = "bold"),
    axis.title.x = element_blank()) +
  #geom_text(aes(label=MAG_quality),color="white",size=3,position=position_stack(vjust=0.5)) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
  ggtitle("MAGs containing siderophore synthesis and transport genes") + 
  theme(plot.title = element_text(hjust = 0.5)) 

##Funnel plot - both 
plt_1 <- both_mags_1 %>%  # Cast the users table as a number
  mutate(Coverage = as.numeric(Coverage)) %>%
  # Pipe into ggplot
  ggplot(aes(x = reorder(Family,abs(Coverage)), y = ifelse(Group == "LacADN", Coverage, -Coverage), ymin = Coverage - sd, ymax = Coverage + sd,
             fill = Group)) +
  # Plot the bars
  geom_bar(stat="identity", position="identity") +
  # add colours 
  scale_fill_manual(values=c('#999999','#E69F00')) +
  # Shift the y axis
  scale_y_continuous(limits = c(-max(both_mags$Coverage), max(both_mags$Coverage))) +
  # Flip the coordinates
  coord_flip() +
  # Add error bars
  geom_errorbar(width = 0.1) +
  # Add a theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y= "Mean Coverage",
       x = "MAGs")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) 

  plt_1
  
  ## only transporters 
  trans_MAG <- read_excel("antismash/Siderophore.xlsx", 
                          sheet = "only_trans")
  trans_MAG$Coverage <- as.numeric(trans_MAG$Coverage)
  
  trans_MAG$Group <- as.factor(trans_MAG$Group)
  trans_MAG$Group <- factor(trans_MAG$Group, levels=c('LacAMP', 'LacADN'))


  ## Funnel plot - only transporters
plt_2 <- trans_MAG %>%  # Cast the users table as a number
  mutate(Coverage = as.numeric(Coverage)) %>%
  # Pipe into ggplot
  ggplot(aes(x = reorder(Family,abs(Coverage)), y = ifelse(Group == "LacADN", Coverage, -Coverage) , fill = Group)) +
  # Plot the bars
  geom_bar(stat="identity", position="identity") +
  # Shift the y axis
  scale_y_continuous(limits = c(-max(trans_MAG$Coverage), max(trans_MAG$Coverage))) +
  # Flip the coordinates
  coord_flip() +
  # Add a theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y= "Mean Coverage",
       x = "MAGs")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) 

plt_2

## Box plots with error bars facetted by MAG ####
## for both siderophore producers and transporters
mags_both  <- read_excel("antismash/Replicates.xlsx", 
                             sheet = "Both ")
mags_both$Group <- as.factor(mags_both$Group)
mags_both$Group <- factor(mags_both$Group, levels=c('LacAMP', 'LacADN'))
mags_both$Coverage <- as.numeric(mags_both$Coverage)
mags_both_1 <- mags_both_rep %>%
  mutate(Family = fct_reorder(Family, desc(Coverage))) %>%
  ggplot(aes(x=Family, y=Coverage, fill=Group)) + 
  geom_boxplot() +
  facet_wrap(~Family, scale="free") +
  #stat_compare_means( label =  "p.signif", label.x = 1.5) + 
  #guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
   scale_fill_brewer(palette="Paired") +
  # change order of legend
  # Add a theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y= "Mean Coverage of MAGs",
       x = "MAGs")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) 

mags_both_1
## by replicates
mags_both_rep <- read_excel("antismash/Replicates.xlsx", 
                        sheet = "Both ")
mags_both_rep$Replicate <- as.factor(mags_both_rep$Replicate)
mags_both_rep$Replicate <- factor(mags_both_rep$Replicate, levels=c('LacAMPa','LacAMPb', 'LacAMPc', 'LacADNa','LacADNb','LacADNc'))
mags_both_rep$Coverage <- as.numeric(mags_both_rep$Coverage)
mags_both_rep$Family <- as.factor(mags_both_rep$Family)

mags_both_rep_1 <- mags_both_rep %>%
  mutate(Family = fct_reorder(Family, desc(Coverage))) %>%
  ggplot(aes(x=Family, y=Coverage, fill=Replicate)) + 
  geom_bar(position = position_dodge(0.8),
           width = 0.7, stat="identity") + 
  facet_wrap(~Family, scale="free") +
  #stat_compare_means( label =  "p.signif", label.x = 1.5) + 
  #guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
  scale_fill_manual(values=c("#0099FF",
                             "#0099FF",
                             "#0099FF",
                             "#003399",
                             "#003399",
                             "#003399")) +
  # change order of legend
  # Add a theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y= "Coverage of MAGs",
       x = "MAGs")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) 

mags_both_rep_1



## for only transporters
mags_trans <- read_excel("antismash/Replicates.xlsx", 
                        sheet = "Trans")
mags_trans$Group <- as.factor(mags_trans$Group)
mags_trans$Group <- factor(mags_trans$Group, levels=c('LacAMP', 'LacADN'))
mags_trans$Coverage <- as.numeric(mags_trans$Coverage)
mags_trans_1 <-ggplot(mags_trans, aes(x=Family, y=Coverage, fill=Group)) + 
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~Family, scale="free") +
  stat_compare_means( label =  "p.signif", label.x = 1.5) + 
  #guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
  scale_fill_brewer(palette="Paired") +
  # change order of legend
  # Add a theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y= "Mean Coverage of MAGs",
       x = "MAGs")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) 

mags_trans_1

## by replicates
mags_trans_rep <- read_excel("antismash/Replicates.xlsx", 
                            sheet = "Trans")
mags_trans_rep$Replicate <- as.factor(mags_trans_rep$Replicate)
mags_trans_rep$Replicate <- factor(mags_trans_rep$Replicate, levels=c('LacAMPa','LacAMPb', 'LacAMPc', 'LacADNa','LacADNb','LacADNc'))
mags_trans_rep$Coverage <- as.numeric(mags_trans_rep$Coverage)
mags_trans_rep$Family <- as.factor(mags_trans_rep$Family)

mags_trans_rep_1 <- mags_trans_rep %>%
  mutate(Family = fct_reorder(Family, desc(Coverage))) %>%
  ggplot(aes(x=Family, y=Coverage, fill=Replicate)) + 
  geom_bar(position = position_dodge(0.8),
           width = 0.7, stat="identity") + 
  facet_wrap(~Family, scale="free") +
  #stat_compare_means( label =  "p.signif", label.x = 1.5) + 
  #guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
  scale_fill_manual(values=c("#0099FF",
                             "#0099FF",
                             "#0099FF",
                             "#003399",
                             "#003399",
                             "#003399")) +
  # change order of legend
  # Add a theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y= "Coverage of MAGs",
       x = "MAGs")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) 


mags_trans_rep_1


################ test #########
# One box per treatment
p1 <- ggplot(mags_both, aes(x=Family, y=Coverage, fill=Group)) + 
  geom_boxplot() +
  facet_wrap(~Group) +
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  ylab("MAG Abundance  \n") + 
  ggtitle("Family distribution between treatments") + theme(plot.title = element_text(hjust = 0.5))
p1
# one box per variety
p2 <- ggplot(mags_both, aes(x=Family, y=Coverage, fill=Group)) + 
  geom_boxplot() +
  facet_wrap(~Family, scale="free") +
  #stat_compare_means( label =  "p.signif", label.x = 1.5) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
  ylab("Mean Coverage of MAGs") 
p2

#### Bubble plot ####
# Siderophore synthesis 
bb_gc_s <- read_excel("antismash/Replicates.xlsx", 
                        sheet = "GC_syn")
bb_gc_s$MAG <- as.factor(bb_gc_s$MAG)
levels(bb_gc_s$MAG)
#[1] "MAG1_Burkholderiaceae"     "MAG10_Pseudomonadaceae"    "MAG11_Burkholderiaceae"    "MAG5_Methylobacteriaceae" 
#[5] "MAG6_Nitrincolaceae"       "MAG7_Nitrincolaceae"       "MAG8_Propionibacteriaceae" "MAG9_Pseudomonadaceae" 
bb_gc_s$number <- 1

sum_bb_gc_s <- ddply(bb_gc_s,.(MAG, HMM), #find the mean of all data, list all separate groups
                    summarize, count=sum(number)) 
write_csv(sum_bb_gc_s, "sum_MAGs_sid_syn.csv")
sum_bb_gc_s %>%
  arrange(desc(count)) %>%
  mutate(MAG = fct_reorder(MAG, HMM, .fun='length' )) %>%
  ggplot(aes(x=MAG, y=HMM, size=count, fill=HMM)) +
  #scale_fill_discrete(labels=c('High Program', 'Low Program')) +
  geom_point(alpha=0.7) +
  coord_flip() +
  geom_point(fill = "cornflowerblue", color = "black", shape = 21) +
  scale_size(range = c(.1, 24), name="Gene count") + 
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
  ylab("Gene families for siderophore synthesis") +
  guides(fill = FALSE)

#Siderophore transporters
bb_gc_t <- read_excel("antismash/Replicates.xlsx", 
                      sheet = "GC_trans")
bb_gc_t$MAG <- as.factor(bb_gc_t$MAG)
levels(bb_gc_t$MAG)
#[1] "MAG1_Burkholderiaceae"     "MAG10_Pseudomonadaceae"    "MAG11_Burkholderiaceae"    "MAG13_Vibrionaceae"       
#[5] "MAG3_Marinomonadaceae"     "MAG4_Marinomonadaceae"     "MAG5_Methylobacteriaceae"  "MAG6_Nitrincolaceae"      
#[9] "MAG7_Nitrincolaceae"       "MAG8_Propionibacteriaceae" "MAG9_Pseudomonadaceae"  
bb_gc_t$number <- 1

sum_bb_gc_t <- ddply(bb_gc_t,.(MAG, HMM), #find the mean of all data, list all separate groups
                     summarize, count=sum(number)) 

sum_bb_gc_t %>%
  arrange(desc(count)) %>%
  mutate(MAG = fct_reorder(MAG, HMM, .fun='length' )) %>%
  ggplot(aes(x=MAG, y=HMM, size=count, fill=HMM)) +
  #scale_fill_discrete(labels=c('High Program', 'Low Program')) +
  geom_point(alpha=0.7) +
  coord_flip() +
  geom_point(fill = "cornflowerblue", color = "black", shape = 21) +
  scale_size(range = c(.1, 24), name="Gene count") + 
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
  ylab("Gene families for siderophore synthesis") +
  guides(fill = FALSE)

### Heatmap ####
# Siderophore synthesis

set.seed(123)  # for reproducibility
class(sum_bb_gc_s)
sum_bb_gc_s$HMM <- gsub("-siderophore-synthesis", "", sum_bb_gc_s$HMM)
sum_bb_gc_s$HMM <- gsub("-siderophore-utilization", "", sum_bb_gc_s$HMM)

wide_data <- spread(sum_bb_gc_s, MAG, count) ##100个MAGs(TOP10)
wide_data[is.na(wide_data)] <- 0 
write_csv(wide_data, "sid_syn_mag.csv")

sid_syn_mag <- read.csv("~/Desktop/BINGO2/MAGs/sid_syn_mag.csv", row.names=1, sep=";")
#trace <- read.csv2("plot/1.ed-top10-Fe-KEGG-FL.csv", row.names = 1)
sid_syn_mag <- sid_syn_mag[,-1]
sid_syn_mag <- as.data.frame(t(sid_syn_mag))

#wide_data2 <- wide_data
#wide_data2[is.na(wide_data2)] <- ""

# only focused on pvd and pvs gene families
sid_syn_mag1 <- sid_syn_mag %>% select(`Pvd-family`, `Pvs-family`)
sid_syn_mag1 <- as.matrix(sid_syn_mag1)
# Calculate row sums to get the counts
row_sums <- rowSums(sid_syn_mag1)

# Sort rows based on row sums (counts) in descending order
sorted_indices <- order(-row_sums)

# Reorder the data matrix based on the sorted indices
sorted_data <- data[sorted_indices, ]

p1 <- pheatmap(sid_syn_mag1, scale="none", fontsize_col = 10, na_col="white",
               #display_numbers = wide_data2, number_format = "%.0f",
               color = colorRampPalette(c("#DFE8F1", "navy"))(15),
               cluster_cols = FALSE, 
               cluster_rows = FALSE,
               main = "Siderophore synthesis gene families")
            

#complext heat map 
library(ComplexHeatmap)
Heatmap(sid_syn_mag1)

# Siderophore transport

set.seed(123)  # for reproducibility
class(sum_bb_gc_t)
sum_bb_gc_t$HMM <- gsub("-siderophore-receptor", "", sum_bb_gc_t$HMM)
sum_bb_gc_t$HMM <- gsub("-substrate-binding-protein", "", sum_bb_gc_t$HMM)
sum_bb_gc_t$HMM <- gsub("-siderophore-transport", "", sum_bb_gc_t$HMM)
sum_bb_gc_t$HMM <- gsub("-ATP-binding-protein", "", sum_bb_gc_t$HMM)
sum_bb_gc_t$HMM <- gsub("-siderophore receptor", "", sum_bb_gc_t$HMM)
sum_bb_gc_t$HMM <- gsub("-siderophore receptor", "", sum_bb_gc_t$HMM)
sum_bb_gc_t$HMM <- gsub("-sideropore-export", "", sum_bb_gc_t$HMM)
sum_bb_gc_t$HMM <- gsub("-siderophore-export", "", sum_bb_gc_t$HMM)
sum_bb_gc_t$HMM <- gsub("-permease", "", sum_bb_gc_t$HMM)


wide_data_t <- spread(sum_bb_gc_t, MAG, count) ##100个MAGs(TOP10)
wide_data_t[is.na(wide_data_t)] <- 0 
write_csv(wide_data_t, "sid_trans_mag.csv")

sid_trans_mag <- read.csv("~/Desktop/BINGO2/MAGs/sid_trans_mag.csv", row.names=1, sep=";")
#trace <- read.csv2("plot/1.ed-top10-Fe-KEGG-FL.csv", row.names = 1)
sid_trans_mag <- sid_trans_mag[,-1]
sid_trans_mag <- as.data.frame(t(sid_trans_mag))
####将缺失值不显示
#wide_data2 <- wide_data
#wide_data2[is.na(wide_data2)] <- ""

# only focused on pvd and pvs gene families
sid_trans_mag1 <- sid_trans_mag %>% select(`Fpv-family`, `Pvu-family`)
sid_trans_mag1 <- as.matrix(sid_trans_mag1)

# Calculate row sums to get the counts
row_sums1 <- rowSums(sid_trans_mag1)

# Sort rows based on row sums (counts) in descending order
sorted_indices1 <- order(-row_sums1)

# Reorder the data matrix based on the sorted indices
sorted_data1 <- sid_trans_mag1[sorted_indices1, ]

p2 <- pheatmap(sid_trans_mag1, scale="none", fontsize_col = 10, na_col="white",
               #display_numbers = wide_data2, number_format = "%.0f",
               color = colorRampPalette(c("#DFE8F1", "navy"))(10),
               cluster_cols = FALSE, 
               cluster_rows = FALSE, 
               legend_labels = NULL,  # Set legend_labels to NULL
               main = "Siderophore transport gene families")

#complext heat map 
library(ComplexHeatmap)
Heatmap(sid_trans_mag1)
#library(circlize) # >= 0.4.10
#split = factor(sid_syn_mag, levels = sid_syn_mag[1:8])
#col_fun1 = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
#circos.heatmap(sid_syn_mag, split = split, col = col_fun1)

## combine heatmaps
lm <- rbind(c(1,2))
grid.arrange(grobs = list(p1[[4]],
                          p2[[4]]),
             layout_matrix = lm)

## Heatmap based on type of siderophore ####
# Sets the minimum (0), the maximum (15), and the increasing steps (+1) for the color scale
# Note: if some of your genes are outside of this range, they will appear white on the heatmap
breaksList = seq(0, 15, by = 1)
pheatmap(expressionData[1:10, ], # Plots the first 10 genes of the dataset
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList) # Sets the breaks of the color scale as in breaksList


## Pyoverdine 
sid <- read.csv("~/Desktop/BINGO2/MAGs/Siderophore_grouped.csv", row.names=1, sep=";")
sid <- sid[,-1]
sid <- as.data.frame(t(sid))
# only focused on Pvd and Fpv families
sid1 <- sid %>% select(`Pvd-family `, `Fpv-family`)
sid1 <- as.matrix(sid1)

# Calculate row sums to get the counts
row_sums_1 <- rowSums(sid1)

# Sort rows based on row sums (counts) in descending order
sorted_indices_1 <- order(-row_sums_1)

# Reorder the data matrix based on the sorted indices
sorted_data_1 <- sid1[sorted_indices_1, ]

sid1_ordered_both <- sid1[order(rowMeans(sid1), decreasing = TRUE), 
                          order(colMeans(sid1), decreasing = TRUE)]

hp_1 <- pheatmap(sid1_ordered_both, scale="none", fontsize_col = 12, 
                 #na_col="white",
               #display_numbers = wide_data2, number_format = "%.0f",
               color = colorRampPalette(c("white", "#4292c6", "#2171b5" ,"#084594"))(25),
               breaks = seq(0, 25, by = 1),   # Breaks from 0 to 25
               legend_breaks = seq(0, 25, by = 5),
               border_color = "black",
               cluster_cols = FALSE, 
               cluster_rows = FALSE, 
               legend_labels = NULL,  # Set legend_labels to NULL
               main = "Pyoverdine")
hp_1

## Vibrioferrin 
sid <- read.csv("~/Desktop/BINGO2/MAGs/Siderophore_grouped.csv", row.names=1, sep=";")
sid <- sid[,-1]
sid <- as.data.frame(t(sid))
# only focused on Pvd and Fpv families
sid2 <- sid %>% select(`Pvs-family `, `Pvu-family`)
sid2 <- sid2[rownames(sid2) != "MAG5_Methylobacteriaceae", ]
sid2 <- as.matrix(sid2)

# Calculate row sums to get the counts
row_sums_2 <- rowSums(sid2)

# Sort rows based on row sums (counts) in descending order
sorted_indices_2 <- order(-row_sums_2)

# Reorder the data matrix based on the sorted indices
sorted_data_2 <- sid2[sorted_indices_2, ]

sid2_ordered_both <- sid2[order(rowMeans(sid2), decreasing = TRUE), 
                        order(colMeans(sid2), decreasing = TRUE)]


hp_2 <- pheatmap(sid2_ordered_both, scale="none", fontsize_col = 12, na_col="white",
                 #display_numbers = wide_data2, number_format = "%.0f",
                 color = colorRampPalette(c("white", "#c6dbef", "#9ecae1","#6baed6", "#4292c6", "#2171b5", "#084594" ))(10),
                 cluster_cols = FALSE, 
                 cluster_rows = FALSE, 
                 border_color = "black", 
                 legend_labels = NULL,  # Set legend_labels to NULL
                 main = "Vibrioferrin")

## combine heatmaps
lm_1 <- rbind(c(1,2))
grid.arrange(grobs = list(hp_1[[4]],
                          hp_2[[4]]),
             layout_matrix = lm_1)

