
#####################. 
####### Use Primer clipped files from data
##############
library(dada2)
library(pals)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(forcats)
library(readxl)
library(phyloseq)
library(RColorBrewer)
library(vegan)
library(reshape2)
library(reshape)
library(knitr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(hrbrthemes)
library(Biostrings)


##########

##### on terminal #### cp */*_R1.fastq.bz2 Reads ####
path <- "/Users/rmthoppil/Desktop/BINGO1/16s-bingo/NGS2481-2021-06-10-1257/PrimeClipped_Rhea/Reads" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
####
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq.bz2", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.bz2", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

##### Change the forward and reverse base reads based on read quality plots 
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
#### Sample inference 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
#######
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
## Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
#### Construct sequence table 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#### Remove chimaeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
###### Has to be around 0.9
sum(seqtab.nochim)/sum(seqtab)

ncol(seqtab.nochim)
write.csv(seqtab.nochim, "count_data.csv")
 write(seqtab.nochim, "count_data.txt")

 
 write.table(seqtab.nochim, file = "count_data.txt", sep="\t", quote=F)
 
#### Track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

###### Assign taxonomy using SILVA
taxa <- assignTaxonomy(seqtab.nochim, "/Users/rmthoppil/Desktop/SIDEMAR/SIDEMAR_2/ASVs/GTDB/GTDB_bac-arc_ssu_r86.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/rmthoppil/Desktop/SIDEMAR/SIDEMAR_2/ASVs/GTDB/GTDB_dada2_assignment_species.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa.print, "taxa.print.csv")

#### Phyloseq for metabarcoding data ####
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

theme_set(theme_bw())

samples.out <- rownames(seqtab.nochim)

samples_inf <- read.csv("~/Desktop/BINGO1/ASVs/ASV_bingo-result/bingo_sample_informations-1.28.csv", row.names=1, sep=";")

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samples_inf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "515FY-926Rjed-LGC-Positivkontrolle-Zymo", ps) # Remove mock sample


dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

library(phyloseq)
otu_table(ps)
tax_table(ps)
sample_data(ps)
refseq(ps)


genusabundance <- ps %>%
  tax_glom(taxrank = "Genus") %>%                        # Set to smallest taxonomic level you are interested in
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()  %>%                                             # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Genus) 
head(genusabundance)

Kingdomabundance <- ps %>%
  tax_glom(taxrank = "Kingdom") %>%                        # Set to smallest taxonomic level you are interested in
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()  %>%                                             # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Kingdom) 
head(Kingdomabundance)


# We can use the `sample_sums()` function of phyloseq to add up total read counts for each sample.
sample_sum_df <- data.frame(sum = sample_sums(ps))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

## ASV - all ####
genusabundance_1 <- ps %>%
  tax_glom(taxrank = "Genus") %>%                        # Set to smallest taxonomic level you are interested in
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()  %>%                                             # Melt to long format
  arrange(Genus) 

head(genusabundance_1)
write.csv(genusabundance_1, "ASV_all.csv")


#### ASV > 1 % abundance ####
write.csv(genusabundance, "Relative_abundance_Genus_BINGO1.csv")

Relative_abundance_Genus_BINGO1 <- read_excel("Relative_abundance_Genus_BINGO1.xlsx", 
                                                   sheet = "clean")
rb <- subset(Relative_abundance_Genus_BINGO1, Group == "Glacial" | Group == "Non-glacial" )
rb$Abundance <- as.numeric(rb$Abundance)
# Define the number of colors you want
nb.cols <- 16
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
# Create a ggplot with 16 colors 

#rb2 <- rb1 %>%
      #mutate(name = fct_relevel(Sample, 
                          #  "LacAMPa", "LacAMPb", "LacAMPc", 
                           # "LacADNa", "LacADNb", "LacADNc"))

rb1$Group <- as.factor(rb1$Group)

rb3 <- ggplot(rb, aes(x = Sample, y = Abundance, fill = Family)) + 
      geom_bar(position = "stack", stat = "identity")+
      theme_bw() + 
      scale_fill_manual(values = mycolors) +
      scale_x_discrete(drop = TRUE)+
      facet_wrap(~Group, ncol = 2, dir = "v", scales = "free_x") +
      theme(
      strip.text.x = element_text(size = 15, color = "black", face = "bold"),
      strip.text.y = element_text(size = 15, color = "black", face = "bold"),
      strip.background = element_rect(color="black", fill="white", linewidth=2, linetype="solid")) +
      theme(text = element_text(size = 15, color = "black"), 
            axis.text.x = element_text(size = 15, face = "plain") , 
            axis.text.y = element_text(size = 15, face = "plain")) +
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
      ylab("Relative Abundance (> 1%) \n") + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      guides(fill = guide_legend(ncol = 2))
rb3

## BINGO1 Relative abundance + control ####
rc <- subset(Relative_abundance_BINGO1, Type == "glacial" | Type == "non-glacial" | Sample == "Ctr-tf-IIb" )
write.csv(rc, "Relative_control.csv")

rc1 <- read_excel("Relative_control.xlsx")
rc1$Abundance <- as.numeric(rc1$Abundance)
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 17
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

rc2 <- rc1 %>%
  mutate(name = fct_relevel(Sample, 
                            "Control_1","LacAMPa", "LacAMPb", "LacAMPc", 
                            "LacADNa", "LacADNb", "LacADNc"))

rc3 <- ggplot(rc2, aes(x = name, y = Abundance, fill = Family)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = mycolors) +
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
  ylab("Relative Abundance (Family > 1%) \n") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_legend(ncol = 2))
rc3

## BINGO1 Relative abundance + Table bay ####
rt <- subset(Relative_abundance_BINGO1, Type == "glacial" | Type == "non-glacial" | Sample == "BdT0-2" )
write.csv(rt, "Relative_tb.csv")

rt1 <- read_excel("Relative_tb.xlsx")
rt1$Abundance <- as.numeric(rt1$Abundance)
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

rt2 <- rt1 %>%
  mutate(name = fct_relevel(Sample, 
                            "LacAMPa", "LacAMPb", "LacAMPc", 
                            "LacADNa", "LacADNb", "LacADNc","Baie de la Table"))

rt3 <- ggplot(rt2, aes(x = Sample, y = Abundance, fill = Family)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = mycolors) +
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
  ylab("Relative Abundance (> 1%) \n") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_legend(ncol = 2))
rt3

## BINGO1 Relative abundance + Table bay + Control_1 ####
ry <- subset(Relative_abundance_Genus_BINGO1, Group == "Glacial" | Group == "Non-glacial" | Sample == "BdT0-2" | Sample == "Control_1")
ry$Abundance <- as.numeric(ry$Abundance)
# Define the number of colors you want
nb.cols <- 23
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

#ry2 <- ry %>%
 # mutate(name = fct_relevel(Sample, 
  #                          "LacAMPa", "LacAMPb", "LacAMPc", 
   #                         "LacADNa", "LacADNb", "LacADNc","Baie de la Table", "Control_1"))

ry3 <- ggplot(ry, aes(x = Sample, y = Abundance, fill = Family)) + 
  theme_classic() +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = mycolors) +
  scale_x_discrete(drop = TRUE)+
  theme(text = element_text(size = 15, color = "black"), 
        axis.text.x = element_text(size = 15, face = "plain") , 
        axis.text.y = element_text(size = 15, face = "plain"))  + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
  ylab("Relative Abundance (> 1%) \n") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_legend(ncol = 2))
ry3

#### NMDS + ANOSIM - BINGO1 #####
#rarefy_even_depth downsamples/normalizes all samples to the same depth and prunes OTUs that disappear from all samples as a result
ps1 <- subset_samples(ps, Project =="BINGO1")#Keep only samples to be analyzed

sample_sum_df1 <- data.frame(sum = sample_sums(ps1))
physeq1.rarefied = rarefy_even_depth(ps1,sample.size=min(sample_sums(ps1)))
##1590 OTU's removed
#record 

physeq.pop <- transform_sample_counts(physeq1.rarefied, function(otu) otu/sum(otu)) # Transform data to proportions as appropriate for Bray-Curtis distances
ord.nmds.bray <- ordinate(physeq.pop, method="NMDS", distance="bray")

samples <- sample_data(samples_inf)
colourCount = length(unique(samples$Lake))
getPalette = colorRampPalette(brewer.pal(12,"Paired"))

p <- plot_ordination(physeq.pop, ord.nmds.bray, 
                     type="samples", 
                     color="Lake", 
                     shape = "Type", 
                     title="Bray NMDS-BINGO") + theme_bw()
p                                                                      
#
p1 <- p + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = getPalette(colourCount)) +
  geom_point(size = 5)
p1

### only glacial and non-glacial samples
sub1_physeq1.rarefied <- subset_samples(ps1, Lake =="Lake ADN" | Lake == "LakeAMP")
ps_otu<-otu_table(sub1_physeq1.rarefied)
ps_otu_transformed<-decostand(ps_otu, method = "hellinger")
ps_otu_transformed<-ps_otu_transformed[, colSums(ps_otu_transformed)!=0]
ncol(ps_otu_transformed) 

# calculate Bray-Curtis dissimilarity
ps_otu_dist<-vegdist(ps_otu_transformed, method = "bray")
ps_otu_dist_clust<-hclust(ps_otu_dist,method="average")
plot(ps_otu_dist_clust, ylab="ps_otu dist.", hang=-1,cex = 0.7)
ps_otu_NMDS<-monoMDS(ps_otu_dist, k=2, model = "hybrid")
NMDS1 =ps_otu_NMDS$points[,1]
NMDS2 =ps_otu_NMDS$points[,2]

NMDS = data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, 
                  Lake=sample_data(sub1_physeq1.rarefied)$Lake, 
                  ID=sample_data(sub1_physeq1.rarefied)$Number,
                  Shape=sample_data(sub1_physeq1.rarefied)$Type, 
                  Label=sample_data(sub1_physeq1.rarefied)$Number)
head(NMDS)

# global test
anosim_filter <- anosim(ps_otu_dist, NMDS$Lake, permutations = 9999, distance = "bray")
anosim_filter #ANOSIM statistic R: 0.5185, Significance: 0.1

## PREVIOUS ANOSIM - Anosim results are: ANOSIM statistic R: 0.5455,  P- Significance: 0.023

### SIMPER - Rhea ####
ps_otu_transformed_a <- as.data.frame(ps_otu_transformed)
write_csv(ps_otu_transformed_a, "ps_otu.csv")
ps_otu <- read_delim("ps_otu.csv", delim = ";", 
                     escape_double = FALSE, trim_ws = TRUE)
ps_otu_1 <- ps_otu %>% remove_rownames %>% column_to_rownames(var="Replicate")

simp_2 <- simper(ps_otu_1)
summary(simp_2)




### Flow cytometry ####
s1 <- read_excel("FlowCytometrie_Incub_BINGO_1_Rhea.xlsx", 
                     sheet = "Clean")
summary(s1)
s1 <- as.data.frame(s1)
s1_1 <- melt(setDT(s1[2:6]), id.vars = c("Samples"), variable.name = "Days")
 # Create the plot
ggplot(s1_1, aes(x = Days, y = value, color = factor(Samples))) +
  geom_line(aes(group = Samples), size = 1) +
  scale_color_manual(values = c("#ca7dcc","#ca7dcc","#ca7dcc",
                                "#cc0000","#cc0000","#cc0000",
                                "#1b98e0","#1b98e0","#1b98e0")) +
  geom_point(aes(shape = factor(Samples), size = 2) , fill="white", color="black") +
  scale_shape_manual(values = c(19, 19, 19, 19, 19, 19, 19, 19, 19)) +
  scale_size_identity() +
  labs(x = "Time (d)", y = "Prokaryotic abundance") +
  theme_minimal() +
  theme(legend.title = element_blank()) 

##  plot with group
write.csv(s1_1, "flow_cytometry.csv")
s1_2 <- read_delim("flow_cytometry.csv", 
                       delim = ";", escape_double = FALSE, trim_ws = TRUE)
s1_2$Days <- factor(x = s1_2$Days, 
                      levels = c("D0","D1","D2","D3", "D4", "D5","D6","D7","D8","D9"))

# Create a line plot faceted by bacterial species
s1_3 <- ggplot(s1_2, aes(x = Days, y = value, color = factor(Samples))) +
  geom_line(aes(group = Samples), size = 1) +
  scale_color_manual(values = c("#ca7dcc","#ca7dcc","#ca7dcc",
                                "#cc0000","#cc0000","#cc0000",
                                "#1b98e0","#1b98e0","#1b98e0")) +
  geom_point(shape=21, color="black", fill="#69b3a2", size=3) +
  #scale_shape_manual(values = c(19, 19, 19, 19, 19, 19, 19, 19, 19)) +
  scale_size_identity() +
  labs(x = "Time (d)", y = "Prokaryotic abundance") +
  facet_wrap(~Group) +
  #theme_ipsum() +
  theme(legend.title = element_blank())

s1_3

## dot plots ####
n1 <- read_excel("FlowCytometrie_Incub_BINGO_1_Rhea.xlsx", 
                 sheet = "new")
summary(n1)
n1 <- as.data.frame(n1)
n1_1 <- melt(setDT(n1[2:12]), id.vars = c("Samples"), variable.name = "Days")
ggplot(n1_1, aes(x = Days, y = value, shape = Samples)) +
  geom_point(color="darkblue", size=2) +          # Set the size of the dots
  scale_shape_manual(values = c(16,16,16,
                                17,17, 17,
                                18, 18, 18)) +  # Different shapes for each group
  theme_minimal() +
  #facet_wrap(~Group) +
  labs( x = "Time (d)",  y = "Prokaryotic abundance")

# facet by group
write.csv(n1_1, "flow_cytometry_new.csv")
n1_2 <- read.csv2("~/Desktop/BINGO1/ASVs/ASV_bingo-result/flow_cytometry_new.csv")
# Values you want to remove (e.g., removing values equal to 15 and 40)
filtered_n1 <- n1_2[!n1_2$value %in% c(0), ]
n1_2$Days <- factor(x = n1_2$Days, 
                    levels = c("D0","D1","D2","D3", "D4", "D5","D6","D7","D8","D9"))
ggplot(n1_2, aes(x = Days, y = ifelse(value >= 1, value, NA), shape = Replicates)) +
  geom_point(color="darkblue", size=4) +          # Set the size of the dots
  scale_shape_manual(values = c(16,16,16,
                                17,17, 17,
                                18, 18, 18)) +  # Different shapes for each group
  #theme_minimal() +
  facet_wrap(~Group) +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  labs( x = "Time (d)",  y = "Prokaryotic abundance (cells/ml)")

## line plot
n1_3 <- ggplot(n1_2, aes(x = Days, y = ifelse(value >= 1, value, NA), color = factor(Replicates))) +
  geom_line(aes(group = Replicates), size = 1) +
  scale_color_manual(values = c("#ca7dcc","#ca7dcc","#ca7dcc",
                                "#cc0000","#cc0000","#cc0000",
                                "#1b98e0","#1b98e0","#1b98e0")) +
  geom_point(shape=21, color="black", fill="#69b3a2", size=3) +
  #scale_shape_manual(values = c(19, 19, 19, 19, 19, 19, 19, 19, 19)) +
  scale_size_identity() +
  labs(x = "Time (d)", y = "Prokaryotic abundance") +
  facet_wrap(~Group) +
  #theme_ipsum() +
  theme(legend.title = element_blank())

s1_3
#statistics ####
## focus on the the D9 (last day prokaryotic abundance was measured)
library(stats)
s1_4 <- subset(s1_2, Days == "D9")

anova_1<-aov(value~Group, data=s1_4)
summary(anova_1)
lm_1 <-lm(value~Group,data=s1_4)
summary(lm_1)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)        660320      50908  12.971 1.29e-05 ***
#  Group-Glacial       212420      71995   2.950   0.0256 *  
#  Group-Non-glacial    72164      71995   1.002   0.3549 

s1_4$Group <- as.factor(s1_4$Group)
contrasts(s1_4$Group)

# Run model without intercept
model3<-lm(value~Group-1, data=s1_4)
summary(model3)

par(mfrow=c(2,2))
plot(anova_1)
plot(lm_1)


##Check normality and equal variances (homocedasticity)
##Shapiro test
shapiro.test(s1_4$value[s1_4$Group=="Glacial"])
shapiro.test(s1_4$value[s1_4$Group=="Non-glacial"])

##Barlett's test for equal variances
bartlett.test(s1_4$value~s1_4$Group, data=s1_4)

#### t-test - group
attach(s1_4)
s1_4 <- as.data.frame(s1_4)

s1_4 %>% t_test(value ~ Group,
                        var.equal = TRUE,
                        detailed = TRUE) %>%
  glimpse()

# Statistical test - T-test
stat.test_1 <- s1_4 %>%
  t_test(value ~ Group) %>%
  add_significance()
stat.test_1

cm <- compare_means(value ~ Group, data = s1_4, method = "t.test")
cm
#.y.   group1  group2           p p.adj p.format p.signif method
#<chr> <chr>   <chr>        <dbl> <dbl> <chr>    <chr>    <chr> 
#1 value Control Glacial     0.0294 0.088 0.029    *        T-test
#2 value Control Non-glacial 0.396  0.4   0.396    ns       T-test
#3 value Glacial Non-glacial 0.166  0.33  0.166    ns       T-test
