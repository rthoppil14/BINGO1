
#####################. 
####### Use Primer clipped files from data
##############
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.16", force = TRUE)
library(dada2)
library(pals)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(forcats)
library(readxl)
library(vegan)
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
taxa <- assignTaxonomy(seqtab.nochim, "/Users/rmthoppil/Desktop/BINGO2_DNA_extraction/SILVA/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
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
head(samples.out)
ncol(samples.out)
samples <- as.data.frame(samples.out)
sample <- sapply(strsplit(samples.out, "D"), `[`, 1)
treat <- substr(sample,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Sample=sample, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
samples_df <- read.csv("~/Desktop/BINGO1/ASV_bingo-result/bingo_sample_informations-1.28.csv", sep=";")
samples_df <- samples_df %>% 
  column_to_rownames(var = "Samples")
ncol(samples_df)
rownames(samples_df) <- samples.out

sample_names(samples_df)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samples_df), 
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

## ASV - all ####
genusabundance_1 <- ps %>%
  tax_glom(taxrank = "Genus") %>%                        # Set to smallest taxonomic level you are interested in
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()  %>%                                             # Melt to long format
  arrange(Genus) 
head(genusabundance_1)
df <- subset(genusabundance_1, Type == "glacial" | Type == "non-glacial" | Sample == "BdT0-2")
write.csv(df, "ASV_all.csv")


#### ASV > 1 % abundance ####
colourCount = length(unique(genusabundance$Family))
getPalette = colorRampPalette(brewer.pal(12,"Paired"))

write.csv(genusabundance, "Relative_abundance_Genus_BINGO1.csv")

Relative_abundance_BINGO1 <- read.csv("~/Desktop/BINGO1/ASV_bingo-result/Relative_abundance_Genus_BINGO1.csv")
Relative_abundance_BINGO1$Abundance <- as.numeric(Relative_abundance_BINGO1$Abundance)
rb <- subset(Relative_abundance_BINGO1, Type == "glacial" | Type == "non-glacial" )
write.csv(rb, "Relative_treatment.csv")

rb1 <- read_excel("Relative_treatment.xlsx")
rb1$Abundance <- as.numeric(rb1$Abundance)
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 17
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
# Create a ggplot with 17 colors 

#rb2 <- rb1 %>%
      #mutate(name = fct_relevel(Sample, 
                          #  "LacAMPa", "LacAMPb", "LacAMPc", 
                           # "LacADNa", "LacADNb", "LacADNc"))

rb1$Group <- as.factor(rb1$Group)

rb3 <- ggplot(rb1, aes(x = Sample, y = Abundance, fill = Family)) + 
      geom_bar(position = "stack", stat = "identity")+
      theme_bw() + 
      scale_fill_manual(values = mycolors) +
      scale_x_discrete(drop = TRUE)+
      facet_wrap(~Group, ncol = 2, dir = "v", scales = "free_x") +
      theme(
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_text(size = 10, color = "black", face = "bold"),
      strip.background = element_rect(color="black", fill="white", size=2, linetype="solid")) +
      theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
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
                            "Table_B","LacAMPa", "LacAMPb", "LacAMPc", 
                            "LacADNa", "LacADNb", "LacADNc"))

rt3 <- ggplot(rt2, aes(x = name, y = Abundance, fill = Family)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = mycolors) +
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
  ylab("Relative Abundance (Family > 1%) \n") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_legend(ncol = 2))
rt3

#### NMDS + ANOSIM - BINGO1 #####

#rarefy_even_depth downsamples/normalizes all samples to the same depth and prunes OTUs that disappear from all samples as a result
physeq1.rarefied = rarefy_even_depth(ps,sample.size=min(sample_sums(ps)))
#record 

sub1_physeq1.rarefied <- subset_samples(physeq1.rarefied, Treat =="incubation experiments")#Keep only samples to be analyzed
sub1_physeq1.rarefied1 <- subset_samples(physeq1.rarefied, Type =="glacial"  | Type == "non-glacial")#Keep only samples to be analyzed

physeq.pop <- transform_sample_counts(sub1_physeq1.rarefied, function(otu) otu/sum(otu)) # Transform data to proportions as appropriate for Bray-Curtis distances
ord.nmds.bray <- ordinate(physeq.pop, method="NMDS", distance="bray")

p <- plot_ordination(physeq.pop, ord.nmds.bray, type="samples", color="Fraction", shape="Lake", title="Bray NMDS-SWINGS") + theme_bw()
p                                                                      
#
p1 <- p + theme(plot.title = element_text(hjust = 0.5))
p1

ps_otu<-otu_table(sub1_physeq1.rarefied1)
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

#sample_data(sub1_physeq1.rarefied)$Treat <- factor(sample_data(sub1_physeq1.rarefied)$Treat, 
#                                                   levels = c("STSW","SASW","PFSW","ASW",
#                                                          "STSW/STUW","WW","ASLOW","STMW","AAIW"))
                                                                                                                                  
NMDS = data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, 
                  Treat=sample_data(sub1_physeq1.rarefied1)$Type, 
                  ID=sample_data(sub1_physeq1.rarefied1)$Number,
                  Shape=sample_data(sub1_physeq1.rarefied1)$Type, 
                  Label=sample_data(sub1_physeq1.rarefied1)$Number)
head(NMDS)

# global test
anosim_filter <- anosim(ps_otu_dist, NMDS$Treat, permutations = 9999, distance = "bray")
anosim_filter #ANOSIM statistic R: 0.5455,  P- Significance: 0.023

