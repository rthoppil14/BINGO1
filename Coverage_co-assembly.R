#Coverage - co-assembly
library(readxl)
library(tidyverse)
library(plyr)
library(ggplot2)
library(viridis)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(hrbrthemes)
library(ggpubr)
library(patchwork)
library(cowplot)
library(ggpubr)
library(rstatix)
library(pals)
library(forcats)

### Merge output from name_conversions, Fe_Genie output and coverage file ####
name_conversions <- read.delim("~/Desktop/BINGO2/Co-assembly/Gene abundance/Coverage_zip/name_conversions.txt", header=FALSE)
colnames(name_conversions)[1] <- "contig"
FeGenie.geneSummary <- read.csv("~/Desktop/BINGO2/Co-assembly/FeGenie/FeGenie-geneSummary.csv", sep=";")
colnames(FeGenie.geneSummary)[3] <- "contig"



###### ADN samples - Sample 64 ####
ADN_64profile <- read_delim("~/Desktop/BINGO2/Co-assembly/Gene abundance/Coverage_zip/ADN_64profile.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
co_64 <- merge(name_conversions, ADN_64profile, by="contig")
quant_64 <- read.delim("~/Desktop/BINGO2/Co-assembly/Salmon/quant_64.sf")
colnames(quant_64)[1] <- "contig"
co_64_fe <- merge(FeGenie.geneSummary, quant_64, by="contig")
co_64_fe$TPM <- as.numeric(co_64_fe$TPM)

colnames(co_64_fe)[14] <- "TPM_64"
sid_syn_tax_fe1 <- merge(co_64_fe, sid_syn_tax, by="contig")
co64 <- sid_syn_tax_fe %>%                                 # Group data
  group_by(HMM.x, Taxa) %>%
  dplyr::summarize(sumTPM = sum(TPM_64)) %>% 
  as.data.frame()
co64$HMM.x <- gsub("-siderophore-synthesis", "", co64$HMM.x)
co64_1 <- ggplot(co64, aes(x = HMM.x, y = sumTPM, fill = Taxa)) +
  geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values=as.vector(cols25(n=14))) + coord_flip()
co64_1 

###### ADN samples - Sample 65
ADN_65profile <- read_delim("~/Desktop/BINGO2/Co-assembly/Gene abundance/Coverage_zip/ADN_65profile.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
co_65 <- merge(name_conversions, ADN_65profile, by="contig")
quant_65 <- read.delim("~/Desktop/BINGO2/Co-assembly/Salmon/quant_65.sf")
colnames(quant_65)[1] <- "contig"
co_65_fe <- merge(FeGenie.geneSummary, quant_65, by="contig")
co_65_fe$TPM <- as.numeric(co_65_fe$TPM)
colnames(co_65_fe)[14] <- "TPM_65"
sid_syn_tax_fe2 <- merge(co_65_fe, sid_syn_tax, by="contig")
co65 <- sid_syn_tax_fe2 %>%                                 # Group data
  group_by(HMM.x, Taxa) %>%
  dplyr::summarize(sumTPM = sum(TPM_65)) %>% 
  as.data.frame()
co65$HMM.x <- gsub("-siderophore-synthesis", "", co65$HMM.x)
co65_1 <- ggplot(co65, aes(x = HMM.x, y = sumTPM, fill = Taxa)) +
  geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values=as.vector(cols25(n=14))) + coord_flip()
co65_1

###### ADN samples - Sample 66
ADN_66profile <- read_delim("~/Desktop/BINGO2/Co-assembly/Gene abundance/Coverage_zip/ADN_66profile.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
co_66 <- merge(name_conversions, ADN_66profile, by="contig")
quant_66 <- read.delim("~/Desktop/BINGO2/Co-assembly/Salmon/quant_66.sf")
colnames(quant_66)[1] <- "contig"
co_66_fe <- merge(FeGenie.geneSummary, quant_66, by="contig")
co_66_fe$TPM <- as.numeric(co_66_fe$TPM)
colnames(co_66_fe)[14] <- "TPM_66"
sid_syn_tax_fe3 <- merge(co_66_fe, sid_syn_tax, by="contig")
co66 <- sid_syn_tax_fe3 %>%                                 # Group data
  group_by(HMM.x, Taxa) %>%
  dplyr::summarize(sumTPM = sum(TPM_66)) %>% 
  as.data.frame()
co66$HMM.x <- gsub("-siderophore-synthesis", "", co66$HMM.x)
co66_1 <- ggplot(co65, aes(x = HMM.x, y = sumTPM, fill = Taxa)) +
  geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values=as.vector(cols25(n=14))) + coord_flip()
co66_1

######### AMP samples #####
###### AMP samples - Sample 67
ADN_67profile <- read_delim("~/Desktop/BINGO2/Co-assembly/Gene abundance/Coverage_zip/AMP_67profile.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
co_67 <- merge(name_conversions, ADN_67profile, by="contig")
quant_67 <- read.delim("~/Desktop/BINGO2/Co-assembly/Salmon/quant_67.sf")
colnames(quant_67)[1] <- "contig"
co_67_fe <- merge(FeGenie.geneSummary, quant_67, by="contig")
co_67_fe$TPM <- as.numeric(co_67_fe$TPM)
co67 <- co_67_fe %>%                                 # Group data
  group_by(HMM, category) %>%
  dplyr::summarize(sumTPM = sum(TPM)) %>% 
  as.data.frame()
colnames(co67)[3] <- "sumTPM_67"

write_csv(co67, "co67.csv")

###### AMP samples - Sample 68
ADN_68profile <- read_delim("~/Desktop/BINGO2/Co-assembly/Gene abundance/Coverage_zip/AMP_68profile.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
co_68 <- merge(name_conversions, ADN_68profile, by="contig")
quant_68 <- read.delim("~/Desktop/BINGO2/Co-assembly/Salmon/quant_68.sf")
colnames(quant_68)[1] <- "contig"
co_68_fe <- merge(FeGenie.geneSummary, quant_68, by="contig")
co_68_fe$TPM <- as.numeric(co_68_fe$TPM)
co68 <- co_68_fe %>%                                 # Group data
  group_by(HMM, category) %>%
  dplyr::summarize(sumTPM = sum(TPM)) %>% 
  as.data.frame()
colnames(co68)[3] <- "sumTPM_68"

write_csv(co68, "co68.csv")

###### AMP samples - Sample 69
ADN_69profile <- read_delim("~/Desktop/BINGO2/Co-assembly/Gene abundance/Coverage_zip/AMP_69profile.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
co_69 <- merge(name_conversions, ADN_69profile, by="contig")
quant_69 <- read.delim("~/Desktop/BINGO2/Co-assembly/Salmon/quant_69.sf")
colnames(quant_69)[1] <- "contig"
co_69_fe <- merge(FeGenie.geneSummary, quant_69, by="contig")
co_69_fe$TPM <- as.numeric(co_69_fe$TPM)
co69 <- co_69_fe %>%                                 # Group data
  group_by(HMM, category) %>%
  dplyr::summarize(sumTPM = sum(TPM)) %>% 
  as.data.frame()
colnames(co69)[3] <- "sumTPM_69"

write_csv(co69, "co69.csv")

#### For FeGenie_metagenome community analysis #######
## Siderophore synthesis 
sid_syn_tax <- read.csv("~/Desktop/BINGO2/Co-assembly/FeGenie/sid_syn_tax.csv", sep=";")
sid_syn_tax_fe <- merge(co_64_fe, sid_syn_tax, by="contig")
sid_syn_tax_fe <- merge(sid_syn_tax_fe, co_65_fe, by = "contig")
sid_syn_tax_fe <- merge(sid_syn_tax_fe, co_66_fe, by = "HMM")
sid_syn_tax_fe <- merge(sid_syn_tax_fe, co_67_fe, by = "HMM")
sid_syn_tax_fe <- merge(sid_syn_tax_fe, co_68_fe, by = "HMM")
sid_syn_tax_fe <- merge(sid_syn_tax_fe, co_69_fe, by = "HMM")

sid_syn_tax_fe <- sid_syn_tax_fe %>% select(-category.x, -category.y, -contig)
write_csv(sid_syn_tax_fe, "sid_syn_fe_TPM.csv")


## Siderophore transport 
sid_trans_tax <- read.csv("~/Desktop/BINGO2/Co-assembly/FeGenie/sid_trans_tax.csv", sep=";")
sid_trans_tax_fe <- merge(co64, sid_trans_tax, by="HMM")
sid_trans_tax_fe <- merge(sid_trans_tax_fe, co65, by = "HMM")
sid_trans_tax_fe <- merge(sid_trans_tax_fe, co66, by = "HMM")
sid_trans_tax_fe <- merge(sid_trans_tax_fe, co67, by = "HMM")
sid_trans_tax_fe <- merge(sid_trans_tax_fe, co68, by = "HMM")
sid_trans_tax_fe <- merge(sid_trans_tax_fe, co69, by = "HMM")
sid_trans_tax_fe <- sid_trans_tax_fe %>% select(-category.x, -category.y, -contig)

write_csv(sid_trans_tax_fe, "sid_trans_fe_TPM.csv")

########## based on categories wihtout 69 ##################
Heme_T <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories.xlsx", 
 sheet = "Heme transport")
Heme_T$sumTPM <- as.numeric(Heme_T$sumTPM)
Heme_T$HMM <-gsub("-substrate-binding-protein","",as.character(Heme_T$HMM))
Heme_T$HMM <-gsub("-ATP-binding-protein","",as.character(Heme_T$HMM))
Heme_T <- ggplot(Heme_T) +
  aes(x = Sample, y = sumTPM, fill = HMM) +
  geom_bar(position = "stack", stat = "identity", show.legend = TRUE,  width = 0.98) +
  scale_fill_brewer(palette="Spectral") +
  ggtitle("Heme Transport") +
  labs(y = "Genes per million") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.6)) +
  theme(axis.title.y = element_text(size = 13L, face = "bold")) +
  theme(axis.title.x = element_text(size = 13L,
                                    face = "bold"))
Heme_T

Iron_T <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories.xlsx", 
                     sheet = "Iron transport")
Iron_T$sumTPM <- as.numeric(Iron_T$sumTPM)
Iron_T <- ggplot(Iron_T) +
  aes(x = Sample, y = sumTPM, fill = HMM) +
  geom_bar(position = "stack", stat = "identity", show.legend = TRUE,  width = 0.98) +
  scale_fill_viridis(discrete = T) +
  ggtitle("Iron Transport") +
  labs(y = "Genes per million") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.6)) +
  theme(axis.title.y = element_text(size = 13L, face = "bold")) +
  theme(axis.title.x = element_text(size = 13L,
                                    face = "bold"))
Iron_T

Sid_syn <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories.xlsx", 
                     sheet = "Siderophore synthesis")
Sid_syn$sumTPM <- as.numeric(Sid_syn$sumTPM)
Sid_syn$HMM <-gsub("-siderophore-synthesis","",as.character(Sid_syn$HMM))
Sid_syn <- ggplot(Sid_syn) +
  aes(x = Sample, y = sumTPM, fill = HMM) +
  geom_bar(position = "stack", stat = "identity", show.legend = TRUE,  width = 0.98) +
  scale_fill_viridis(discrete = T) +
  ggtitle("Siderophore Synthesis") +
  labs(y = "Genes per million") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.6)) +
  theme(axis.title.y = element_text(size = 13L, face = "bold")) +
  theme(axis.title.x = element_text(size = 13L,
                                    face = "bold"))
Sid_syn

  
  Sid_T <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories.xlsx", 
                        sheet = "Siderophore transport")
  Sid_T$sumTPM <- as.numeric(Sid_T$sumTPM)
  Sid_T$HMM <-gsub("-siderophore-export","",as.character(Sid_T$HMM))
  Sid_T$HMM <-gsub("-siderophore export","",as.character(Sid_T$HMM))
  Sid_T$HMM <-gsub("-sideropore-export","",as.character(Sid_T$HMM))
  Sid_T$HMM <-gsub("-siderophore-transport","",as.character(Sid_T$HMM))
  Sid_T$HMM <-gsub("-substrate-binding-protein","",as.character(Sid_T$HMM))
  Sid_T$HMM <-gsub("-siderophore-receptor","",as.character(Sid_T$HMM))
  Sid_T$HMM <-gsub("-siderophore receptor","",as.character(Sid_T$HMM))
  Sid_T$HMM <-gsub("-ATP-binding-protein","",as.character(Sid_T$HMM))
  Sid_T$HMM <-gsub("-siderophore-utilization","",as.character(Sid_T$HMM))
  Sid_T$HMM <-gsub("-permease","",as.character(Sid_T$HMM))
  Sid_T$HMM <-gsub("-iron-reductase","",as.character(Sid_T$HMM))
  Sid_T$HMM <-gsub("-hypothetical","",as.character(Sid_T$HMM))
  Sid_T <- ggplot(Sid_T) +
    aes(x = Sample, y = sumTPM, fill = HMM) +
    geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
    scale_fill_viridis(discrete = T) +
    ggtitle("Siderophore Transport") +
    labs(y = "Genes per million") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.6)) +
    theme(axis.title.y = element_text(size = 13L, face = "bold")) +
    theme(axis.title.x = element_text(size = 13L,
                                      face = "bold"))
  Sid_T

  Sid_TP <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories.xlsx", 
                      sheet = "Siderophore transport potential")
  Sid_TP$sumTPM <- as.numeric(Sid_TP$sumTPM)
  Sid_TP <- ggplot(Sid_TP) +
    aes(x = Sample, y = sumTPM, fill = HMM) +
    geom_bar(position = "stack", stat = "identity", show.legend = TRUE, width = 0.98) +
    scale_fill_brewer(palette = "Paired") +
  ggtitle("Siderophore Transport potential") +
    labs(y = "Genes per million") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.6)) +
    theme(axis.title.y = element_text(size = 13L, face = "bold")) +
    theme(axis.title.x = element_text(size = 13L,
                                      face = "bold"))
  Sid_TP
  
  Iron_GR <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories.xlsx", 
                       sheet = "Iron gene regulation")
  Iron_GR$sumTPM <- as.numeric(Iron_GR$sumTPM)
  Iron_GR <- ggplot(Iron_GR) +
    aes(x = Sample, y = sumTPM, fill = HMM) +
    geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
    scale_fill_brewer(palette = "Paired") +
    ggtitle("Iron gene regulation") +
    labs(y = "Genes per million") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.6)) +
    theme(axis.title.y = element_text(size = 13L, face = "bold")) +
    theme(axis.title.x = element_text(size = 13L,
                                      face = "bold"))
  Iron_GR

  Iron_S <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories.xlsx", 
                        sheet = "Iron storage")
  Iron_S$sumTPM <- as.numeric(Iron_S$sumTPM)
  Iron_S <- ggplot(Iron_S) +
    aes(x = Sample, y = sumTPM, fill = HMM) +
    geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
    scale_fill_brewer(palette = "Spectral") +
    ggtitle("Iron Storage") +
    labs(y = "Genes per million") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.6)) +
    theme(axis.title.y = element_text(size = 13L, face = "bold")) +
    theme(axis.title.x = element_text(size = 13L,
                                      face = "bold"))
  Iron_S
  
# all plots ###
  ## all plots in one
  plot_grid(Heme_T , Iron_T, Iron_GR, Iron_S, labels = c('A', 'B', 'C', 'D'), 
            label_size = 8)
  
  plot_grid(Sid_syn , Sid_T, Sid_TP, labels = c('A', 'B', 'C'), label_size = 8)
  
  
  ######  OR 
  
  figure <- ggarrange(Heme_T, Iron_T, Iron_S, Iron_GR,
                      ncol = 2, nrow = 2)
  figure1 <- ggarrange(Sid_syn, Sid_T, Sid_TP, 
                       ncol = 2, nrow = 2)
  
  
  
######### with sample 69 ##### 
  ##.n   69
  Heme_Ta <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_with69.xlsx", 
                       sheet = "Heme transport")
  Heme_Ta$sumTPM <- as.numeric(Heme_Ta$sumTPM)
  Heme_Ta$HMM <-gsub("-substrate-binding-protein","",as.character(Heme_Ta$HMM))
  Heme_Ta$HMM <-gsub("-ATP-binding-protein","",as.character(Heme_Ta$HMM))
  Heme_Ta <- ggplot(Heme_Ta) +
    aes(x = Sample, y = sumTPM, fill = HMM) +
    geom_bar(position = "stack", stat = "identity", show.legend = TRUE,  width = 0.98) +
    scale_fill_brewer(palette="Spectral") +
    ggtitle("Heme Transport") +
    labs(y = "Genes per million") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.6)) +
    theme(axis.title.y = element_text(size = 13L, face = "bold")) +
    theme(axis.title.x = element_text(size = 13L,
                                      face = "bold"))
  
  Iron_Ta <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_with69.xlsx", 
                       sheet = "Iron transport")
  Iron_Ta$sumTPM <- as.numeric(Iron_Ta$sumTPM)
  Iron_Ta <- ggplot(Iron_Ta) +
    aes(x = Sample, y = sumTPM, fill = HMM) +
    geom_bar(position = "stack", stat = "identity", show.legend = TRUE,  width = 0.98) +
    scale_fill_viridis(discrete = T) +
    ggtitle("Iron Transport") +
    labs(y = "Genes per million") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.6)) +
    theme(axis.title.y = element_text(size = 13L, face = "bold")) +
    theme(axis.title.x = element_text(size = 13L,
                                      face = "bold"))
  
  Sid_syn_a <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_with69.xlsx", 
                        sheet = "Siderophore synthesis")
  Sid_syn_a$sumTPM <- as.numeric(Sid_syn_a$sumTPM)
  Sid_syn_a <-  Sid_syn_a %>% mutate(sumTPM = replace(sumTPM, sumTPM<0.01, NA))
  Sid_syn_a <- drop_na(Sid_syn_a)
  Sid_syn_a$HMM <-gsub("-siderophore-synthesis","",as.character(Sid_syn_a$HMM))
  Sid_syn_a <- ggplot(Sid_syn_a) +
    aes(x = Sample, y = sumTPM, fill = HMM) +
    geom_bar(position = "stack", stat = "identity", show.legend = TRUE,  width = 0.98) +
    scale_fill_viridis(discrete = T) +
    ggtitle("Siderophore Synthesis") +
    labs(y = "Genes per million") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.6)) +
    theme(axis.title.y = element_text(size = 13L, face = "bold")) +
    theme(axis.title.x = element_text(size = 13L,
                                      face = "bold"))
  
  
  Sid_Ta <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_with69.xlsx", 
                      sheet = "Siderophore transport")
  Sid_Ta$sumTPM <- as.numeric(Sid_Ta$sumTPM)
  Sid_Ta$HMM <-gsub("-siderophore-export","",as.character(Sid_Ta$HMM))
  Sid_Ta$HMM <-gsub("-siderophore export","",as.character(Sid_Ta$HMM))
  Sid_Ta$HMM <-gsub("-sideropore-export","",as.character(Sid_Ta$HMM))
  Sid_Ta$HMM <-gsub("-siderophore-transport","",as.character(Sid_Ta$HMM))
  Sid_Ta$HMM <-gsub("-substrate-binding-protein","",as.character(Sid_Ta$HMM))
  Sid_Ta$HMM <-gsub("-siderophore-receptor","",as.character(Sid_Ta$HMM))
  Sid_Ta$HMM <-gsub("-siderophore receptor","",as.character(Sid_Ta$HMM))
  Sid_Ta$HMM <-gsub("-ATP-binding-protein","",as.character(Sid_Ta$HMM))
  Sid_Ta$HMM <-gsub("-siderophore-utilization","",as.character(Sid_Ta$HMM))
  Sid_Ta$HMM <-gsub("-permease","",as.character(Sid_Ta$HMM))
  Sid_Ta$HMM <-gsub("-iron-reductase","",as.character(Sid_Ta$HMM))
  Sid_Ta$HMM <-gsub("-hypothetical","",as.character(Sid_Ta$HMM))
  Sid_Ta <- ggplot(Sid_Ta) +
    aes(x = Sample, y = sumTPM, fill = HMM) +
    geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
    scale_fill_viridis(discrete = T) +
    ggtitle("Siderophore Transport") +
    labs(y = "Genes per million") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.6)) +
    theme(axis.title.y = element_text(size = 13L, face = "bold")) +
    theme(axis.title.x = element_text(size = 13L,
                                      face = "bold"))
  
  Sid_TPa <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_with69.xlsx", 
                       sheet = "Siderophore transport potential")
  Sid_TPa$sumTPM <- as.numeric(Sid_TPa$sumTPM)
  Sid_TPa <- ggplot(Sid_TPa) +
    aes(x = Sample, y = sumTPM, fill = HMM) +
    geom_bar(position = "stack", stat = "identity", show.legend = TRUE, width = 0.98) +
    scale_fill_brewer(palette="Spectral") +
    ggtitle("Siderophore Transport potential") +
    labs(y = "Genes per million") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.6)) +
    theme(axis.title.y = element_text(size = 13L, face = "bold")) +
    theme(axis.title.x = element_text(size = 13L,
                                      face = "bold"))
  
  Iron_GRa <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_with69.xlsx", 
                        sheet = "Iron gene regulation")
  Iron_GRa$sumTPM <- as.numeric(Iron_GRa$sumTPM)
  Iron_GRa <- ggplot(Iron_GRa) +
    aes(x = Sample, y = sumTPM, fill = HMM) +
    geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
    scale_fill_brewer(palette="Spectral") +
    ggtitle("Iron gene regulation") +
    labs(y = "Genes per million") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.6)) +
    theme(axis.title.y = element_text(size = 13L, face = "bold")) +
    theme(axis.title.x = element_text(size = 13L,
                                      face = "bold"))
  
  Iron_Sa <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_with69.xlsx", 
                       sheet = "Iron storage")
  Iron_Sa$sumTPM <- as.numeric(Iron_Sa$sumTPM)
  Iron_Sa <- ggplot(Iron_Sa) +
    aes(x = Sample, y = sumTPM, fill = HMM) +
    geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
    scale_fill_brewer(palette="Spectral") +
    ggtitle("Iron Storage") +
    labs(y = "Genes per million") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.6)) +
    theme(axis.title.y = element_text(size = 13L, face = "bold")) +
    theme(axis.title.x = element_text(size = 13L,
                                      face = "bold"))
  
  figure_A <- ggarrange(Heme_Ta, Iron_Ta, Iron_Sa,
                      ncol = 2, nrow = 2)
  figure_A
  
  figure_B <- ggarrange(Iron_GRa, Sid_TPa,
                        ncol = 2, nrow = 2)
  figure_B
  
  figure_C <- ggarrange(Sid_syn_a, Sid_Ta, 
                       ncol = 2, nrow = 2)
  figure_C
  
  
  
# Siderophore synthesis ##########  
 #### statistical analysis 
  Gene_abun_SS <- read_excel("Gene_Categories_with69.xlsx", 
           sheet = "Siderophore synthesis")
  Gene_abun_SS$sumTPM <- as.numeric(Gene_abun_SS$sumTPM)
  Gene_abun_SS <-  Gene_abun_SS %>% mutate(sumTPM = replace(sumTPM, sumTPM<0.01, NA))
  
  summary(Gene_abun_SS)
 summary(Gene_abun_SS$sumTPM) 
 Gene_abun_SS$sumTPM <- round(Gene_abun_SS$sumTPM, digits=2)
 
 
 library(naniar)
 
 vis_miss(Gene_abun_SS)
 Gene_abun_SS %>% drop_na() %>% vis_miss()
 
 boxplot(Gene_abun_SS$sumTPM)
 qqnorm(Gene_abun_SS$sumTPM)
 attach(Gene_abun_SS)
 tapply(sumTPM,Group,mean, na.rm=T)
 tapply(sumTPM,Group,sd, na.rm=T)
 tapply(sumTPM,Group,length)
 
 xbar <- tapply(sumTPM,Group,mean, na.rm=T)
 s <- tapply(sumTPM,Group,sd, na.rm=T)
 n <- tapply(sumTPM,Group,length)
 cbind(mean=xbar, std.dev=s, n=n)
 
 #    mean   std.dev  n.    NA=17
  #LacADN  3.102105  4.663154 78
 #LacAMP 15.470635 20.572130 78
 
 Gene_abun_SS <- drop_na(Gene_abun_SS)
 ##for LacADN
 par(mfrow=c(2,1))
 hist(Gene_abun_SS$sumTPM[Gene_abun_SS$Group=="LacADN"], col="gold", probability=TRUE,
      xlab="Genes per million", main="LacADN", breaks=seq(0.01,35))
 lines(density(Gene_abun_SS$sumTPM[Gene_abun_SS$Group=="LacADN"]), col="red", lwd=2)
 ###for LacAMP
 hist(Gene_abun_SS$sumTPM[Gene_abun_SS$Group=="LacAMP"], col="gold", probability=TRUE,
      xlab="Genes per million", main="LacAMP", breaks=seq(0,100))
 lines(density(Gene_abun_SS$sumTPM[Gene_abun_SS$Group=="LacAMP"]), col="red", lwd=2)
 
 ###Q2 Hypothesis testing
 ###Test the hypothesis, does LacAMP have higher gene abundance than LacADN for siderophore
 ## synthesis genes?
 # Hypothesis testing
 #Using both ANOVA and linear models to test this
 library(stats)
 anova_gene<-aov(sumTPM~Group, data=Gene_abun_SS)
 summary(anova_gene)
 lm_gene <-lm(sumTPM~Group,data=Gene_abun_SS)
 summary(lm_gene)
 
 # Estimate Std. Error t value Pr(>|t|)    
 # (Intercept)    3.102      1.636   1.896   0.0601
 ## GroupLacAMP   12.369      2.430   5.090 1.16e-06 ***
 
  Gene_abun_SS$Group <- as.factor(Gene_abun_SS$Group)
 contrasts(Gene_abun_SS$Group)

 ### Is it different based on sample type? 
 anova_gene_sa <-aov(sumTPM~Sample, data=Gene_abun_SS)
 summary(anova_gene_sa)
 lm_gene_sa <-lm(sumTPM~Sample,data=Gene_abun_SS)
 summary(lm_gene_sa)
 
 # Plotting the intercepts and slope pf the lm regression
 par(mfrow = c(1, 2))
 plot(Gene_abun_SS$sumTPM ~ Gene_abun_SS$Group, xlab="", ylab = "Genes per million")
 AdnAmp <- which(Gene_abun_SS$Group != "Amp")
 GroupM <- ifelse(Gene_abun_SS$Group == "Amp", 1, 0)
 plot(Gene_abun_SS$sumTPM[AdnAmp] ~ GroupM[AdnAmp], xlab = "Group", ylab = "Genes per million")
 abline(lm(Gene_abun_SS$sumTPM[AdnAmp] ~ GroupM[AdnAmp]))
 abline(v = 0, lty = "dotted")
 text(0.3, 15, "Intercept\n(GroupF)")
 
 # Run model without intercept
 model3<-lm(sumTPM~Group-1, data=Gene_abun_SS)
 summary(model3)
 
 #Coefficients:
  # Estimate Std. Error t value Pr(>|t|)    
 #GroupLacADN    3.102      1.636   1.896   0.0601 .  
 #GroupLacAMP   15.471      1.797   8.609 1.55e-14 ***
 
 
 ###Q3: Model fit diagnostics
 ##Plot model diagnostics to check residuals distribution and homocedasticity conditions/outliers
 par(mfrow=c(2,2))
 plot(lm_gene)
 plot(anova_gene)
 plot(model3)
 
 ##Check normality and equal variances (homocedasticity)
 ##Shapiro test
 shapiro.test(Gene_abun_SS$sumTPM[Gene_abun_SS$Group=="LacADN"])
 shapiro.test(Gene_abun_SS$sumTPM[Gene_abun_SS$Group=="LacAMP"])
 
 ##Barlett's test for equal variances
 bartlett.test(Gene_abun_SS$sumTPM~Gene_abun_SS$Group, data=Gene_abun_SS)
 
 #### t-test - group
 attach(Gene_abun_SS)
 Gene_abun_SS <- as.data.frame(Gene_abun_SS)
 
 Gene_abun_SS %>% t_test(sumTPM ~ Group,
                 var.equal = TRUE,
                 detailed = TRUE) %>%
   glimpse()
 
 # p           <dbl> 1.16e-06
 
 # Statistical test
 stat.test_1 <- Gene_abun_SS %>%
   t_test(sumTPM ~ Group) %>%
   add_significance()
 stat.test_1
 #.y.       group1 group2    n1    n2 statistic    df         p p.signif
 #  1 sumTPM LacADN LacAMP    76    63     -4.67  67.3 0.0000147 **** 
 
 compare_means(sumTPM ~ Group, data = Gene_abun_SS, method = "t.test")
 
# .y.    group1 group2          p    p.adj p.format p.signif method  
#   1 sumTPM LacADN LacAMP 0.0000147 0.000015 1.5e-05    ****     t-test
 
 p <- ggboxplot(rhb, x = "Group", y = "sumTPM",
                color = "Group", palette = "jco",
                add = "jitter")
 # Change method
 p + stat_compare_means(method = "t.test")
 
 #### statistical analysis based on subfamilies
 ###################### Based on each family 
 #### Rhb family 
 rhb <- read_excel("Siderophore synthesis .xlsx", 
                        sheet = "Rhb family")
 # Box plot facetted by "HMM"
 rhb1 <- ggboxplot(rhb, x = "Group", y = "sumTPM",
                   color = "Group", palette = "jco",
                   add = "jitter",
                   facet.by = "HMM", short.panel.labs = FALSE) +   theme(legend.position="none")
 # Use only p.format as label. Remove method name.
 rhb1 + stat_compare_means(method = "t.test",label.x = 1.5)
 rhb1 + stat_compare_means( label =  "p.signif", label.x = 1.5)
 
 #### Vab family 
 vab <- read_excel("Siderophore synthesis .xlsx", 
                   sheet = "Vab family")
 # Box plot facetted by "HMM"
 vab1 <- ggboxplot(vab, x = "Group", y = "sumTPM",
                   color = "Group", palette = "jco",
                   add = "jitter",
                   facet.by = "HMM", short.panel.labs = FALSE) +   theme(legend.position="none")
 # Use only p.format as label. Remove method name.
 vab1 + stat_compare_means(method = "t.test",label.x = 1.5)
 vab1 + stat_compare_means( label =  "p.signif", label.x = 1.5)

 #### Pvs family 
 pvs <- read_excel("Siderophore synthesis .xlsx", 
                   sheet = "Pvs family ")
 # Box plot facetted by "HMM"
 pvs1 <- ggboxplot(pvs, x = "Group", y = "sumTPM",
                   color = "Group", palette = "jco",
                   add = "jitter",
                   facet.by = "HMM", short.panel.labs = FALSE) +   theme(legend.position="none")
 # Use only p.format as label. Remove method name.
 pvs1 + stat_compare_means(method = "t.test",label.x = 1.5)
 pvs1 + stat_compare_means( label =  "p.signif", label.x = 1.5)
 
 #### Pvd family 
 pvd_ss <- read_excel("Siderophore synthesis .xlsx", 
                          sheet = "Pvd family ")
 # Box plot facetted by "HMM"
 pvd_ss1 <- ggboxplot(pvd_ss, x = "Group", y = "sumTPM",
                   color = "Group", palette = "jco",
                   add = "jitter",
                   facet.by = "HMM", short.panel.labs = FALSE) +   theme(legend.position="none")
 pvd_ss1 + stat_compare_means(method = "wilcox.test",label.x = 1.5)
 pvd_ss1 + stat_compare_means( label =  "p.signif", label.x = 1.5)
 
 compare_means(sumTPM ~ Group, data = pvd_ss, method = "t.test")
 p1 <- ggboxplot(pvd_ss, x = "Group", y = "sumTPM",
                color = "Group", palette = "jco",
                add = "jitter")
 p1 + stat_compare_means(method = "t.test")
 ## p- value = 0.0000935
 
 
 
 
 
 
 
 ## Siderophore transport ##### 
 Gene_abun_ST <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_with69.xlsx", 
                                          sheet = "Siderophore transport")
 Gene_abun_ST$sumTPM <- as.numeric(Gene_abun_ST$sumTPM)
 Gene_abun_ST <-  Gene_abun_ST %>% mutate(sumTPM = replace(sumTPM, sumTPM<0.01, NA))
 
 Gene_abun_ST$sumTPM <- round(Gene_abun_ST$sumTPM, digits=2)
 library(naniar)
 
 vis_miss(Gene_abun_ST)
 Gene_abun_ST %>% drop_na() %>% vis_miss()
 
 boxplot(Gene_abun_ST$sumTPM)
 qqnorm(Gene_abun_ST$sumTPM)
 attach(Gene_abun_ST)
 tapply(sumTPM,Group,mean, na.rm=T)
 tapply(sumTPM,Group,sd, na.rm=T)
 tapply(sumTPM,Group,length)
 
 xbar <- tapply(sumTPM,Group,mean, na.rm=T)
 s <- tapply(sumTPM,Group,sd, na.rm=T)
 n <- tapply(sumTPM,Group,length)
 cbind(mean=xbar, std.dev=s, n=n)
 
# mean  std.dev  n
 #LacADN  67.13013 109.5218 78
 #LacAMP 106.29545 138.6912 78
 Gene_abun_ST <- drop_na(Gene_abun_ST)
 ##for LacADN
 par(mfrow=c(2,1))
 hist(Gene_abun_ST$sumTPM[Gene_abun_ST$Group=="LacADN"], col="gold", probability=TRUE,
      xlab="Genes per million", main="LacADN", breaks=seq(0, 500))
 lines(density(Gene_abun_ST$sumTPM[Gene_abun_ST$Group=="LacADN"]), col="red", lwd=2)
 ###for LacAMP
 hist(Gene_abun_ST$sumTPM[Gene_abun_ST$Group=="LacAMP"], col="gold", probability=TRUE,
      xlab="Genes per million", main="LacAMP", breaks=seq(0,700))
 lines(density(Gene_abun_ST$sumTPM[Gene_abun_ST$Group=="LacAMP"]), col="red", lwd=2)
 
 ###Q2 Hypothesis testing
 ###Test the hypothesis, does LacAMP have higher gene abundance than LacADN for siderophore
 ## transport genes?
 # Hypothesis testing
 #Using both ANOVA and linear models to test this
 library(stats)
 anova_gene_ST <-aov(sumTPM~Group, data=Gene_abun_ST)
 summary(anova_gene_ST)
 lm_gene_ST <-lm(sumTPM~Group,data=Gene_abun_ST)
 summary(lm_gene_ST)
 
 #Estimate Std. Error t value Pr(>|t|)    
 #(Intercept)    67.13      14.14   4.748 4.68e-06 ***
  # GroupLacAMP    39.17      20.06   1.952   0.0527 .  

 
 Gene_abun_ST$Group <- as.factor(Gene_abun_ST$Group)
 contrasts(Gene_abun_SS$Group)
 
 ### Is it different based on sample type? 
 anova_gene_sa_ST <-aov(sumTPM~Sample, data=Gene_abun_ST)
 summary(anova_gene_sa_ST)
 lm_gene_sa_ST <-lm(sumTPM~Sample,data=Gene_abun_ST)
 summary(lm_gene_sa_ST)
 
 # Plotting the intercepts and slope pf the lm regression
 par(mfrow = c(1, 2))
 plot(Gene_abun_ST$sumTPM ~ Gene_abun_ST$Group, xlab="", ylab = "Genes per million")
 AdnAmp_ST <- which(Gene_abun_ST$Group != "Amp")
 GroupM_ST <- ifelse(Gene_abun_ST$Group == "Amp", 1, 0)
 plot(Gene_abun_ST$sumTPM[AdnAmp] ~ GroupM[AdnAmp], xlab = "Group", ylab = "Genes per million")
 abline(lm(Gene_abun_ST$sumTPM[AdnAmp] ~ GroupM[AdnAmp]))
 abline(v = 0, lty = "dotted")
 text(0.3, 15, "Intercept\n(GroupM_ST)")
 
 # Run model without intercept
 model3_ST <-lm(sumTPM~Group-1, data=Gene_abun_ST)
 summary(model3_ST)
 
 #Coefficients:
 #  Estimate Std. Error t value Pr(>|t|)    
 #GroupLacADN    67.13      14.14   4.748 4.68e-06 ***
 #  GroupLacAMP   106.30      14.23   7.470 5.72e-12 ***
 
 ###Q3: Model fit diagnostics
 ##Plot model diagnostics to check residuals distribution and homocedasticity conditions/outliers
 par(mfrow=c(2,2))
 plot(lm_gene_ST)
 plot(anova_gene_ST)
 plot(model3_ST)
 
 ##Check normality and equal variances (homocedasticity)
 ##Shapiro test
 shapiro.test(Gene_abun_ST$sumTPM[Gene_abun_ST$Group=="LacADN"])
 ### p-value = 3.042e-12
 shapiro.test(Gene_abun_ST$sumTPM[Gene_abun_ST$Group=="LacAMP"])
 ## p-value = 5.471e-11
 
 ##Barlett's test for equal variances
 bartlett.test(Gene_abun_ST$sumTPM~Gene_abun_ST$Group, data=Gene_abun_ST)
 # p-value = 0.04036
 
 #### t-test - group
 attach(Gene_abun_ST)
 Gene_abun_ST <- as.data.frame(Gene_abun_ST)
 
 Gene_abun_ST %>% t_test(sumTPM ~ Group,
                         var.equal = TRUE,
                         detailed = TRUE) %>%
   glimpse()
 
 # p           <dbl> 1.16e-06
 
 # Statistical test
 stat.test_2 <- Gene_abun_ST %>%
   t_test(sumTPM ~ Group) %>%
   add_significance()
 stat.test_2
 #.y.       group1 group2    n1    n2 statistic    df         p p.signif
 #  1 sumTPM LacADN LacAMP    76    63     -4.67  67.3 0.0000147 **** 
 
 compare_means(sumTPM ~ Group, data = Gene_abun_ST, method = "t.test")
 ## P -VALUE = 0.053
 compare_means(sumTPM ~ Group, data = Gene_abun_ST)
 ### P-VALUE = 0.0006
 
 # .y.    group1 group2          p    p.adj p.format p.signif method  
 #   1 sumTPM LacADN LacAMP 0.0000147 0.000015 1.5e-05    ****     t-test
 
 p_ST <- ggboxplot(Gene_abun_ST, x = "Group", y = "sumTPM",
                color = "Group", palette = "jco",
                add = "jitter")
 # Change method
 p_ST + stat_compare_means(method = "wilcox.test")
 
 
 ############ Based on family 
 Sid_trans_fam <- read_excel("Siderophore transport.xlsx", 
                                         sheet = "Clean")
 Sid_trans_fam$HMM <-gsub("-siderophore-export","",as.character(Sid_trans_fam$HMM))
 Sid_trans_fam$HMM <-gsub("-siderophore export","",as.character(Sid_trans_fam$HMM))
 Sid_trans_fam$HMM <-gsub("-sideropore-export","",as.character(Sid_trans_fam$HMM))
 Sid_trans_fam$HMM <-gsub("-siderophore-transport","",as.character(Sid_trans_fam$HMM))
 Sid_trans_fam$HMM <-gsub("-substrate-binding-protein","",as.character(Sid_trans_fam$HMM))
 Sid_trans_fam$HMM <-gsub("-siderophore-receptor","",as.character(Sid_trans_fam$HMM))
 Sid_trans_fam$HMM <-gsub("-siderophore receptor","",as.character(Sid_trans_fam$HMM))
 Sid_trans_fam$HMM <-gsub("-ATP-binding-protein","",as.character(Sid_trans_fam$HMM))
 Sid_trans_fam$HMM <-gsub("-siderophore-utilization","",as.character(Sid_trans_fam$HMM))
 Sid_trans_fam$HMM <-gsub("-permease","",as.character(Sid_trans_fam$HMM))
 Sid_trans_fam$HMM <-gsub("-iron-reductase","",as.character(Sid_trans_fam$HMM))
 Sid_trans_fam$HMM <-gsub("-hypothetical","",as.character(Sid_trans_fam$HMM))

 #Sid_trans_fam <- as.data.frame(Sid_trans_fam)
 #Sid_trans_fam$Sub_family <- paste(Sid_trans_fam$Group, Sid_trans_fam$HMM, sep="_") 
 
 #####Sub_family ######################
 #### Fpt family 
 fpt <- read_excel("Siderophore transport.xlsx", 
                                          sheet = "Fpt family")
 # Box plot facetted by "HMM"
 fpt1 <- ggboxplot(fpt, x = "Group", y = "sumTPM",
                color = "Group", palette = "jco",
                add = "jitter",
                facet.by = "HMM", short.panel.labs = FALSE) +   theme(legend.position="none")
 # Use only p.format as label. Remove method name.
 fpt1 + stat_compare_means(method = "t.test",label.x = 1.5)
 fpt1 + stat_compare_means( label =  "p.signif", label.x = 1.5)
 
 #### Fpv family 
 fpv <- read_excel("Siderophore transport.xlsx", 
                   sheet = "Fpv family")
 # Box plot facetted by "HMM"
 fpv1 <- ggboxplot(fpv, x = "Group", y = "sumTPM",
                   color = "Group", palette = "jco",
                   add = "jitter", 
                   facet.by = "HMM", short.panel.labs = TRUE) +   theme(legend.position="none")
 # Use only p.format as label. Remove method name.
 fpv1 + stat_compare_means(method = "wilcox.test", label.x = 1.5)
 fpv1 + stat_compare_means (label =  "p.signif", label.x = 1.5)
 
 #### Pvd family 
 pvd <- read_excel("Siderophore transport.xlsx", 
                   sheet = "Pvd family")
 # Box plot facetted by "HMM"
 pvd1 <- ggboxplot(pvd, x = "Group", y = "sumTPM",
                   color = "Group", palette = "jco",
                   add = "jitter", 
                   facet.by = "HMM", short.panel.labs = FALSE) +   theme(legend.position="none")
 # Use only p.format as label. Remove method name.
 pvd1 + stat_compare_means(method = "t.test", label.x = 1.5)
 pvd1 + stat_compare_means (label =  "p.signif", label.x = 1.5)
 
 #### Pvu family 
 pvu <- read_excel("Siderophore transport.xlsx", 
                   sheet = "Pvu family")
 # Box plot facetted by "HMM"
 pvu1 <- ggboxplot(pvu, x = "Group", y = "sumTPM",
                   color = "Group", palette = "jco",
                   add = "jitter", 
                   facet.by = "HMM", short.panel.labs = TRUE) +   theme(legend.position="none")
 # Use only p.format as label. Remove method name.
 pvu1 + stat_compare_means(method = "t.test", label.x = 1.5)
 pvu1 + stat_compare_means (label =  "p.signif", label.x = 1.5)
 
 #### Others family 
 others_st <- read_excel("Siderophore transport.xlsx", 
                   sheet = "Others")
 # Box plot facetted by "HMM"
 others_st_1 <- ggboxplot(others_st, x = "Group", y = "sumTPM",
                   color = "Group", palette = "jco",
                   add = "jitter", 
                   facet.by = "HMM", short.panel.labs = TRUE) +   theme(legend.position="none")
 # Use only p.format as label. Remove method name.
 others_st_1 + stat_compare_means(method = "t.test", label.x = 1.5)
 others_st_1 + stat_compare_means (label =  "p.signif", label.x = 1.5)
 
 
### Family level #######
 #### Siderophore synthesis genes
 ## Pvd family 
 pvd_ss <- read_excel("Siderophore synthesis.xlsx", 
                       sheet = "Pvd family ")
 p1 <- ggboxplot(pvd_ss, x = "Group", y = "sumTPM",
                 color = "HMM", palette = "jco",
                 add = "jitter")
 p2 <- p1 + stat_compare_means(method = "t.test")
 
 ## Pch family 
 pch_ss <- read_excel("Siderophore synthesis .xlsx", 
                          sheet = "Pch family")
 pc1 <- ggboxplot(pch_ss, x = "Group", y = "sumTPM",
                 color = "HMM", palette = "jco",
                 add = "jitter")
 pc1 + stat_compare_means(method = "t.test")
 
 
 ## Luc family 
 luc_ss <- read_excel("Siderophore synthesis .xlsx", 
                          sheet = "luc Family ")
 l1 <- ggboxplot(luc_ss, x = "Group", y = "sumTPM",
                  color = "HMM", palette = "jco",
                  add = "jitter")
 l1 + stat_compare_means(method = "t.test")
 
 
 ## Pvs family 
 pvs_ss <- read_excel("Siderophore synthesis .xlsx", 
                      sheet = "Pvs family ")
 pvs1 <- ggboxplot(pvs_ss, x = "Group", y = "sumTPM",
                 color = "HMM", palette = "jco",
                 add = "jitter")
 pvs2 <- pvs1 + stat_compare_means(method = "t.test")
 
 ## Vab family 
 vab_ss <- read_excel("Siderophore synthesis .xlsx", 
                      sheet = "Vab family ")
 v1 <- ggboxplot(vab_ss, x = "Group", y = "sumTPM",
                 color = "HMM", palette = "jco",
                 add = "jitter")
 v1 + stat_compare_means(method = "t.test")
 
 ## Pvd and PvS family 
pvd_pvs_fam <- ggarrange(p1 , pvs1,
                         ncol=2, nrow =1)
pvd_pvs_fam

 ##############  Siderophore transport genes
 ### Fpt family
 fpt_st <- read_excel("Siderophore transport.xlsx", 
                          sheet = "Fpt family")
 f1 <- ggboxplot(fpt_st, x = "Group", y = "sumTPM",
                 color = "HMM", palette = "jco",
                 add = "jitter")
 f1 + stat_compare_means(method = "t.test")
 
 ### Fpv family
 fpv_st <- read_excel("Siderophore transport.xlsx", 
                      sheet = "Fpv family")
 fpv1 <- ggboxplot(fpv_st, x = "Group", y = "sumTPM",
                 color = "HMM", palette = "jco",
                 add = "jitter")
 fpv1 + stat_compare_means(method = "t.test")
 
 ### Pvd family
 pvd_st <- read_excel("Siderophore transport.xlsx", 
                      sheet = "Pvd family")
 pvd_st1 <- ggboxplot(pvd_st, x = "Group", y = "sumTPM",
                   color = "HMM", palette = "jco",
                   add = "jitter")
 pvd_st1 + stat_compare_means(method = "t.test")

 ### Pvu family
 pvu_st <- read_excel("Siderophore transport.xlsx", 
                      sheet = "Pvu family")
 pvu1 <- ggboxplot(pvu_st, x = "Group", y = "sumTPM",
                   color = "HMM", palette = "jco",
                   add = "jitter")
 pvu1 + stat_compare_means(method = "t.test") 
 
 ### Plots - Gene categories based on type of samples ####
 heme_t <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_Samples.xlsx", 
                       sheet = "Heme transport")
 heme_t$sumTPM <- as.numeric(heme_t$sumTPM)
 heme_t$HMM <-gsub("-substrate-binding-protein","",as.character(heme_t$HMM))
 heme_t$HMM <-gsub("-ATP-binding-protein","",as.character(heme_t$HMM))
 heme_t <- ggplot(heme_t) +
   aes(x = Sample, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE,  width = 0.98) +
   scale_fill_brewer(palette="Spectral") +
   ggtitle("Heme Transport") +
   labs(y = "Genes per million") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 heme_t
 
 iron_t <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_Samples.xlsx", 
                       sheet = "Iron transport")
 iron_t$sumTPM <- as.numeric(iron_t$sumTPM)
 iron_t <- ggplot(iron_t) +
   aes(x = Sample, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE,  width = 0.98) +
   scale_fill_viridis(discrete = T) +
   ggtitle("Iron Transport") +
   labs(y = "Genes per million") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 iron_t
 
 sid_tranp_a <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_Samples.xlsx", 
                       sheet = "Siderophore transport potential")
 sid_tranp_a$sumTPM <- as.numeric(sid_tranp_a$sumTPM)
 sid_tranp_a <- ggplot(sid_tranp_a) +
   aes(x = Sample, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE, width = 0.98) +
   scale_fill_brewer(palette="Spectral") +
   ggtitle("Siderophore Transport potential") +
   labs(y = "Genes per million") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 sid_tranp_a
 
 iron_gener_a <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_Samples.xlsx", 
                        sheet = "Iron gene regulation")
 iron_gener_a$sumTPM <- as.numeric(iron_gener_a$sumTPM)
 iron_gener_a <- ggplot(iron_gener_a) +
   aes(x = Sample, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
   scale_fill_brewer(palette="Spectral") +
   ggtitle("Iron gene regulation") +
   labs(y = "Genes per million") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 iron_gener_a
 
 iron_sa <- read_excel("~/Desktop/BINGO2/Co-assembly/Gene_Categories_Samples.xlsx", 
                       sheet = "Iron storage")
 iron_sa$sumTPM <- as.numeric(iron_sa$sumTPM)
 iron_sa <- ggplot(iron_sa) +
   aes(x = Sample, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
   scale_fill_brewer(palette="Spectral") +
   ggtitle("Iron Storage") +
   labs(y = "Genes per million") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 iron_sa
 
 figure_A <- ggarrange(heme_t, iron_t, iron_sa, 
                       ncol = 2, nrow = 2)
 figure_A
 figure_B <- ggarrange(Sid_syn_a, Sid_Ta, Sid_TPa, 
                       ncol = 2, nrow = 2)
 figure_B
 
 
 ### Plots - Gene categories based on type of samples - REPLICATES ####
 c1 <- read_excel("Gene_Categories_with69.xlsx", 
                      sheet = "Heme transport")
 c1 <- c1 %>%
   reorder_levels("Sample", order = c("LacAMPa", "LacAMPb", "LacAMPc",
                                      "LacADNa", "LacADNb", "LacADNc"))
 levels(c1$Sample)
 
 c1$sumTPM <- as.numeric(c1$sumTPM)
 c1$HMM <-gsub("-substrate-binding-protein","",as.character(c1$HMM))
 c1$HMM <-gsub("-ATP-binding-protein","",as.character(c1$HMM))
 c1 <- ggplot(c1) +
   aes(x = Sample, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE,  width = 0.98) +
   scale_fill_brewer(palette="Spectral") +
   ggtitle("Heme Transport") +
   labs(y = "Genes per million") +
   theme_minimal() +
 theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 c1
 
 c2 <- read_excel("Gene_Categories_with69.xlsx", 
                      sheet = "Iron transport")
 c2 <- c2 %>%
   reorder_levels("Sample", order = c("LacAMPa", "LacAMPb", "LacAMPc",
                                      "LacADNa", "LacADNb", "LacADNc"))
 levels(c2$Sample)
 c2$sumTPM <- as.numeric(c2$sumTPM)
 c2 <- ggplot(c2) +
   aes(x = Sample, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE,  width = 0.98) +
   scale_fill_viridis(discrete = T) +
   ggtitle("Iron Transport") +
   labs(y = "Genes per million") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 c2
 
 c3 <- read_excel("Gene_Categories_with69.xlsx",
                            sheet = "Iron gene regulation")
 c3$sumTPM <- as.numeric(c3$sumTPM)
 c3 <- c3 %>%
   reorder_levels("Sample", order = c("LacAMPa", "LacAMPb", "LacAMPc",
                                      "LacADNa", "LacADNb", "LacADNc"))
 levels(c3$Sample)
 c3 <- ggplot(c3) +
   aes(x = Sample, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
   scale_fill_brewer(palette="Spectral") +
   ggtitle("Iron gene regulation") +
   labs(y = "Genes per million") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 c3
 
 c4 <- read_excel("Gene_Categories_with69.xlsx",
                       sheet = "Iron storage")
 c4$sumTPM <- as.numeric(c4$sumTPM)
 c4 <- c4 %>%
   reorder_levels("Sample", order = c("LacAMPa", "LacAMPb", "LacAMPc",
                                      "LacADNa", "LacADNb", "LacADNc"))
 levels(c4$Sample)
 c4 <- ggplot(c4) +
   aes(x = Sample, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
   scale_fill_brewer(palette="Spectral") +
   ggtitle("Iron Storage") +
   labs(y = "Genes per million") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 
 c4
 
 d1 <- read_excel("Gene_Categories_with69.xlsx",
                           sheet = "Siderophore synthesis")
 d1$sumTPM <- as.numeric(d1$sumTPM)
 d1$HMM <-gsub("-siderophore-synthesis","",as.character(d1$HMM))
 d1 <- d1 %>%
   reorder_levels("Sample", order = c("LacAMPa", "LacAMPb", "LacAMPc",
                                      "LacADNa", "LacADNb", "LacADNc"))
 levels(d1$Sample)
 d1 <- ggplot(d1) +
   aes(x = Sample, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE, width = 0.98) +
   scale_fill_viridis(discrete = T) +
   ggtitle("Siderophore synthesis") +
   labs(y = "Genes per million") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 d1
 
 d2 <- read_excel("Gene_Categories_with69.xlsx",
                           sheet = "Siderophore transport")
 d2$sumTPM <- as.numeric(d2$sumTPM)
 d2$HMM <-gsub("-siderophore-export","",as.character(d2$HMM))
 d2$HMM <-gsub("-siderophore export","",as.character(d2$HMM))
 d2$HMM <-gsub("-sideropore-export","",as.character(d2$HMM))
 d2$HMM <-gsub("-siderophore-transport","",as.character(d2$HMM))
 d2$HMM <-gsub("-substrate-binding-protein","",as.character(d2$HMM))
 d2$HMM <-gsub("-siderophore-receptor","",as.character(d2$HMM))
 d2$HMM <-gsub("-siderophore receptor","",as.character(d2$HMM))
 d2$HMM <-gsub("-ATP-binding-protein","",as.character(d2$HMM))
 d2$HMM <-gsub("-siderophore-utilization","",as.character(d2$HMM))
 d2$HMM <-gsub("-permease","",as.character(d2$HMM))
 d2$HMM <-gsub("-iron-reductase","",as.character(d2$HMM))
 d2$HMM <-gsub("-hypothetical","",as.character(d2$HMM))
 d2 <- d2 %>%
   reorder_levels("Sample", order = c("LacAMPa", "LacAMPb", "LacAMPc",
                                      "LacADNa", "LacADNb", "LacADNc"))
 levels(d2$Sample)
 d2 <- ggplot(d2) +
   aes(x = Sample, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE, width = 0.98) +
   scale_fill_viridis(discrete = T) +
   ggtitle("Siderophore Transport") +
   labs(y = "Genes per million") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 d2
 
 d3 <- read_excel("Gene_Categories_with69.xlsx",
                           sheet = "Siderophore transport potential")
 d3$sumTPM <- as.numeric(d3$sumTPM)
 d3 <- d3 %>%
   reorder_levels("Sample", order = c("LacAMPa", "LacAMPb", "LacAMPc",
                                      "LacADNa", "LacADNb", "LacADNc"))
 levels(d3$Sample)
 d3 <- ggplot(d3) +
   aes(x = Sample, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE, width = 0.98) +
   scale_fill_brewer(palette="Spectral") +
   ggtitle("Siderophore Transport potential") +
   labs(y = "Genes per million") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 d3
 
 
 figure_C <- ggarrange(c1, c2, c3, c4 ,
                       ncol = 2, nrow = 2)
 figure_C
 figure_D <- ggarrange(d1, d2, d3, 
                       ncol = 2, nrow = 2)
 figure_D
 
 
 
 
 library(ggstatsplot)
 ## Plots with p- values for other gene families ####
 ## Heme transport ##
 ht <- read_excel("Gene_Categories_coassembly.xlsx", 
                        sheet = "Heme transport")
 ht$sumTPM <- as.numeric(ht$sumTPM)
 hist(ht$sumTPM)
 
 model_ht <- aov(sumTPM~Group, data=ht)
 summary(model_ht)
 
 ht$Group <- as.factor(ht$Group)
 # Box plot facetted by "HMM"
 ht_1 <- ggboxplot(ht, x = "Group", y = "sumTPM", 
                      order = c("LacAMP", "LacADN"), 
                      color = "Group", 
                      palette = "lancet", 
                      add = "jitter", 
                      title = "Heme Transport", 
                      label.rectangle = TRUE, 
                      ggtheme = theme_bw(), 
                      facet.by = "HMM", short.panel.labs = TRUE)
 # Use only p.format as label. Remove method name.
 ht_1 + stat_compare_means(label = "p.signif", label.x = 1.5)
 
 #ht_1 <- ggboxplot(ht, x = "Group", y = "sumTPM",
                  #color = "HMM", palette = "jco",
                  #add = "jitter")
 #ht_1 + stat_compare_means(method = "t.test", label.x = 1.5)
 #ht_1 + stat_compare_means(label = "p.signif", label.x = 1.5)
 

 ### Iron transport ##
 it <- read_excel("Gene_Categories_coassembly.xlsx", 
                      sheet = "Iron Transport")
 it$sumTPM <- as.numeric(it$sumTPM)
 # Box plot facetted by "HMM"
 it_1 <- ggboxplot(it, x = "Group", y = "sumTPM", 
                   order = c("LacAMP", "LacADN"), 
                   color = "Group", 
                   palette = "lancet", 
                   add = "jitter", 
                   title = "Iron Transport", 
                   label.rectangle = TRUE, 
                   ggtheme = theme_bw(), 
                   facet.by = "HMM", short.panel.labs = TRUE)
 # Use only p.format as label. Remove method name.
 it_1 + stat_compare_means(label = "p.signif", label.x = 1.5)

 
 ### Iron gene regulation  ##
 igr <- read_excel("Gene_Categories_coassembly.xlsx", 
                  sheet = "Iron gene regulation")
 igr$sumTPM <- as.numeric(igr$sumTPM)
 # Box plot facetted by "HMM"
 igr_1 <- ggboxplot(igr, x = "Group", y = "sumTPM", 
                   order = c("LacAMP", "LacADN"), 
                   color = "Group", 
                   palette = "lancet", 
                   add = "jitter", 
                   title = "Iron Gene Regulation", 
                   label.rectangle = TRUE, 
                   ggtheme = theme_bw(), 
                   facet.by = "HMM", short.panel.labs = TRUE)
 # Use only p.format as label. Remove method name.
 igr_1 + stat_compare_means(label = "p.signif", label.x = 1.5)
 
 ### Iron storage ##
 ist <- read_excel("Gene_Categories_coassembly.xlsx", 
                   sheet = "Iron storage")
 ist$sumTPM <- as.numeric(ist$sumTPM)
 # Box plot facetted by "HMM"
 ist_1 <- ggboxplot(ist, x = "Group", y = "sumTPM", 
                    order = c("LacAMP", "LacADN"), 
                    color = "Group", 
                    palette = "lancet", 
                    add = "jitter", 
                    title = "Iron Gene Regulation", 
                    label.rectangle = TRUE, 
                    ggtheme = theme_bw(), 
                    facet.by = "HMM", short.panel.labs = TRUE)
 # Use only p.format as label. Remove method name.
 ist_1 + stat_compare_means(label = "p.signif", label.x = 1.5)
 
 ### Siderophore transport potential ##
 sid_pot <- read_excel("Gene_Categories_coassembly.xlsx", 
                   sheet = "Siderophore transport potential")
 sid_pot$sumTPM <- as.numeric(sid_pot$sumTPM)
 # Box plot facetted by "HMM"
 sid_pot_1 <- ggboxplot(sid_pot, x = "Group", y = "sumTPM", 
                    order = c("LacAMP", "LacADN"), 
                    color = "Group", 
                    palette = "lancet", 
                    add = "jitter", 
                    title = "Siderophore Transport Potential", 
                    label.rectangle = TRUE, 
                    ggtheme = theme_bw(), 
                    facet.by = "HMM", short.panel.labs = TRUE)
 # Use only p.format as label. Remove method name.
 sid_pot_1 + stat_compare_means(label = "p.signif", label.x = 1.5)
 
 

 #### Combination of gene families - Siderophore ##
 ## Siderophore synthesis - some genes ####
 pvd_pvs <- read_excel("Siderophore synthesis.xlsx", 
                       sheet = "Pvd_Pvs")
 pvd_pvs$Group <- as.factor(pvd_pvs$Group)
 pvd_pvs$Group <- factor(pvd_pvs$Group, levels=c('LacAMP', 'LacADN'))
 pvd_pvs_1 <- ggboxplot(pvd_pvs, x = "Sample", y = "sumTPM",
                        color = "HMM", palette = "jco",
                        add = "jitter")
 pvd_pvs_1 + stat_compare_means(method = "t.test")
 
 # Box plot facetted by "HMM"
 tiff("test.tiff", units="in", width=5, height=5, res=300)
 sid_pvd_pvs <- ggboxplot(pvd_pvs, x = "Group", y = "sumTPM",
                          color = "Group",
                          add = "jitter",
                          facet.by = "HMM", short.panel.labs = FALSE) +  
   scale_fill_brewer(palette="BuPu") +
   theme(legend.position="none")
 # Use only p.format as label. Remove method name.
 sid_pvd_pvs + stat_compare_means(method = "t.test",label.x = 1.5)
 sid_pvd_pvs + stat_compare_means( label =  "p.signif", label.x = 1.5) + 
   guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
   ylab("GPM (Genes per kilobase million)") 
 dev.off()
 
 
 tiff("pvd_pvs.tiff", units="in", width=5, height=5, res=300)
 sid_pvd_pvs_1 <- ggplot(pvd_pvs, aes(x=Group, y=sumTPM, fill=Group)) + 
   geom_boxplot(outlier.shape = NA) +
   coord_cartesian(ylim = c(0, 70))+
   facet_wrap(~HMM, scale="free_x") +
   stat_compare_means(paired = TRUE, label =  "p.format", label.x = 1.5) + 
   #guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
   ylab("GPM (Genes per kilobase million)") +
   scale_fill_brewer(palette="Paired") +
   # change order of legend
   # Add a theme
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   theme(axis.text.x=element_blank(),
         axis.ticks.x=element_blank() ) 
 sid_pvd_pvs_1
 dev.off()
 
 ## Siderophore transport ####
 fpv_pvu <- read_excel("Siderophore transport.xlsx", 
                  sheet = "Fpv_Pvu")
 fpv_pvu$sumTPM <- as.numeric(fpv_pvu$sumTPM)
 fpv_pvu$HMM <-gsub(" - siderophore transport","",as.character(fpv_pvu$HMM))
 fpv_pvu$Group <- as.factor(fpv_pvu$Group)
 fpv_pvu$Group <- factor(fpv_pvu$Group, levels=c('LacAMP', 'LacADN'))

 # Box plot facetted by "HMM"
 sid_fpv_pvu <- ggboxplot(fpv_pvu, x = "Group", y = "sumTPM",
                          color = "Group",
                          add = "jitter",
                          facet.by = "HMM", short.panel.labs = FALSE) +  
   scale_fill_brewer(palette="BuPu") +
   theme(legend.position="none")
 # Use only p.format as label. Remove method name.
 sid_fpv_pvu + stat_compare_means(method = "t.test",label.x = 1.5)
 sid_fpv_pvu + stat_compare_means( label =  "p.signif", label.x = 1.5) + 
   guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
   ylab("GPM (Genes per kilobase million)") 
 dev.off()
 
 
 tiff("fpv_pvu.tiff", units="in", width=5, height=5, res=300)
 sid_fpv_pvu_1 <-ggplot(fpv_pvu, aes(x=Group, y=sumTPM, fill=Group)) + 
   geom_boxplot(outlier.shape = NA) +
   coord_cartesian(ylim = c(0, 300))+
   facet_wrap(~HMM, scale="free_x") +
   stat_compare_means(paired = TRUE, label =  "p.signif", label.x = 1.5, label.y = 270) + 
   #guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
   ylab("GPM (Genes per kilobase million)") +
   scale_fill_brewer(palette="Paired") +
   # change order of legend
   # Add a theme
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   theme(axis.text.x=element_blank(),
         axis.ticks.x=element_blank())
 sid_fpv_pvu_1
 dev.off()
 
 library(patchwork)
 combined <- sid_pvd_pvs_1 + sid_fpv_pvu_1 & theme(legend.position = "bottom")
 combined + plot_layout(guides = "collect")
 
 
## pyoverdine and vibrioferrin #####
 all_1 <- read_excel("Siderophore transport.xlsx", 
                       sheet = "Combined")
 all_1$sumTPM <- as.numeric(all_1$sumTPM)
 all_1$Group <- as.factor(all_1$Group)
 all_1$Group <- factor(all_1$Group, levels=c('LacAMP', 'LacADN'))
 
 # Box plot facetted by "HMM"
 all_2 <- ggboxplot(all_1, x = "Group", y = "sumTPM",
                          color = "Group",
                          add = "jitter",
                          facet.by = "HMM", short.panel.labs = FALSE) +  
   scale_fill_brewer(palette="BuPu") +
   theme(legend.position="none")
 # Use only p.format as label. Remove method name.
 all_2 + stat_compare_means(method = "t.test",label.x = 1.5)
 all_2 + stat_compare_means( label =  "p.signif", label.x = 1.5) + 
   guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
   ylab("GPM (Genes per kilobase million)") 
 dev.off()
 

 all_3 <-ggplot(all_1, aes(x=Group, y=sumTPM, fill=Group)) + 
   geom_boxplot(outlier.shape = NA) +
   facet_wrap(~HMM, scale="free_x") +
   stat_compare_means(paired = TRUE, label =  "p.signif", label.x = 1.5) + 
   #guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
   ylab("GPM (Genes per kilobase million)") +
   scale_fill_brewer(palette="Paired") +
   # change order of legend
   # Add a theme
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   theme(axis.text.x=element_blank(),
         axis.ticks.x=element_blank()) 
 
 all_3

 ## Siderophore synthesis genes #####
 all_syn <- read_excel("Siderophore synthesis.xlsx", 
                          sheet = "All_stats")
 all_syn$Group <- as.factor(all_syn$Group)
 all_syn$Group <- factor(all_syn$Group, levels=c('LacAMP', 'LacADN'))

 # Box plot facetted by "HMM"
 #tiff("test.tiff", units="in", width=5, height=5, res=300)
 all_syn_1 <- ggboxplot(all_syn, x = "Group", y = "sumTPM",
                          color = "Group",
                          add = "jitter",
                          facet.by = "HMM", short.panel.labs = FALSE) +  
   scale_fill_brewer(palette="BuPu") +
   theme(legend.position="none")
 # Use only p.format as label. Remove method name.
 all_syn_1 + stat_compare_means(method = "t.test",label.x = 1.5)
 all_syn_1 + stat_compare_means( label =  "p.signif", label.x = 1.5) + 
   guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
   ylab("GPM (Genes per kilobase million)") 
 dev.off()
 
 all_syn_1 <- ggplot(all_syn, aes(x=Group, y=sumTPM, fill=Group)) + 
   geom_boxplot(outlier.shape = NA) +
   facet_wrap(~HMM, scale="free_x") +
   stat_compare_means(paired = TRUE, label =  "p.signif", label.x = 1.5) + 
   #guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
   ylab("GPM (Genes per kilobase million)") +
   scale_fill_brewer(palette="Paired") +
   # change order of legend
   # Add a theme
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   theme(axis.text.x=element_blank(),
         axis.ticks.x=element_blank() ) 
 
 all_syn_1
 
 ## Pvd + pvs with subunits ####
 all_syn_pvd_pvs <- read_excel("Siderophore synthesis.xlsx", 
                       sheet = "Pvd_Pvs_sub")
 all_syn_pvd_pvs$Group <- as.factor(all_syn_pvd_pvs$Group)
 all_syn_pvd_pvs$Group <- factor(all_syn_pvd_pvs$Group, levels=c('LacAMP', 'LacADN'))
 all_syn_pvd_pvs$sumTPM <- as.numeric(all_syn_pvd_pvs$sumTPM)
 
# all_syn_pvd_pvs_1 <- ggplot(all_syn_pvd_pvs, aes(x=Group, y=sumTPM, fill=Group)) + 
   #geom_boxplot(outlier.shape = NA) +
   #facet_wrap(~HMM, scale="free_x") +
   #stat_compare_means(paired = TRUE, label =  "p.signif", label.x = 1.5) + 
   #ylab("GPM (Genes per kilobase million)") +
   #scale_fill_brewer(palette="Paired") +
   #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   #theme(axis.text.x=element_blank(),
    #     axis.ticks.x=element_blank() ) 
 
# all_syn_pvd_pvs_1
 
 all_syn_pvd_pvs_1 <- ggplot(all_syn_pvd_pvs) +
   aes(x = Group, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
   scale_fill_brewer(palette="Paired") +
   ggtitle("Pvd & Pvs sub-families of siderophore synthesis genes") +
   labs(y = "GPM (Genes per kilobase million)") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 all_syn_pvd_pvs_1
 
 ## Siderophore transport genes ####
 all_trans <- read_excel("Siderophore transport.xlsx", 
                       sheet = "All_stats")
 all_trans$Group <- as.factor(all_trans$Group)
 all_trans$Group <- factor(all_trans$Group, levels=c('LacAMP', 'LacADN'))
 
 # Box plot facetted by "HMM"
 #tiff("test.tiff", units="in", width=5, height=5, res=300)
 all_trans_1 <- ggboxplot(all_trans, x = "Group", y = "sumTPM",
                        color = "Group",
                        add = "jitter",
                        facet.by = "HMM", short.panel.labs = FALSE) +  
   scale_fill_brewer(palette="BuPu") +
   theme(legend.position="none")
 all_trans_1
 # Use only p.format as label. Remove method name.
 all_trans_1 + stat_compare_means(method = "t.test",label.x = 1.5)
 all_trans_1 + stat_compare_means( label =  "p.signif", label.x = 1.5) + 
   guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
   ylab("GPM (Genes per kilobase million)") 
 dev.off()
 
 all_trans_1 <- ggplot(all_trans, aes(x=Group, y=sumTPM, fill=Group)) + 
   geom_boxplot(outlier.shape = NA) +
   facet_wrap(~HMM, scale="free_x") +
   stat_compare_means(paired = TRUE, label =  "p.signif", label.x = 1.5) + 
   #guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  
   ylab("GPM (Genes per kilobase million)") +
   scale_fill_brewer(palette="Paired") +
   # change order of legend
   # Add a theme
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   theme(axis.text.x=element_blank(),
         axis.ticks.x=element_blank() ) 
 
 all_trans_1
 
 ## Fpv + Pvu with subunits ####
 all_trans_fpv_pvu <- read_excel("Siderophore transport.xlsx", 
                               sheet = "Fpv_Pvu_sub")
 all_trans_fpv_pvu$Group <- as.factor(all_trans_fpv_pvu$Group)
 all_trans_fpv_pvu$Group <- factor(all_trans_fpv_pvu$Group, levels=c('LacAMP', 'LacADN'))
 all_trans_fpv_pvu$sumTPM <- as.numeric(all_trans_fpv_pvu$sumTPM)
 
 # all_syn_pvd_pvs_1 <- ggplot(all_syn_pvd_pvs, aes(x=Group, y=sumTPM, fill=Group)) + 
 #geom_boxplot(outlier.shape = NA) +
 #facet_wrap(~HMM, scale="free_x") +
 #stat_compare_means(paired = TRUE, label =  "p.signif", label.x = 1.5) + 
 #ylab("GPM (Genes per kilobase million)") +
 #scale_fill_brewer(palette="Paired") +
 #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
 #     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
 #theme(axis.text.x=element_blank(),
 #     axis.ticks.x=element_blank() ) 
 
 # all_syn_pvd_pvs_1
 
 all_trans_fpv_pvu_1 <- ggplot(all_trans_fpv_pvu) +
   aes(x = Group, y = sumTPM, fill = HMM) +
   geom_bar(position = "stack", stat = "identity", show.legend = TRUE) +
   scale_fill_brewer(palette="Paired") +
   ggtitle("Fpv & Pvu sub-families of siderophore transport genes") +
   labs(y = "GPM (Genes per kilobase million)") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.6)) +
   theme(axis.title.y = element_text(size = 13L, face = "bold")) +
   theme(axis.title.x = element_text(size = 13L,
                                     face = "bold"))
 all_trans_fpv_pvu_1
 
 ## combine plots into one (siderophore synthesis + transport)
  fig1 <- ggarrange(all_syn_pvd_pvs_1,all_trans_fpv_pvu_1,
                    ncol=2, nrow=1)
  fig1
  
