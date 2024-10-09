
ssh rthoppil@core.cluster.france-bioinformatique.fr

#######################################################
Co-assembly 

ls Adapter_trimmed/*fastq.gz

### environmental variables

R1s=`ls Adapter_trimmed/*R1* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
R2s=`ls Adapter_trimmed/*R2* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`

echo $R1s
echo $R2s


######## 
Based on parameter <1000

nano slurm_co-assembly_megahit.sh

#!/bin/bash
#
#SBATCH -o Adapter_trimmed/AFQFco_megahit.out
#SBATCH -e Adapter_trimmed/AFQFco_megahit.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition long
#SBATCH --cpus-per-task 40
#SBATCH --mem 200GB

R1s=`ls Adapter_trimmed/*R1* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
R2s=`ls Adapter_trimmed/*R2* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
megahit -1 $R1s -2 $R2s --min-contig-len 1000 -o megahit/Co-assembly/ -t 40

sbatch -A bingo2_rmt slurm_co-assembly_megahit.sh

###### re-format fasta.   - minimum length at 2000

anvi-script-reformat-fasta final.contigs.fa -o contigs-fixed.fa --min-len 2000 --simplify-names --report name_conversions.txt


###################################### FeGenie 

nano co_fegenie_iron.sh

#!/bin/bash
#
#SBATCH -o iron_genie/iron_genie.out
#SBATCH -e iron_genie/iron_genie.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition long
#SBATCH --cpus-per-task 40
#SBATCH --mem 200GB

module load fegenie
FeGenie.py -bin_dir megahit/ -bin_ext fa -out iron_genie/co-assembly_iron/ --meta --norm

sbatch -A bingo2_rmt co_fegenie_iron.sh


######################
nano co_fegenie_iron_gene_calls.sh

#!/bin/bash
#
#SBATCH -o iron_genie/iron_genie_gene_calls.out
#SBATCH -e iron_genie/iron_genie_gene_calls.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition long
#SBATCH --cpus-per-task 40
#SBATCH --mem 200GB

conda activate fegenie
FeGenie.py -bin_dir prodigal/ -bin_ext fna -out iron_genie/gene_calls/ --meta --norm -t 16

sbatch -A bingo2_rmt co_fegenie_iron_gene_calls.sh

########################################################

################ MAPPING

nano bowtie_build.sh

#!/bin/bash
#
#SBATCH -o 04_MAPPING/bowtie_build.out
#SBATCH -e 04_MAPPING/bowtie_build.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

conda activate anvio-7.1
bowtie2-build Co-assembly/contigs-fixed.fa 04_MAPPING/contigs

sbatch -A bingo2_rmt bowtie_build.sh

##################################
nano anvio_bowtie2.sh

#!/bin/bash
#
#SBATCH -o errors/bowtie2.out
#SBATCH -e errors/bowtie2.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition long
#SBATCH --cpus-per-task 40
#SBATCH --mem 200GB

conda activate anvio-7.1

bowtie2 -x /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/contigs -1 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-64_R1.fastq.gz -2 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-64_R2.fastq.gz -S /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_64.sam
bowtie2 -x /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/contigs -1 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-65_R1.fastq.gz -2 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-65_R2.fastq.gz -S /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_65.sam
bowtie2 -x /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/contigs -1 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-66_R1.fastq.gz -2 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-66_R2.fastq.gz -S /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_66.sam
bowtie2 -x /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/contigs -1 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-67_R1.fastq.gz -2 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-67_R2.fastq.gz -S /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_67.sam
bowtie2 -x /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/contigs -1 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-68_R1.fastq.gz -2 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-68_R2.fastq.gz -S /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_68.sam
bowtie2 -x /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/contigs -1 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-69_R1.fastq.gz -2 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-69_R2.fastq.gz -S /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_69.sam

sbatch -A bingo2_rmt anvio_bowtie2.sh

rthoppil14
#########################################
##########
nano samtools.sh

#!/bin/bash
#
#SBATCH -o errors/samtools.out
#SBATCH -e errors/samtools.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 40
#SBATCH --mem 200GB

module load samtools/1.15.1 

samtools view -F 4 -bS /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_64.sam> /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/raw_bam/ADN_64-RAW.bam
samtools view -F 4 -bS /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_65.sam> /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/raw_bam/ADN_65-RAW.bam
samtools view -F 4 -bS /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_66.sam> /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/raw_bam/ADN_66-RAW.bam
samtools view -F 4 -bS /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_67.sam > /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/raw_bam/AMP_67-RAW.bam
samtools view -F 4 -bS /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_68.sam > /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/raw_bam/AMP_68-RAW.bam
samtools view -F 4 -bS /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_69.sam > /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/raw_bam/AMP_69-RAW.bam

sbatch -A bingo2_rmt samtools.sh



##################### bedtools

nano bamtobed_all.sh

#!/bin/bash
#
#SBATCH -o errors/bamtobed_all.out
#SBATCH -e errors/bamtobed_all.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 16
#SBATCH --mem 100GB

module load bedtools/2.30.0

bedtools bamtobed -i /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/ADN_64.bam 
bedtools bamtobed -i /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/ADN_65.bam
bedtools bamtobed -i /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/ADN_66.bam
bedtools bamtobed -i /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/AMP_67.bam
bedtools bamtobed -i /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/AMP_68.bam
bedtools bamtobed -i /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/AMP_69.bam



sbatch -A bingo2_rmt bamtobed_all.sh


#########################################
nano anvi_bam.sh

#!/bin/bash
#
#SBATCH -o errors/anvi_bam.out
#SBATCH -e errors/anvi_bam.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 40
#SBATCH --mem 200GB

module load anvio
anvi-init-bam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/raw_bam/ADN_64-RAW.bam -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_64.bam
anvi-init-bam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/raw_bam/ADN_65-RAW.bam -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_65.bam
anvi-init-bam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/raw_bam/ADN_66-RAW.bam -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_66.bam
anvi-init-bam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/raw_bam/AMP_67-RAW.bam -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_67.bam
anvi-init-bam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/raw_bam/AMP_68-RAW.bam -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_68.bam
anvi-init-bam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/raw_bam/AMP_69-RAW.bam -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_69.bam

sbatch -A bingo2_rmt anvi_bam.sh

####################################

rm /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_64.sam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_64.sam-RAW.bam
rm /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_65.sam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_65.sam-RAW.bam
rm /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_66.sam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ADN_66.sam-RAW.bam
rm /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_67.sam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_67.sam-RAW.bam
rm /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_68.sam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_68.sam-RAW.bam
rm /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_69.sam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/AMP_69.sam-RAW.bam

###################################
nano anvio_db.sh

#!/bin/bash
#
#SBATCH -o errors/anvio_db.out
#SBATCH -e errors/anvio_db.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 16
#SBATCH --mem 100GB

conda activate anvio-7.1
anvi-gen-contigs-database -f '/shared/projects/bingo2_rmt/co_assembly/megahit/contigs-fixed.fa' -o '/shared/projects/bingo2_rmt/co_assembly/04_MAPPING/nr_contigs-fixed.db' 

sbatch -A bingo2_rmt anvio_db.sh

################ get contig- or gene-level coverage and detection stats 
nano anvio_blitz.sh

#!/bin/bash
quant_64
#SBATCH -o errors/anvio_blitz.out
#SBATCH -e errors/anvio_blitz.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 40
#SBATCH --mem 200GB

conda activate anvio-7.1

anvi-profile-blitz /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/ADN_64.bam  -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/nr_contigs-fixed.db  --gene-mode -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/coverage/ADN_64profile.txt
anvi-profile-blitz /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/ADN_65.bam  -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/nr_contigs-fixed.db  --gene-mode -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/coverage/ADN_65profile.txt
anvi-profile-blitz /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/ADN_66.bam  -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/nr_contigs-fixed.db  --gene-mode -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/coverage/ADN_66profile.txt
anvi-profile-blitz /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/AMP_67.bam  -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/nr_contigs-fixed.db  --gene-mode -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/coverage/AMP_67profile.txt
anvi-profile-blitz /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/AMP_68.bam  -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/nr_contigs-fixed.db  --gene-mode -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/coverage/AMP_68profile.txt
anvi-profile-blitz /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/AMP_69.bam  -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/nr_contigs-fixed.db  --gene-mode -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/coverage/AMP_69profile.txt

sbatch -A bingo2_rmt anvio_blitz.sh

############ without using ANVIO

nano bedtools.sh

#!/bin/bash
#
#SBATCH -o errors/bedtools.out
#SBATCH -e errors/bedtools.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

module load bedtools 

bedtools genomecov -ibam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/ADN_64.bam > /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bedtools/64_histogram.tab
bedtools genomecov -ibam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/ADN_65.bam > /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bedtools/65_histogram.tab
bedtools genomecov -ibam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/ADN_66.bam > /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bedtools/66_histogram.tab
bedtools genomecov -ibam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/AMP_67.bam > /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bedtools/67_histogram.tab
bedtools genomecov -ibam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/AMP_68.bam > /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bedtools/68_histogram.tab
bedtools genomecov -ibam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/AMP_69.bam > /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bedtools/69_histogram.tab

sbatch -A bingo2_rmt bedtools.sh

#############################
nano salmon_align64.sh

#!/bin/bash
#
#SBATCH -o align_salmon.out
#SBATCH -e align_salmon.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

salmon quant -t /shared/projects/bingo2_rmt/cd_hit/nr_whole_gene_anvio.fna -l A -a /shared/projects/bingo2_rmt/anvio/anvio_mapping/bedtools/Sample_64.bam -o /shared/projects/bingo2_rmt/salmon_alignment/salmon_quant 

sbatch -A bingo2_rmt salmon_align64.sh



#######################################
nano prodigal_normal_co.sh

#!/bin/bash
#
#SBATCH -o errors/prodigal_normal_co.out
#SBATCH -e errors/prodigal_normal_co.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

prodigal -i /shared/projects/bingo2_rmt/co_assembly/megahit/contigs-fixed.fa -o /shared/projects/bingo2_rmt/co_assembly/prodigal/AFQF_co.genes -d /shared/projects/bingo2_rmt/co_assembly/prodigal/AFQF_co.fna -a /shared/projects/bingo2_rmt/co_assembly/prodigal/AFQF_co.faa -p meta 

sbatch -A bingo2_rmt prodigal_normal_co.sh



###################
Diamond for blast 

diamond makedb --in /shared/projects/bingo2_rmt/co_assembly/prodigal/AFQF_co.faa -d /shared/projects/bingo2_rmt/co_assembly/diamond_blast -p 16

nano diamond_fe.sh

#!/bin/bash
#
#SBATCH -o errors/diamond_fe.out
#SBATCH -e errors/diamond_fe.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 200GB

module load diamond

diamond blastp --more-sensitive --threads 30 -e 1e-5 -d /shared/projects/bingo2_rmt/co_assembly/diamond_blast -q /shared/projects/bingo2_rmt/co_assembly/iron_genie/co-assembly_iron/ORF_calls/contigs-fixed.fa-proteins.faa -o /shared/projects/bingo2_rmt/co_assembly/diamond_blast/co_assembly_fe --max-target-seqs 1 -c1 --outfmt 6 qseqid sseqid qlen slen pident length gaps mismatch qstart qend sstart send evalue bitscore stitle 

sbatch -A bingo2_rmt diamond_fe.sh


######################
Salmon

nano salmon_index.sh

#!/bin/bash
#
#SBATCH -o errors/salmon_index.out
#SBATCH -e errors/salmon_index.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 200GB

module load salmon
salmon index -t /shared/projects/bingo2_rmt/co_assembly/prodigal/AFQF_co.fna -i /shared/projects/bingo2_rmt/co_assembly/salmon/indexes/transcripts --gencode

sbatch -A bingo2_rmt salmon_index.sh

##########################
nano salmon_quant_meta.sh

#!/bin/bash
#
#SBATCH -o errors/salmon_quant_meta.out
#SBATCH -e errors/salmon_quant_meta.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 16
#SBATCH --mem 100GB

module load salmon

salmon quant i /shared/projects/bingo2_rmt/co_assembly/salmon/indexes/transcripts -l A -1 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-64_R1.fastq.gz -2 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-64_R2.fastq.gz -o /shared/projects/bingo2_rmt/co_assembly/salmon/meta/transcripts_quant_64 —meta

salmon quant -i /shared/projects/bingo2_rmt/co_assembly/salmon/indexes/transcripts -l A -1 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-65_R1.fastq.gz -2 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-65_R2.fastq.gz -o /shared/projects/bingo2_rmt/co_assembly/salmon/meta/transcripts_quant_65 —meta

salmon quant -i /shared/projects/bingo2_rmt/co_assembly/salmon/indexes/transcripts -l A -1 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-66_R1.fastq.gz -2 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-66_R2.fastq.gz -o /shared/projects/bingo2_rmt/co_assembly/salmon/meta/transcripts_quant_66 —meta

salmon quant  -i /shared/projects/bingo2_rmt/co_assembly/salmon/indexes/transcripts -l A -1 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-67_R1.fastq.gz -2 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-67_R2.fastq.gz -o /shared/projects/bingo2_rmt/co_assembly/salmon/meta/transcripts_quant_67 —meta

salmon quant -i /shared/projects/bingo2_rmt/co_assembly/salmon/indexes/transcripts -l A -1 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-68_R1.fastq.gz -2 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-68_R2.fastq.gz -o /shared/projects/bingo2_rmt/co_assembly/salmon/meta/transcripts_quant_68 —meta

salmon quant -i /shared/projects/bingo2_rmt/co_assembly/salmon/indexes/transcripts -l A -1 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-69_R1.fastq.gz -2 /shared/projects/bingo2_rmt/Adapter_trimmed/AFQF-69_R2.fastq.gz -o /shared/projects/bingo2_rmt/co_assembly/salmon/metal/transcripts_quant_69 —meta

sbatch -A bingo2_rmt salmon_quant_meta.sh


#####################################
Feature counts

nano feature_counts.sh

#!/bin/bash
#
#SBATCH -o errors/feature_counts.out
#SBATCH -e errors/feature_counts.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 16
#SBATCH --mem 200GB

conda activate featurecounts 

featureCounts -a /shared/projects/bingo2_rmt/co_assembly/prodigal/gff_mode/AFQF_co.gtf -o /shared/projects/bingo2_rmt/co_assembly/feature_counts/MetaG_Counts.txt -F GTF -t transcript -g gene_id -O --fracOverlap 0.25 -Q 1 --primary --ignoreDup -s 0 -p -B -P -C -T 48 /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/ADN_64.bam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/ADN_65.bam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/ADN_66.bam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/AMP_67.bam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/AMP_68.bam /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/AMP_69.bam 

sbatch -A bingo2_rmt  feature_counts.sh




############################


