
anvi-gen-contigs-database -f contigs.fa -o contigs.db -n "my metagenome"

# /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/nr_contigs-fixed.db

anvi-run-hmms -c nr_contigs-fixed.db -T 16 --just-do-it

anvi-run-ncbi-cogs -c nr_contigs-fixed.db  -T 16

anvi-get-sequences-for-gene-calls -c nr_contigs-fixed.db -o /shared/projects/bingo2_rmt/co_assembly/anvio/gene_calls.fa

###########
centrifuge -f -x /media/eclipse/centrifuge_db/nt/nt /shared/projects/bingo2_rmt/co_assembly/anvio/gene_calls.fa -S /shared/projects/bingo2_rmt/co_assembly/anvio/centrifuge_hits.tsv -p 20
anvi-import-taxonomy-for-genes -c contigs.db -i centrifuge_report.tsv centrifuge_hits.tsv -p centrifuge
##########

# for i in $(cat samples.txt) 
# do  
#   anvi-profile -i "$i".bam -c contigs.db -T 4 
# done

nano anvi-profile.sh

#!/bin/bash
#
#SBATCH -o anvi-profile.out
#SBATCH -e anvi-profile.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

conda activate anvio-8

anvi-profile -i ADN_64.bam -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db  -T 16 
anvi-profile -i ADN_65.bam -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db  -T 16
anvi-profile -i ADN_66.bam -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db  -T 16
anvi-profile -i AMP_67.bam -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db  -T 16
anvi-profile -i AMP_68.bam -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db  -T 16
anvi-profile -i AMP_69.bam -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db  -T 16

sbatch -A bingo2_rmt anvi-profile.sh

############# Cluster-contigs

nano anvi-profile_cluster.sh

#!/bin/bash
#
#SBATCH -o anvi-profile_cluster.out
#SBATCH -e anvi-profile_cluster.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

module load anvio/7.1

anvi-profile -i ADN_64.bam -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db  -T 16 -o ADN_64_cluster --cluster-contigs
anvi-profile -i ADN_65.bam -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db  -T 16 -o ADN_65_cluster --cluster-contigs
anvi-profile -i ADN_66.bam -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db  -T 16 -o ADN_66_cluster --cluster-contigs
anvi-profile -i AMP_67.bam -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db  -T 16 -o AMP_67_cluster --cluster-contigs
anvi-profile -i AMP_68.bam -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db  -T 16 -o AMP_68_cluster --cluster-contigs
anvi-profile -i AMP_69.bam -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db  -T 16 -o AMP_69_cluster --cluster-contigs

sbatch -A bingo2_rmt anvi-profile_cluster.sh

# anvi-merge */PROFILE.db -o merged_profile -c contigs.db

###############

 nano anvi-merge_all.sh

#!/bin/bash
#
#SBATCH -o anvi-merge_all.out
#SBATCH -e anvi-merge_all.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

conda activate anvio-7.1

anvi-merge ADN_64/PROFILE.db ADN_65/PROFILE.db ADN_66/PROFILE.db AMP_67/PROFILE.db AMP_68/PROFILE.db AMP_69/PROFILE.db -o merged_profile_all -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/ nr_contigs-fixed.db --enforce-hierarchical-clustering 

sbatch -A bingo2_rmt anvi-merge_all.sh


###################

nano anvi-merge_all_cluster.sh

#!/bin/bash
#
#SBATCH -o anvi-merge_all_cluster.out
#SBATCH -e anvi-merge_all_cluster.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

module load anvio/7.1


anvi-merge ADN_64/PROFILE.db ADN_65/PROFILE.db ADN_66/PROFILE.db AMP_67/PROFILE.db AMP_68/PROFILE.db AMP_69PROFILE.db -o merged_profile_all -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/nr_contigs-fixed.db --enforce-hierarchical-clustering 

sbatch -A bingo2_rmt anvi-merge_all_cluster.sh

###############################. ADN
nano anvi-merge_ADN.sh

#!/bin/bash
#
#SBATCH -o anvi-merge_ADN.out
#SBATCH -e anvi-merge_ADN.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

conda activate anvio-7.1

anvi-merge ADN_64/PROFILE.db ADN_65/PROFILE.db ADN_66/PROFILE.db  -o merged_profile_ADN -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/nr_contigs-fixed.db --enforce-hierarchical-clustering 
 
sbatch -A bingo2_rmt anvi-merge_ADN.sh



###################     AMP
nano anvi-merge_AMP.sh

#!/bin/bash
#
#SBATCH -o anvi-merge_AMP.out
#SBATCH -e anvi-merge_AMP.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

conda activate anvio-7.1

anvi-merge AMP_67/PROFILE.db AMP_68/PROFILE.db AMP_69/PROFILE.db -o merged_profile_AMP -c /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/nr_contigs-fixed.db --enforce-hierarchical-clustering 

sbatch -A bingo2_rmt anvi-merge_AMP.sh


anvi-run-scg-taxonomy -T 10 -c nr_contigs-fixed.db


####VISUALISATION ####

 
anvi-interactive  -p  merged_profile_ADN/PROFILE.db  -c   nr_contigs-fixed.db

anvi-interactive  -p  merged_profile_AMP/PROFILE.db  -c   nr_contigs-fixed.db

anvi-interactive -p merged_profile_all/PROFILE.db -c    nr_contigs-fixed.db
 
anvi-interactive -p merged_profile_all_cluster/PROFILE.db   -c nr_contigs-fixed.db
 
 
anvi-show-collections-and-bins p PROFILE.db --server-only -P 8080

############## 
anvi-estimate-scg-taxonomy -c nr_contigs-fixed.db --metagenome-mode outl0strokewidth0 --output-file  taxa_MAGs_single.txt

anvi-estimate-scg-taxonomy -c nr_contigs-fixed.db -p merged_profile_all/PROFILE.db --metagenome-mode --compute-scg-coverages --output-file taxonomy_MAGs.txt

anvi-estimate-scg-taxonomy -c nr_contigs-fixed.db -p merged_profile_all/PROFILE.db --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy

anvi-interactive merged_profile_all/PROFILE.db -c nr_contigs-fixed.db --list-collections

anvi-show-collections-and-bins -p merged_profile_all/PROFILE.db

anvi-refine -p merged_profile_all/PROFILE.db -c nr_contigs-fixed.db Vibrio_MAG8


########################################## AUTOMATIC BINNING - LOCAL

anvi-show-collections-and-bins -p merged_profile_all/PROFILE.db
anvi-summarize -p merged_profile_all/PROFILE.db -c nr_contigs-fixed.db -C All -o Bins 

anvi-import-collection -c nr_contigs-fixed.db -p merged_profile_all/PROFILE.db --contigs-mode -C ALL merged_all

anvi-estimate-genome-completeness -p merged_profile_all/PROFILE.db -c nr_contigs-fixed.db -C All

anvi-estimate-scg-taxonomy -p merged_profile_all/PROFILE.db -c nr_contigs-fixed.db -C All --output-file taxa_MAGs_all.txt

################# Phylogenomics - BACTERIA 71

anvi-get-sequences-for-hmm-hits -c nr_contigs-fixed.db -p merged_profile_all/PROFILE.db -o phylogenomics.fa --list-hmm-sources

anvi-get-sequences-for-hmm-hits -c nr_contigs-fixed.db
fs24 outl0strokewidth0 strokec12  strokec16 -p merged_profile_all/PROFILE.db
-o seqs-for-phylogenomics.fa --hmm-sourcestrokec12  Bacteria_71  --list-available-gene-names


anvi-get-sequences-for-hmm-hits  -c nr_contigs-fixed.db -p merged_profile_all/PROFILE.db \
    -o seqs-for-phylogenomics.fa --hmm-source  Bacteria_71  -C All --gene-names  Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6

anvi-get-sequences-for-hmm-hits  -c nr_contigs-fixed.db -p merged_profile_all/PROFILE.db \
    -o  seqs-for-phylogenomics.fa --hmm-source  Bacteria_71  -C All  --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
     --concatenate-genes   --return-best-hit

anvi-get-sequences-for-hmm-hits  -c nr_contigs-fixed.db -p merged_profile_all/PROFILE.db \
    -o  seqs-for-phylogenomics.fa --hmm-source Bacteria_71 -C  All --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
    --concatenate-genes --return-best-hit --get-aa-sequences

anvi-gen-phylogenomic-tree -f seqs-for-phylogenomics.fa -o  phylogenomic-tree.txt

anvi-interactive --tree phylogenomic-tree.txt -p temp-profile.db --title "Phylogenomics of BINGO2 Bins" --manual

anvi-interactive -p merged_profile_all/PROFILE.db -c nr_contigs-fixed.db -C All --tree phylogenomic-tree.txt

########################### Manual binning 
Metabat

 nano metabat_bins.sh

#!/bin/bash
#
#SBATCH -o anvi-merge_all.out
#SBATCH -e anvi-merge_all.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

module load metabat2/2.15

jgi_summarize_bam_contig_depths --outputDepth Mags/depth.txt *.bam

sbatch -A bingo2_rmt metabat_bins.sh

metabat2 -i /shared/projects/bingo2_rmt/co_assembly/megahit/contigs-fixed.fa \
    -a /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/depth.txt \
    -t 4 -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/Binning/Bin

grep -c ">" Binning/*.fa | awk -F: '{ s+=$2 } END { print s }'

grep -v ">" Binning/*.fa |wc -m

################### Maxbin

jgi_summarize_bam_contig_depths  --outputDepth  /Mags/depth_maxbin.txt  *.bam

# Use same depth.txt 
# Format the depth file into MaxBin 
# I don't recommend using only "totalAvgDepth" for MaxBin as written in the Pelikan et al., 2020
# I prefer treat each sample separately

# MaxBin Usage: Run maxbin

nano maxbin_bins.sh

#!/bin/bash
#
#SBATCH -o errors/maxbin_all.out
#SBATCH -e errors/maxbin_all.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

module load maxbin2/2.2.7

run_MaxBin.pl  -contigs /shared/projects/bingo2_rmt/co_assembly/megahit/contigs-fixed.fa \
     -out /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/maxbin \
     -abund_list /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/maxbin/abund_list  -markerset  -prob_threshold

sbatch -A bingo2_rmt maxbin_bins.sh

# Rename the fasta files in $WORKDIR, and move them to $WORKDIR/fasta_bins 
$WORKDIR /maxbin.*.fasta
strokec33 dostrokec27 
     $NEW = $( echo  $OLD   |  sed 's/.*/maxbin./maxbin.bin_/g'  |  sed 's/.fasta$/.fa/g')
    mv $OLD $WORKDIR/fasta_bins/$NEW


 ###### Concoct
cut_up_fasta.py contigs-fixed.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

nano concoct_bins.sh

#!/bin/bash
#
#SBATCH -o errors/concoct_all.out
#SBATCH -e errors/concoct_all.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

module load concoct/1.1.0

concoct_coverage_table.py /shared/projects/bingo2_rmt/co_assembly/megahit/contigs_10K.bed /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/*.bam > /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/coverage_table.tsv

concoct --composition_file /shared/projects/bingo2_rmt/co_assembly/megahit/contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/

merge_cutup_clustering.py /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/concoct/clustering_gt1000.csv > /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/concoct/clustering_merged.csv

sbatch -A bingo2_rmt concoct_bins.sh

mkdir concoct_output/fasta_bins
extract_fasta_bins.py original_contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins

###### Vamb 

vamb --outdir /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/vamb fasta 
 /shared/projects/bingo2_rmt/co_assembly/megahit --bamfiles /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/*.bam -o C


########### DASTOOL
module load das_tool/1.1.3


f2 cf33 cb28 outl0strokewidth0  for7 cb28   cf32 strokec32 OLD7   cf33 cb28 strokec33 in7 cb28 strokec27  cf32 strokec32 /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags7 strokec27 /maxbin.*.fastacb1 
cf33 cb28 strokec33 do7 cb1 strokec27 
pardpardeftab720partightenfactor0
7 cb28     cf32 strokec32 NEW7 strokec27 =9 strokec29 $(cf35 cb28 strokec35 echo9 cb28 strokec29  cf32 strokec32 $OLD9 strokec29  7 strokec27 |9 strokec29  sed '92s/.*/Mags/maxbin.bin_/g'92 7 strokec27 |9 strokec29  sed 's/.fasta$/.fa/g')7 cb1 strokec27 
cb28     mv cf32 strokec32 $OLD7 strokec27  cf32 strokec32 /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/maxbin9 strokec29 /cf32 strokec32 $NEW7 cb1 strokec27 
pardpardeftab720partightenfactor0
cf33 cb28 strokec33 done7 cb1 strokec27 
f0 cf4  outl0strokewidth0 

pardpardeftab720partightenfactor0

f2 cf33 cb28 outl0strokewidth0 strokec33 for7 cb28 strokec27  cf32 strokec32 BIN7 strokec27  cf33 cb28 strokec33 in7 cb28 strokec27  cf32 strokec32 /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags7 strokec27 /*.fastacb1 
cf33 cb28 strokec33 do7 cb1 strokec27 
pardpardeftab720partightenfactor0
7 cb28     cf32 strokec32 BIN_ID7 strokec27 =9 strokec29 $(cf35 cb28 strokec35 echo9 cb28 strokec29  cf32 strokec32 $BIN9 strokec29  7 strokec27 |9 strokec29  sed '92s/.*Mags///g'92 7 strokec27 |9 strokec29  sed 's/.fasta$//g')7 cb1 strokec27 
cb28     cf33 cb28 strokec33 while7 cb28 strokec27  cf35 cb28 strokec35 read7 cb28 strokec27  9 strokec29 LINE7 cb1 strokec27 
cb28     cf33 cb28 strokec33 do7 cb1 strokec27 
cb28         cf33 cb28 strokec33 if7 cb28 strokec27  [[ cf32 strokec32 $LINE7 strokec27  =~ ^9 strokec29 ">"7 strokec27  ]]; cf33 cb28 strokec33 then7 cb1 strokec27 
cb28             cf32 strokec32 CONTIG_ID7 strokec27 =9 strokec29 $(cf35 cb28 strokec35 echo9 cb28 strokec29  cf32 strokec32 $LINE9 strokec29  7 strokec27 |9 strokec29  sed 's/^>//g' 7 strokec27 |9 strokec29  sed 's/ .*//g')7 cb1 strokec27 
cb28             cf35 cb28 strokec35 echo7 cb28 strokec27  9 strokec29 -e7 strokec27  9 strokec29 "cf32 strokec32 $CONTIG_ID9 strokec29 tcf32 strokec32 $BIN_ID9 strokec29 "7 strokec27  >> cf36 7 outl0strokewidth0 /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/dastool9 cb28 outl0strokewidth0 strokec29 /maxbin_contigs2bin.tsv7 cb1 strokec27 
cb28         cf33 cb28 strokec33 fi7 cb1 strokec27 
cb28     cf33 cb28 strokec33 done7 cb28 strokec27  < cf32 strokec32 $BIN7 cb1 strokec27 
pardpardeftab720partightenfactor0
cf33 cb28 strokec33 done7 cb1 strokec27 
f0 cf4  outl0strokewidth0 
pardpardeftab720partightenfactor0

f2 cf10 cb11 outl0strokewidth0 strokec33 forstrokec27  strokec32 BINstrokec27  strokec33 instrokec27  strokec32 /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/
f5fs27fsmilli13600 cf10 outl0strokewidth0 Binning/*.fa
f2fs24 cf10 outl0strokewidth0 strokec27 
strokec33 dostrokec27 
    #strokec32 BIN_IDstrokec27 =strokec29 $(strokec35 echostrokec29  strokec32 $BINstrokec29  strokec27 |strokec29  sed '93s/.*fasta_bins///g'94 strokec27 |strokec29  sed 's/.fa$//g')
    BIN_ID=${BIN%.fa}strokec27 
    strokec33 whilestrokec27  strokec35 readstrokec27  strokec29 LINEstrokec27 
    strokec33 dostrokec27 
        strokec33 ifstrokec27  [[ strokec32 $LINEstrokec27  =~ ^strokec29 ">"strokec27  ]]; strokec33 thenstrokec27 
            strokec32 CONTIG_IDstrokec27 =strokec29 $(strokec35 echostrokec29  strokec32 $LINEstrokec29  strokec27 |strokec29  sed 's/^>//g' strokec27 |strokec29  sed 's/ .*//g')strokec27 
            strokec35 echostrokec27  strokec29 -estrokec27  strokec29 "strokec32 $CONTIG_IDstrokec29 tstrokec32 $BIN_IDstrokec29 "strokec27  >> strokec32 /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/dastoolstrokec29 /metabat_contigs2bin.tsvstrokec27 
        
strokec33 		fistrokec27 
    strokec33 donestrokec27  < strokec32 $BINstrokec27 
strokec33 done

# Download and extract DASTool.zip archive:
unzip 1.1.1.zip
cd ./DAS_Tool-1.1.1

# Unzip SCG database:
unzip ./db.zip -d db

# Run DAS Tool:
./DAS_Tool -h

DAS_Tool -i /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/dastool/maxbin_contigs2bin.tsv /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/dastool/metabat_contigs2bin.tsv \
    -l maxbin,metabat -c /shared/projects/bingo2_rmt/co_assembly/megahit/contigs-fixed.fa \
    -o /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/dastool/das_tool --write_bins \
    --score_threshold -t

export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/Users/username/Programs"

anvi-interactive -p merged_profile/PROFILE.db -c nr_contigs-fixed.db --server-only -P 8080

/shared/ifbstor1/software/miniconda/envs/das_tool-1.1.3/bin/DAS_Tool
######  CheckM

pip3 install numpy
pip3 install matplotlib
pip3 install pysam
pip3 install checkm-genome

 conda create -n checkm python=3.9
conda activate checkm
conda install -c bioconda numpy matplotlib pysam
conda install -c bioconda hmmer prodigal pplacer
pip3 install checkm-genome

conda create -n checkm checkm-genome=1.2.1

checkm data setRoot /shared/projects/bingo2_rmt/miniconda3/mambaforge/envs/checkm
export CHECKM_DATA_PATH=/shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/checkm

nano checkm_maxbin_bins.sh

#!/bin/bash
#
#SBATCH -o checkm_maxbin_all.out
#SBATCH -e checkm_maxbin_all.err
#SBATCH --mail-type END
#SBATCH --mail-user rhea.thoppil@obs-banyuls.fr
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 30
#SBATCH --mem 100GB

conda activate checkm

checkm lineage_wf -t 12 -x fasta /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/ /shared/projects/bingo2_rmt/co_assembly/04_MAPPING/bam/Mags/checkm/maxbin


sbatch -A bingo2_rmt checkm_maxbin_bins.sh

