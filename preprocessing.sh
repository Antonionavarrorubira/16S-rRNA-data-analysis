#!/bin/bash

TMPDIR=/scratch/$USER/$SLURM_JOBID # temporal directory for the script
export TMPDIR # TMPDIR directory variable is exported for availability from other commands
mkdir -p $TMPDIR # temporal directory is created

# Activating qiime2 conda enviroment

conda activate qiime2

# Importing data into .qza format 

qiime tools import \
  --type EMPPairedEndSequences \
  --input-path paired-end-sequences \
  --output-path paired-end-sequences.qza

# generating .qzv file to visualize quality plot and infer filtering parameters for further steps

qiime demux summarize \
  --i-data paired-end-sequences.qza \
  --o-visualization paired-end-sequences.qzv

# denoising the sequences to generate ASVs
# trimming and truncation parameters are 

qiime dada2 denoise-paired \
--i-demultiplexed-seqs paired-end-raw.qza \
--p-trim-left-f 17 \ # 17nt forward illumina v3-v4 primer lenght
--p-trim-left-r 21 \ # 21nt reverse illumina v3-v4 primer lenght
--p-trunc-len-f 230 \
--p-trunc-len-r 200 \
--o-table ./results/tableNoFilt.qza \                      # ASV count table
--o-representative-sequences ./results/repseqsNoFilt.qza \ #the representative sequences for in each read
--o-denoising-stats ./results/denoising-statsNoFilt.qza \  #the status of the denoising process in a full catalogue
--p-n-threads 10                                            #number of the processors (depending on the machine the script is run)

# training silva-138.1 non-redundant database in order to generate a taxonomical classifier for the sequences.

qiime feature-classifier extract-reads \
    --i-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
    --p-f-primer CCTACGGGNGGCWGCAG \
    --p-r-primer GACTACHVGGGTATCTAATCC \
    --p-n-jobs 24 \
    --p-read-orientation 'forward' \
    --o-reads silva-138.1-ssu-nr99-seqs-515f-806r.qza

qiime rescript dereplicate \
    --i-sequences silva-138.1-ssu-nr99-seqs-515f-806r.qza \
    --i-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-515f-806r-uniq.qza \
    --o-dereplicated-taxa  silva-138.1-ssu-nr99-tax-515f-806r-derep-uniq.qza

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138.1-ssu-nr99-seqs-515f-806r-uniq.qza \
    --i-reference-taxonomy silva-138.1-ssu-nr99-tax-515f-806r-derep-uniq.qza \
    --o-classifier silva138_AB_V3-V4_classifier.qza

# using the trained classifier (silva138_AB_v3-V4_classifier.qza) to do the taxonomic classification of the repSeqs. 

qiime feature-classifier classify-sklearn \ #scikit learn package for ml and classification in python
--i-reads results/repseqsNoFilt.qza \ # non filtered representative ASV nucleotide sequences as the input
--i-classifier silva138_AB_V3-V4_classifier.qza \ # SILVA 138 trained classifier (v3-v4)
--o-classification taxonomy/taxonomyNoFilt_test.qza
--p-n-jobs 10 # taxonomic annotation of the repseqs

qiime fragment-insertion sepp 
--i-representative-sequences results/repseqsNoFilt.qza 
--i-reference-database sepp-refs-gg-13-8.qza 
--o-tree taxonomy/treeNoFilt.qza 
--o-placements taxonomy/tree-placementsNoFilt.qza 
--p-threads 10

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID

