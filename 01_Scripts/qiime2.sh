#!/bin/bash

# 1. Import Data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest202406.csv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33

# 2. Cutadapt (Remove Primers)
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux.qza \
  --p-front-f NNNNTGACTNNNNAACMGGATTAGATACCCKG \
  --p-front-r NNTGACTNNNCGTCATCCMCACCTTCCTC \
  --p-error-rate 0.1 \
  --p-discard-untrimmed \
  --o-trimmed-sequences trimmed-seqs.qza \
  --verbose

# 3. DADA2 (Denoise)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-seqs.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --p-n-threads 0

# 4. Decontam (Identify & Remove Contaminants)
qiime quality-control decontam-identify \
  --i-table table.qza \
  --m-metadata-file metadata.tsv \
  --p-method prevalence \
  --p-prev-control-column sample-type \
  --p-prev-control-indicator negative \
  --o-decontam-scores decontam-scores.qza

qiime quality-control decontam-remove \
  --i-decontam-scores decontam-scores.qza \
  --i-table table.qza \
  --p-threshold 0.10 \
  --o-filtered-table table-nocontam.qza

qiime feature-table filter-seqs \
  --i-data rep-seqs.qza \
  --i-table table-nocontam.qza \
  --o-filtered-data rep-seqs-nocontam.qza

# 5. Taxonomy Classification (Your Requirement: Annotate ALL sequences)
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138.1-ssu-nr99-V5-V7-classifier.qza \
  --i-reads rep-seqs-nocontam.qza \
  --o-classification taxonomy-nocontam.qza

# 6. Export RAW Tables (Your Requirement: ASV, Family, Genus)

# 6.1 Export Raw ASV Table
qiime tools export \
  --input-path table-nocontam.qza \
  --output-path exported-feature-table

biom convert \
  -i exported-feature-table/feature-table.biom \
  -o feature-table-asv-raw.tsv \
  --to-tsv

# 6.2 Collapse to Family and Export
qiime taxa collapse \
  --i-table table-nocontam.qza \
  --i-taxonomy taxonomy-nocontam.qza \
  --p-level 5 \
  --o-collapsed-table table-family-raw.qza

qiime tools export \
  --input-path table-family-raw.qza \
  --output-path exported-family-raw

biom convert \
  -i exported-family-raw/feature-table.biom \
  -o feature-table-family-raw.tsv \
  --to-tsv

# 6.3 Collapse to Genus and Export
qiime taxa collapse \
  --i-table table-nocontam.qza \
  --i-taxonomy taxonomy-nocontam.qza \
  --p-level 6 \
  --o-collapsed-table table-genus-raw.qza

qiime tools export \
  --input-path table-genus-raw.qza \
  --output-path exported-genus-raw

biom convert \
  -i exported-genus-raw/feature-table.biom \
  -o feature-table-genus-raw.tsv \
  --to-tsv