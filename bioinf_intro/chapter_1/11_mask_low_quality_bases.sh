# masks low-quality bases inplace of a quality filter
# note: this is typically done as most aligners do NOT support processing paired-end FASTQC files with singletons (the fastq_quality_filter usually results in singletons/reads without pairs)
cd preprocessing
# run the masking of the original data
# mask bases with a Phred score <20
fastq_masker \
  -q 20 \   -i clipped.fastq \
  -o masked.fastq \
  -Q33
# generated the the fastqc report on the masked data
fastqc --limits ../limits.txt masked.fastq
