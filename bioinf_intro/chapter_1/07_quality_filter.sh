# filters the base reads that have:
# - 80% of the bases
# - which have a quality score of 28 or greater
cd preprocessing
# run the filtering
fastq_quality_filter \
  -i bad.fastq \
  -q 28 \ # specify the minimum Phred quality threshold
  -p 80 \ # % of bases with at least the 28 quality threshold
  -o filtered.fastq \
  -Q33 # specify this as Phred+33, default is Phred+64
# generate a fastqc report on the filtered data
fastqc --limits ../limits.txt filtered.fastq
