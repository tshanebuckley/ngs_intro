# trims the low-quality bases from the ends of the reads
# note: cuts only the bases whose qualities are below the given threshold
cd preprocessing
# run the trimming
fastq_quality_trimmer \
  -i filtered.fastq \
  -t 28 \
  -o trimmed.fastq \
  -Q33 # specify this as Phred+33
# generate the fastqc report on the trimmed data
fastqc --limits ../limits.txt trimmed.fastq
