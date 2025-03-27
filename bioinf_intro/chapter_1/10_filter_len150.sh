# filter the clipped reads whose len <150 so that all reads are a length of 150
# note: we know that 150 was the max len from the clipped.fastq report
cd preprocessing
# filter out the shorter than 150bp reads
paste <(cat clipped.fastq | paste - - - -) \
  | awk -v FS='\t' 'length($2) >= 150 && length($4) >= 150' \
  | tee >(cut -f 1-4 | tr '\t' '\n'>len150.fastq)
# generate the fastqc report on the filtered reads
fastqc --limits ../limits.txt len150.fastq
