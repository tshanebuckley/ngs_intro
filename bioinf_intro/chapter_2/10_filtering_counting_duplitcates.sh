# in this script, we will do some basic methods that we may perform on the bam file
cd processed

# 1. Extracting a chromosome
samtools view SRR769545_mem_sorted.bam NC_000011.10 > chr11_human.sam

# 2. Count the chimeric reads
touch chimeric_count.txt
samtools view SRR769545_mem_sorted.bam | grep 'SA:' | wc -l > chimeric_count.txt

# 3. Count the number of reads
touch read_count.txt
samtools view -c SRR769545_mem_sorted.bam > read_count.txt

# 4. Count the mapped reads
touch mapped_read_count.txt
samtools view -F 0x4 SRR769545_mem_sorted.bam > mapped_read_count.txt

# 5. Count the mapped reads
touch unmapped_read_count.txt
samtools view -f 0x4 SRR769545_mem_sorted.bam > unmapped_read_count.txt

# 6. Count the insertions and deletions from the mapped reads
# NOTE: you could count just insertion or deletions using 'I' or 'D', respectively
touch insertions_and_deletions.txt
samtools view -F 0x4 SRR769545_mem_sorted.bam \
  | cut -f 6 \
  | grep -P '[ID]' \
  | tr -cd '[ID]' \
  | wc -c \
  > insertions_and_deletions.txt

# 7. Remove duplicates from the read
samtools rmdup \
  SRR769545_mem_sorted.bam \
  SRR769545_mem_rmdup.bam 2> SRR769545_mem_rmdup.log
