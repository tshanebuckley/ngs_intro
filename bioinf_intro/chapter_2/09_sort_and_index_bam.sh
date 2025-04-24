# in this script, we will use samtools to perform some extra processing steps on produced bam file.
# specifically, we will sort the alignments, and then index them
mkdir -p processed

# sort the alignments by coordinate order
samtools sort \
  -@ 4 \
  -T processed/mem.tmp.sort \
  -o processed/SRR769545_mem_sorted.bam \
  bam/SRR769545_mem.bam

# index the sorted alignments
cd processed
samtools index -@ 4 SRR769545_mem_sorted.bam
