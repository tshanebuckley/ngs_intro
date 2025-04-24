# in this script, we will convert all of our sam files to bam file, the binary version of the sam file
mkdir -p bam
# convert the BWA-MEM sam result to bam
samtools view \
  -uS \
  -o bam/SRR769545_mem.bam \
  sam/SRR769545_mem.sam

# convert the BWA-SW sam result to bam
samtools view \
  -uS \
  -o bam/SRR769545_bwasw.bam \
  sam/SRR769545_bwasw.sam

# convert the BWA-Backtrack sam results to bam 
samtools view \
  -uS \
  -o bam/SRR769545_aln.bam \
  sam/SRR769545_aln.sam

# convert the Bowtie sam resuls to bam
samtools view \
  -uS \
  -o bam/SRR769545_bowtie.bam \
  sam/SRR769545_bowtie2.sam

# convert the STAR sam results to bam
samtools view \
  -uS \
  -o bam/SRR769545_star.bam \
  sam/SRR769545_starAligned.out.sam

