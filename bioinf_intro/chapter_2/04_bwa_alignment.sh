# in this script, we will be indexing using the 3 options from bwa
cd refgenome
# index using bwa
bwa index GRCh38.p13_ref.fna
# return to the parent directory
cd ..

# bwa has 3 potential algorithms that can be used for alignment

# 1. BWA-MEM algorithm
bwa mem \
  -t 4 \
  refgenome/GRCh38.p13_ref.fna \
  data/SRR769545_1.fastq.gz \
  data/SRR769545_2.fastq.gz \
  > sam/SRR769545_mem.sam 2> sam/SRR769545_mem.log

# 2. BWA-SW algorithm
bwa bwasw \
  -t 4 \
  refgenome/GRCh38.p13_ref.fna \
  data/SRR769545_1.fastq.gz \
  data/SRR769545_2.fastq.gz \
  > sam/SRR769545_bwasw.sam 2> sam/SRR769545_bwasw.log

3. BWA-Backtrack algorithm
NOTE: with this algorithm we operate on the forward and reverse reads individually
bwa aln \
  refgenome/GRCh38.p13_ref.fna \
  data/SRR769545_1.fastq.gz \
  > data/SRR769545_1.sai
bwa aln \
  refgenome/GRCh38.p13_ref.fna \
  data/SRR769545_2.fastq.gz \
  > data/SRR769545_2.sai
# use "sampe" to generate the alignments, "samse" is another option. As seen with trimmomatic in the previous section, sampe is for paired-end reads while samse is for single-end reads
bwa sampe \
  refgenome/GRCh38.p13_ref.fna \
  data/SRR769545_1.sai data/SRR769545_2.sai \
  data/SRR769545_1.fastq.gz data/SRR769545_2.fastq.gz \
  > sam/SRR769545_aln.sam 2> sam/SRR769545_aln.log

