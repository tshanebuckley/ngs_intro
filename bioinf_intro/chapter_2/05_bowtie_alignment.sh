# in this script, we will use bowtie to index and align our reference genome
mkdir -p refgenome/bowtie2
# index using bowtie
bowtie2-build \
  --threads 2 \
  refgenome/GRCh38.p13_ref.fna \
  refgenome/bowtie2
# run the alignment using bowtie2
bowtie2 -x refgenome/bowtie2 \
  -1 data/SRR769545_1.fastq.gz \
  -2 data/SRR769545_2.fastq.gz \
  -S sam/SRR769545_bowtie2.sam
