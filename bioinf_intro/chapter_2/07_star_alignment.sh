# in this script, we will use STAR to index and align an RNA-Seq reads
mkdir -p ucscref/star
# run the indexing
STAR --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir ucscref/star \
  --genomeFastaFiles ucscref/hg38.fa \
  --sjdbGTFfile ucscref/hg38.ncbiRefSeq.gtf \
  --sjdbOverhang 100
# run the alignment
STAR --runThreadN 4 \
  --runMode alignReads \
  --readFilesCommand zcat \
  --genomeDir ucscref/star \
  --outFileNamePrefix sam/SRR769545_star \
  --readFilesIn data/SRR769545_1.fastq.gz data/SRR769545_2.fastq.gz
