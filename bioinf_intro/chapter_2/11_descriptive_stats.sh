# in this script, we will use samtools to generate some summary statistics of our bam file
mkdir -p stats
samtools flagstat processed/SRR769545_mem_sorted.bam > stats/flag.txt
samtools coverage processed/SRR769545_mem_sorted.bam > stats/coverage.txt
samtools depth processed/SRR769545_mem_sorted.bam > stats/depth.txt
