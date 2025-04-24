# in this script, we will fetch the rna data we will need to run STAR
mkdir -p ucscref
wget \
  -O "ucscref/hg38.fa.gz" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz"
wget \
  -O "ucscref/hg38.ncbiRefSeq.gtf.gz" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz"
# unzip the data
gzip -d ucscref/hg38.fa.gz
gzip -d ucscref/hg38.ncbiRefSeq.gtf.gz
