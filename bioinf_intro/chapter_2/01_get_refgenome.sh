# this script is used to fetch the reference human genome from the ncbi
mkdir -p refgenome
# we will make this sam directory for later scripts as well
mkdir -p sam
wget -O "refgenome/GRCh38.p13_ref.fna.gz" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
cd refgenome
# unzip the data
gunzip -d GRCh38.p13_ref.fna.gz
# use samtools to index the genome
samtools faidx GRCh38.p13_ref.fna
