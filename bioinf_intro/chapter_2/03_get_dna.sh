# simple script to fetch data as before
mkdir -p data
cd data
fasterq-dump --verbose SRR769545
# compress the files
gzip SRR769545_1.fastq
gzip SRR769545_2.fastq

# this file includes some bootstrap code to get the fastqc config from github
FASTQC_GH="https://raw.githubusercontent.com/s-andrews/FastQC/master/Configuration"

FASTQC_ADAPTER=$FASTQC_GH/adapter_list.txt
FASTQC_CONTAMINANT=$FASTQC_GH/contaminant_list.txt

FASTQC_LIMITS=$FASTQC_GH/limits.txt

wget $FASTQC_ADAPTER
wget $FASTQC_CONTAMINANT
wget $FASTQC_LIMITS

# generate the fastqc reports
fastqc --limits limits.txt SRR769545_1.fastq.gz SRR769545_2.fastq.gz
