# make a preprocessing directory to load our data into
mkdir -p preprocessing && cd preprocessing
fasterq-dump --verbose SRR957824
# remove file 1 as we are focusing on file 2
rm SRR957824_1.fastq
# rename file 2
mv SRR957824_2.fastq bad.fastq
# run fastq report generation for the file
fqfile=$(ls *.fastq)
fastqc --limits ../limits.txt $fqfile
