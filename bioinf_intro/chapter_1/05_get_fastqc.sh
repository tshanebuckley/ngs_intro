# in this script, we will actually generate the html qc report files
QC_OUTPUT=$PWD/qc
QC_LIMITS=$PWD/limits.txt
mkdir -p qc
cd ./ecoli/fastQC
filenames=$(ls *.fastq)
fastqc --limits $QC_LIMITS $filenames \
  --outdir $QC_OUTPUT \
  --threads 3
