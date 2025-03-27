# run the trimmomatic preprocessing
cd trimmomatic
trimmomatic \
  PE SRR957824_1.fastq SRR957824_2.fastq \
  out_PE_SRR957824_1_min150.fastq out_UPE_SRR957824_1_min150.fastq \
  out_PE_SRR957824_2_min150.fastq out_UPE_SRR957824_2_min150.fastq \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True \
  LEADING:3 \
  TRAILING:3 \
  ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 \
  SLIDINGWINDOW:5:30 \
  MINLEN:150
# generate the fastqc reports
fastqc out_PE_SRR957824_1_min150.fastq out_PE_SRR957824_2_min150.fastq
