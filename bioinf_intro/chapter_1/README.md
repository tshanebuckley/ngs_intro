# CH1: Sequencing and Raw Sequence Data Quality Control

For the most part, you may can execute this chapter by running through the scripts.

Here, we will focus any pre/post steps that you may want to execute at each script:

- 01;post: you may want to run the following command to preview the first few lines of the file:
```bash
head -15 SRR030834.fastq
```

- 02;post: you may want to also preview our generated fasta file:
```bash
head -15 SRR030834.fasta
```

- 05;pre: you should update the limits.txt file:
```
set "kmer ignore 1" to "kmer ignore 0"
```

- 05;post: navigate to the "qc" directory and use the "openbrowser" command to open any of the html files generated from our run of fastqc:
```bash
openbrowser SRR576933_fastqc.html
```

NOTE: script 05 and above are a simply intro of fetching fastq files and running fastqc report generation. The book referenced explains the significance of and how to read each graph generated.

After 05, we move into a basic example of preprocessing SRR957824.

- 06;post: open the qc report for the bad fastq file:
```bash
openbrowser bad_fastqc.html
```

- 07;post: open the filtered fastqc report:
```bash
openbrowser filtered_fastqc.html
```

- 08;post: open the trimmed fastqc report:
```bash
openbrowser trimmed_fastqc.html
```

- 09;post: open the clipped fastqc report:
```bash
openbrowser clipped_fastqc.html
```

- 10;post: open the len150 filtered fastqc report:
```bash
openbrowser len150_fastqc.html
```

NOTE: 11 is an alternative approach to script 10 that masks the low quality reads instead of actually removing the ones that are low quality. Script 10 is what is known as hard-masking, whereas script 11 is referred to as soft-masking.

- 11;post: open the masked fastqc report:
```bash
openbrowser masked_fastqc.html
```

NOTE: the last two scripts are exampled of using an alternative to fastx-tools called trimmomatic.

NOTE: SRR957824_1.fastq is the forward read and SRR957824_1.fastq is the reverse read. When using both of these reads, we may refer to this as the "two paired-end FASTQ files". Trimmomatic has a "PE" signature for paired-end reads and "SE" for single-end reads.

NOTE: script 13 is simply removing any reads shorter than 35 base pairs while script 14 uses the same hard-masking strategy as before by only keeping reads of the max length.

- 13;post: open the results of the trimmomatic preprocessing:
```bash
openbrowser out_PE_SRR957824_1_min35_fastqc.html out_PE_SRR957824_2_min35_fastqc.html
```

- 14;post: open the results of the trimmomatic preprocessing with a hard-mask applied:
```bash
openbrowser out_PE_SRR957824_1_min150_fastqc.html out_PE_SRR957824_2_min150_fastqc.html
```

