# Mapping of Sequence Reads to the Reference Genomes

Here, we will focus any pre/post steps that you may want to execute at each script:

- 01;post: you may want to count the number of FASTA sequences in the file:

```bash
grep -c ">" GRCh38.p13_ref.fna
```

Or count the total number of bases:
```bash
grep -v ">" GRCh38.p13_ref.fna | wc | awk '{print $3-$1}'
```

NOTE: script 02 is not from the original book as it seemed to not work. Also note that while files for primary assemblies are created, there are many extra files including Scaffolds, Patches, Fixes, and Novel sequences.

- 03;post: view the fastqc reports:
```bash
openbrowser SRR769545_1_fastq.html SRR769545_2_fastq.html
```

The final script to run in this section is the 'genome_assembly.nf' script in the 'genome_assembly' directory. This is a nexflow script than can be run using the following command from the 'genome_assembly' directory:
```bash
nextflow genome_assembly.nf -resume
```

