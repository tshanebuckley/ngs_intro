include { get_fastq } from './modules/get_fastq.nf'
include { get_ref } from './modules/get_ref.nf'
include { index_refgenome } from './modules/index_refgenome.nf'
include { run_alignment } from './modules/run_alignment.nf'
include { process_bam } from './modules/process_bam.nf'
include { detect_variants } from './modules/detect_variants.nf'
include { create_consensus_seq } from './modules/create_consensus_seq.nf'

params.per_id = 'SRR769545'
params.btprefix = 'bowtie2'

workflow {
  // get the fastq files for the target pair-end read
  get_fastq(params.per_id)
  // get the reference genome
  get_ref()
  // index the reference genome
  index_refgenome(get_ref.out.refgenome, params.btprefix)
  // align the reads
  run_alignment(get_fastq.out.file1, get_fastq.out.file2, index_refgenome.out.bt2files, params.btprefix, params.per_id)
  // process the bam file
  process_bam(run_alignment.out, params.per_id)
  // collapse the pileup of reads, create the variant information file (bcf), convert/compress/index the vcf file
  detect_variants(get_ref.out.refgenome, process_bam.out.sorted, params.per_id)
  // create the final consensus sequence
  create_consensus_seq(get_ref.out.refgenome, detect_variants.out.vcfgz, detect_variants.out.index, params.per_id)
}
