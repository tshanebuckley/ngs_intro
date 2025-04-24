process run_alignment {

  input:
    path read1
    path read2
    path bt2files
    val btprefix
    val pre_id

  output:
    path "${pre_id}.bam"

  script:
  """
  bowtie2 -x ${btprefix} -1 ${read1} -2 ${read2} | samtools view -Sb - > ${pre_id}.bam
  """

}
