process process_bam {

  input:
    path bam_file
    val pre_id

  output:
    path "${pre_id}_sorted.bam", emit: sorted
    path "${pre_id}_sorted.bam.bai", emit: index

  script:
  """
  samtools sort -T tmp.sort -o ${pre_id}_sorted.bam ${bam_file}
  samtools index ${pre_id}_sorted.bam
  """
  
}
