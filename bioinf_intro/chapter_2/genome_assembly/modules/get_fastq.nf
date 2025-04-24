process get_fastq {

  input:
    //paired-end reads id
    val per_id

  output:
    path "${per_id}_1.fastq.gz", emit: file1
    path "${per_id}_2.fastq.gz", emit: file2

  script:
  """
  fasterq-dump --verbose ${per_id}
  gzip ${per_id}_1.fastq
  gzip ${per_id}_2.fastq
  """
}
