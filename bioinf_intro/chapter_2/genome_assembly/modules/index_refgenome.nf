process index_refgenome {

  input:
    path refgenome
    val btprefix

  output:
    path "${btprefix}.*.bt2", emit: bt2files

  script:
  """
  bowtie2-build ${refgenome} ${btprefix}
  """

}
