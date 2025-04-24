process get_ref {

  output:
    path "hg38.fa", emit: refgenome

  script:
  """
  wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip -d hg38.fa.gz
  """
  
}
