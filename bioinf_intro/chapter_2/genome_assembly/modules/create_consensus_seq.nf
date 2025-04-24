process create_consensus_seq {

  input:
    path refgenome
    path vcfgz
    path index
    val per_id

  output:
    path "${per_id}_genome.fasta", emit: genome

  script:
  """
  bcftools consensus -f ${refgenome} ${vcfgz} > ${per_id}_genome.fasta
  """
    
}
