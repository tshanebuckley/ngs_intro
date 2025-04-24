process detect_variants {

  input:
    path refgenome
    path bamfile
    val per_id

  output:
    path "${per_id}.bcf", emit: bcf
    path "${per_id}.vcf", emit: vcf
    path "${per_id}.vcf.gz", emit: vcfgz
    path "${per_id}.vcf.gz.tbi", emit: index

  script:
  """
  # collapse the pileup of reads aligned in to the reference genome and write out the variant info
  bcftools mpileup -f ${refgenome} \
    ${bamfile} \
    | bcftools call -mv -Ob \
    -o ${per_id}.bcf
  # convert the binary variant call file to a vcf file
  bcftools convert -O v \
    -o ${per_id}.vcf \
    ${per_id}.bcf
  # compress and index the vcf file
  bgzip -c ${per_id}.vcf > ${per_id}.vcf.gz
  tabix -p vcf ${per_id}.vcf.gz
  """

}
