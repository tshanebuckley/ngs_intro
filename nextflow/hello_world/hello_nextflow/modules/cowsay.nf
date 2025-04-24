#!/usr/bin/env nextflow

// Generate ASCII art with cowpy
process cowsay {

  publishDir 'results', mode: 'copy'
  container 'grycap/cowsay'

  input:
    path input_file
    val character

  output:
    path "cowsay-${input_file}"

  script:
  """
  cat $input_file | /usr/games/cowsay -f "$character" > cowsay-${input_file}
  """

}

