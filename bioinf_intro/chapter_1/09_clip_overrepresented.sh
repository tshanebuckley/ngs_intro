# clips (removes) the overrepresented sequences
# specifically removes the adaptor noted by the -a flag
cd preprocessing
# runs the clipping
fastx_clipper \
  -a ATCGGGAGAGGGGCGGGGAGGGGAAGAGGGGAGAATTCGGGGGGGGCCGG \
  -i trimmed.fastq \
  -o clipped.fastq \
  -v \
  -Q33
# generated the fastqc report on the clipped data
fastqc --limits ../limits.txt clipped.fastq
