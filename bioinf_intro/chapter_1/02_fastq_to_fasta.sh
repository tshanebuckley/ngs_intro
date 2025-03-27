# move the the directory containing the file
cd fastqs/single
# convert to a tabular format
cat SRR030834.fastq  | paste - - - - \
  >SRR030834_tab.tmp
# remove the '@' symbol
awk '{print $1 "\t" $4}' SRR030834_tab.tmp \
  | sed 's/@//g'>SRR030834_seq.tmp
# add a '>' to the beginning of each line
sed -i 's/^/>/' SRR030834_seq.tmp
# separate the defline and sequence, creatin the fasta file
awk '{print $1, "\n" $2}' SRR030834_seq.tmp \
  >SRR030834.fasta
# remove any tmp files created above
rm *.tmp
