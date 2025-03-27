# we will re-fetch the unmodified data to be preprocessed
mkdir -p trimmomatic && cd trimmomatic
fasterq-dump SRR957824
# we also need to fetch our adapters from github
TRIMMOMATIC_GH="https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters"

wget $TRIMMOMATIC_GH/TruSeq3-PE.fa
wget $TRIMMOMATIC_GH/TruSeq2-PE.fa
