# first, we will create an ecoli directory 
mkdir -p ecoli && cd ecoli
# create a text file to hold the IDs of our SRA files
touch ids.txt && chmod 775 ids.txt
echo "SRR653520" >> ids.txt
echo "SRR653521" >> ids.txt
echo "SRR576933" >> ids.txt
echo "SRR576934" >> ids.txt
echo "SRR576935" >> ids.txt
echo "SRR576936" >> ids.txt
echo "SRR576937" >> ids.txt
echo "SRR576938" >> ids.txt
# now, we download the fastqc data
mkdir fastQC
while read f;
do
  fasterq-dump \
    --outdir fastQC "$f" \
    --progress \
    --threads 4
done<ids.txt
