# example script to split the refgenome into chromosomes
cd refgenome
mkdir -p chromosomes
# Use awk to split the FASTA file

# Input genome file (change this to your actual reference genome)
GENOME="GRCh38.p13_ref.fna"

# Output directory
OUTDIR="chromosomes"
mkdir -p "$OUTDIR"

# Check if genome file exists
if [[ ! -f "$GENOME" ]]; then
    echo "Error: Genome file '$GENOME' not found!"

    exit 1
fi

# Split genome by chromosome
awk '
    /^>/ { 
        if (out) close(out);
        # Remove ">" from chromosome name and create a new file
        out = "'"$OUTDIR"'/" substr($0, 2) ".fa"
    }
    { print >> out }

' "$GENOME"
