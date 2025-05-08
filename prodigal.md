## Install PRODIGAL

conda create --name prodigal
conda activate prodigal
conda install -c bioconda prodigal
prodigal -v

## Run PRODIGAL on MEGAHIT assemblies

mkdir -p /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/prodigal_calls

for f in /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/megahit_assemblies/contigs/*.fa; do
    base=$(basename "$f" .fa)
    prodigal -i "$f" \
             -a /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/prodigal_calls/${base}_ORFs.faa \
             -d /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/prodigal_calls/${base}_nucleotides.fna \
             -o /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/prodigal_calls/${base}.gbk \
             -p meta \
             -f gbk
    done


