## {upload in situ dereplicated MAGs from Rebecca to server}
```
scp -r /Users/migueldesmarais/Downloads/dereplicated_genomes mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/OB_BONCAT-FACS-SEQ
```

## CheckM on them and assign taxonomy with graftM, make table
```
checkm lineage_wf -x fa /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/checkm -t 20
```

## Run PRODIGAL on dereplicated_genomes from Rebecca at Stanford
```
mkdir -p /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/prodigal_calls

for f in /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/*.fa; do
    base=$(basename "$f" .fa)
    prodigal -i "$f" \
             -a /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/prodigal_calls/${base}_ORFs.faa \
             -d /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/prodigal_calls/${base}_nucleotides.fna \
             -o /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/prodigal_calls/${base}.gbk \
             -p meta \
             -f gbk
    done
```

