## {upload in situ dereplicated MAGs from Rebecca to server}
```
scp -r /Users/migueldesmarais/Downloads/dereplicated_genomes mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/OB_BONCAT-FACS-SEQ
```

## CheckM on them and assign taxonomy with graftM, make table
```
checkm lineage_wf -x fa /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/checkm -t 20
```

## then collate checkM resuslt to metadata

