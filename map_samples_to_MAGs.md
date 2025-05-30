
## Merge MAGs
```
cat *.fa > all_MAGs_combined.fa
```

## Create index
```
bowtie2-build all_MAGs_combined.fa all_MAGs_index
```

## MAP
```
for R1 in /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads/*_paired_R1.fastq.gz; do
  R2=${R1/_R1/_R2}
  SAMPLE=$(basename "$R1" | grep -oE 'OB[0-9]+|OBNC')

  echo "Mapping $SAMPLE reads to MAGs..."

  bowtie2 --very-sensitive -p 12 \
    -x /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/all_MAGs_index \
    -1 "$R1" -2 "$R2" \
    -S "/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/MAG_mapping_sam/${SAMPLE}_vs_MAGs.sam"
done
```

## Conver SAM to BAM?
```
for SAM in /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/MAG_mapping_sam/*.sam; do
  BAM=${SAM%.sam}.bam
  samtools view -bS "$SAM" | samtools sort -o "$BAM"
  samtools index "$BAM"
done
```

## Count reads
```
for BAM in *.bam; do
  SAMPLE=$(basename "$BAM" | grep -oE 'OB[0-9]+|OBNC')
  samtools idxstats "$BAM" > "${SAMPLE}_MAG_counts.tsv"
done
```

## MAG abundance
```
coverm genome -1 clean_R1.fastq.gz -2 clean_R2.fastq.gz \
  --genome-fasta-directory MAGs/ --threads 12 --output-file coverm_abundance.tsv

```
