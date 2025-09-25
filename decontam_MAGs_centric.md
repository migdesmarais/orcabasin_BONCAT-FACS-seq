# Conda install
# Concatenate reads (Samples and NC, separately)
```
# where your reads live
BASE=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ
READS=$BASE/reads
OUT=$BASE/coassembly_megahit
mkdir -p "$OUT"
cd "$READS"

# ---- build the two merged libraries ----
# ALL (everything paired, excluding NC_)
ls *_paired_R1.fastq.gz | grep -v '^NC_' > "$OUT/all_R1.list"
sed 's/_R1/_R2/' "$OUT/all_R1.list" > "$OUT/all_R2.list"
cat $(cat "$OUT/all_R1.list") > "$OUT/ALL_R1.fastq.gz"
cat $(cat "$OUT/all_R2.list") > "$OUT/ALL_R2.fastq.gz"

# NC controls
ls OBNC_*_paired_R1.fastq.gz > "$OUT/nc_R1.list"
sed 's/_R1/_R2/' "$OUT/nc_R1.list" > "$OUT/nc_R2.list"
cat $(cat "$OUT/nc_R1.list") > "$OUT/NC_R1.fastq.gz"
cat $(cat "$OUT/nc_R2.list") > "$OUT/NC_R2.fastq.gz"

echo "[done] Merged to:"
ls -lh "$OUT"/ALL_* "$OUT"/NC_*
```

# Coassemble with megahit
```
conda activate megahit
BASE=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ
ASM=$BASE/coassembly_megahit
QC=$BASE/coassembly_megahit/qc
READS=$BASE/reads
mkdir -p $QC
mkdir -p "$ASM"

megahit -1 "$OUT/ALL_R1.fastq.gz" -2 "$OUT/ALL_R2.fastq.gz" -o "$ASM/ALL" -t 24 --min-contig-len 1000
megahit -1 "$OUT/NC_R1.fastq.gz"  -2 "$OUT/NC_R2.fastq.gz"  -o "$ASM/NC"  -t 24 --min-contig-len 1000











conda activate quast_env

# NC
quast.py -t 16 -m 1000 -o $QC/quast_NC  $ASM/NC/final.contigs.fa

grep -E 'contigs \(>= 1000 bp\)|Total length|Largest contig|N50|L50|GC %' \
  $QC/quast_NC/report.txt

# contigs (>= 1000 bp)      70           
Total length (>= 0 bp)      490538       
Total length (>= 1000 bp)   490538       
Total length (>= 5000 bp)   411806       
Total length (>= 10000 bp)  321225       
Total length (>= 25000 bp)  77085        
Total length (>= 50000 bp)  0            
Largest contig              49808        
Total length                490538       
N50                         12586        
L50                         13      

# (when ready) ALL
quast.py -t 16 -m 1000 -o $QC/quast_ALL $ASM/ALL/final.contigs.fa

grep -E 'contigs \(>= 1000 bp\)|Total length|Largest contig|N50|L50|GC %' \
  $QC/quast_ALL/report.txt

```

# Map + bin samples + NC
# NC
```
conda activate mags

# 0) Clean any partial BAMs
BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
READS=$BASE/trimmed_reads
ASM=$BASE/assemblies
RUN=$BASE/mag_results/NC
REF=$ASM/NC/final.contigs.fa
mkdir -p $RUN/map
rm -f $RUN/map/*.bam $RUN/map/*.bai

# 1) Rebuild file lists with FULL paths
ls -1 $READS/NC_*_paired_R1.fastq.gz > $RUN/NC_R1.list
sed 's/_R1/_R2/' $RUN/NC_R1.list > $RUN/NC_R2.list
head $RUN/NC_R1.list   # should show /scratch/... full paths

# 2) Index (message “read 0 ALT contigs” is normal)
bwa index "$REF"

# 3) Map each NC library
paste $RUN/NC_R1.list $RUN/NC_R2.list | while read -r R1 R2; do
  S=$(basename "$R1" | sed 's/_L003.*//')
  echo "[*] Mapping $S"
  bwa mem -t 16 "$REF" "$R1" "$R2" | samtools sort -@8 -o $RUN/map/${S}.bam -
  samtools index $RUN/map/${S}.bam
done

# Sanity check
ls -lh $RUN/map/*.bam
samtools idxstats $RUN/map/*.bam | head

# mapping yield (per-BAM), how many reads mapped?
for B in $BASE/mag_results/NC/map/*.bam; do
  echo "== $(basename "$B") =="
  samtools flagstat "$B" | sed -n '1,8p'
done

# depths
jgi_summarize_bam_contig_depths --outputDepth $RUN/map/depth.txt $RUN/map/*.bam
head -n 3 $RUN/map/depth.txt

# Average contig depth
DEP=$BASE/mag_results/NC/map/depth.txt
awk 'NR>1 {n++; d+=$3} END{print "mean_totalAvgDepth =", (n? d/n : 0)}' "$DEP"

# THE ASSEMBLY IS VERY POOR, SO WHAT WE WILL DO NEXT IS JUST CLASSIFY READS WITH KRAKEN2 ONCE THE DB IS DOWNLOADED (GTDB-TK) AND THEN USE THE NC TAXONOMY TO REMOVE MAGS THAT MATCH IN THE SAMPLES "ALL".
