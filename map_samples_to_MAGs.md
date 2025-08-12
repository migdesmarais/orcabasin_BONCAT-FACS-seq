
## Merge MAG
```
mkdir -p renamed_mags

for MAG in *.fa; do
  ID=$(basename "$MAG" .fa)
  awk -v prefix="$ID" '/^>/{print ">" prefix "_" substr($0,2)} !/^>/' "$MAG" > renamed_mags/"$ID.renamed.fa"
done
```
# Then concatenate:
```
cat *.renamed.fa > all_MAGs_unique.fa
```


## Create index
```
conda activate bowtie2_env2
bowtie2-build all_MAGs_unique.fa all_MAGs_index
```

#!/usr/bin/env bash
set -euo pipefail

# --- paths ---
BASE=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ
READS_DIR=$BASE/reads
OUT=$BASE/magmap_out
THREADS=12

INDEX_DIR=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/renamed_mags/all_MAGs
IDX=$INDEX_DIR/all_MAGs_index                # bowtie2 index prefix (no .bt2)
FASTAS=$INDEX_DIR/all_MAGs_unique.fna       # concatenated MAGs FASTA
REF=$FASTAS                                 # use this in calmd

mkdir -p "$OUT"/{logs,bam,counts}

# ensure reference FASTA index exists
[ -s "$REF" ] || { echo "REF FASTA not found: $REF"; exit 1; }
[ -s "${REF}.fai" ] || samtools faidx "$REF"

# --- map all samples ---
for R1 in "$READS_DIR"/*_paired_R1.fastq.gz; do
  [ -e "$R1" ] || { echo "No R1 files found in $READS_DIR"; exit 1; }
  R2=${R1/_paired_R1/_paired_R2}
  [[ -e "$R2" ]] || { echo "Missing mate for $R1"; exit 1; }

  SAMPLE=$(basename "$R1" | sed -E 's/_paired_R1\.fastq\.gz$//' | grep -oE 'OBNC|OB[0-9]+|[A-Za-z0-9._-]+')
  echo "== Mapping $SAMPLE"

  bowtie2 --very-sensitive -p "$THREADS" \
    --no-unal --no-mixed --no-discordant -k 1 \
    -x "$IDX" -1 "$R1" -2 "$R2" \
    2> "$OUT/logs/${SAMPLE}_bowtie2.log" \
  | samtools view -b -q 30 -F 4 -F 256 -F 2048 \
  | samtools sort -@ "$THREADS" -o "$OUT/bam/${SAMPLE}.q30.primary.bam"

  samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"

  # add MD/NM tags so CoverM can enforce %ID
  samtools calmd -bAr "$OUT/bam/${SAMPLE}.q30.primary.bam" "$REF" > "$OUT/bam/${SAMPLE}.tmp.bam"
  mv -f "$OUT/bam/${SAMPLE}.tmp.bam" "$OUT/bam/${SAMPLE}.q30.primary.bam"
  samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"

  samtools idxstats "$OUT/bam/${SAMPLE}.q30.primary.bam" > "$OUT/counts/${SAMPLE}_idxstats.tsv"
done















# set paths (adjust if yours differ)
BASE=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ
READS_DIR=$BASE/reads
OUT=$BASE/magmap_out
THREADS=12

INDEX_DIR=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/renamed_mags/all_MAGs
IDX=$INDEX_DIR/all_MAGs_index              # prefix only, no .bt2
FASTAS=$INDEX_DIR/all_MAGs_unique.fna     # concatenated MAGs fasta

mkdir -p "$OUT"/{logs,bam,counts}

for R1 in "$READS_DIR"/*_paired_R1.fastq.gz; do
  R2=${R1/_R1/_R2}; [[ -e "$R2" ]] || { echo "Missing mate for $R1"; exit 1; }
  SAMPLE=$(basename "$R1" | sed -E 's/_R1\.fastq\.gz$//' | grep -oE 'OBNC|OB[0-9]+|[A-Za-z0-9._-]+')
  echo "== Mapping $SAMPLE"

  bowtie2 --very-sensitive -p "$THREADS" \
    --no-unal --no-mixed --no-discordant -k 1 \
    -x "$IDX" -1 "$R1" -2 "$R2" \
    2> "$OUT/logs/${SAMPLE}_bowtie2.log" \
  | samtools view -b -q 30 -F 4 -F 256 -F 2048 \
  | samtools sort -@ "$THREADS" -o "$OUT/bam/${SAMPLE}.q30.primary.bam"
  samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"

  samtools calmd -bAr "$OUT/bam/${SAMPLE}.q30.primary.bam" "$REF" > "$OUT/bam/${SAMPLE}.tmp.bam"
  mv "$OUT/bam/${SAMPLE}.tmp.bam" "$OUT/bam/${SAMPLE}.q30.primary.bam"
  samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"

  samtools idxstats "$OUT/bam/${SAMPLE}.q30.primary.bam" > "$OUT/counts/${SAMPLE}_idxstats.tsv"
done










# pick one sample to test
R1="$READS_DIR"/OB129_S51_L003_paired_R1.fastq.gz
R2="$READS_DIR"/OB129_S51_L003_paired_R2.fastq.gz
SAMPLE=$(basename "$R1" | sed -E 's/_R1\.fastq\.gz$//' | grep -oE 'OBNC|OB[0-9]+|[A-Za-z0-9._-]+')
echo "Testing $SAMPLE"

# competitive mapping → primary, high-confidence alignments only
bowtie2 --very-sensitive -p "$THREADS" \
  --no-unal --no-mixed --no-discordant \
  -k 1 \
  -x "$IDX" \
  -1 "$R1" -2 "$R2" \
  2> "$OUT/logs/${SAMPLE}_bowtie2.log" \
| samtools view -b -q 30 -F 4 -F 256 -F 2048 \
| samtools sort -@ "$THREADS" -o "$OUT/bam/${SAMPLE}.q30.primary.bam"

samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"

# add MD/NM tags so CoverM’s %ID filter works
samtools calmd -bAr "$OUT/bam/${SAMPLE}.q30.primary.bam" "$FASTAS" > "$OUT/bam/${SAMPLE}.tmp.bam"
mv "$OUT/bam/${SAMPLE}.tmp.bam" "$OUT/bam/${SAMPLE}.q30.primary.bam"
samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"

# quick sanity checks
grep -E "overall alignment rate|No alignments" "$OUT/logs/${SAMPLE}_bowtie2.log"
samtools idxstats "$OUT/bam/${SAMPLE}.q30.primary.bam" | head
















Cool—back to mapping. Here’s a tight, **step-by-step** run you can paste in and execute. It assumes you’ve already concatenated and indexed the MAGs.

# A) Set paths (edit if needed)

```bash
BASE=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ
READS_DIR=$BASE/trimmed_reads                         # <-- your decontaminated/trimmed reads
REN_DIR=$BASE/dereplicated_genomes/renamed_mags/all_MAGs
IDX=$REN_DIR/all_MAGs_index                           # bowtie2 index prefix (no .bt2)
REF=$REN_DIR/all_MAGs_unique.fna                      # concatenated MAGs fasta
OUT=$BASE/magmap_out
THREADS=12
mkdir -p "$OUT"/{logs,bam,counts,coverm}
```

# B) Quick check: does your BAM need MD/NM tags?

```bash
# after mapping (step C) you can test the BAM with:
samtools view -h "$OUT/bam/TEST.bam" | head -200 | grep -m1 -E 'MD:Z|NM:i' || echo "No MD/NM found"
```

If you see neither `MD:Z` nor `NM:i`, run the `samtools calmd` step below.

# C) Map ONE sample (test), filter, add MD/NM

```bash
R1=$(ls "$READS_DIR"/*_R1.fastq.gz | head -n 1)
R2=${R1/_R1/_R2}
SAMPLE=$(basename "$R1" | sed -E 's/_R1\.fastq\.gz$//' | grep -oE 'OBNC|OB[0-9]+|[A-Za-z0-9._-]+')
echo "Testing $SAMPLE"

bowtie2 --very-sensitive -p "$THREADS" \
  --no-unal --no-mixed --no-discordant -k 1 \
  -x "$IDX" -1 "$R1" -2 "$R2" \
  2> "$OUT/logs/${SAMPLE}_bowtie2.log" \
| samtools view -b -q 30 -F 4 -F 256 -F 2048 \
| samtools sort -@ "$THREADS" -o "$OUT/bam/${SAMPLE}.q30.primary.bam"
samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"

# Ensure CoverM can enforce %ID:
samtools calmd -bAr "$OUT/bam/${SAMPLE}.q30.primary.bam" "$REF" > "$OUT/bam/${SAMPLE}.tmp.bam"
mv "$OUT/bam/${SAMPLE}.tmp.bam" "$OUT/bam/${SAMPLE}.q30.primary.bam"
samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"

# Sanity peek
grep -E "overall alignment rate|No alignments" "$OUT/logs/${SAMPLE}_bowtie2.log"
samtools idxstats "$OUT/bam/${SAMPLE}.q30.primary.bam" | head
```

# D) Map ALL samples (same recipe)

```bash
for R1 in "$READS_DIR"/*_R1.fastq.gz; do
  R2=${R1/_R1/_R2}; [[ -e "$R2" ]] || { echo "Missing mate for $R1"; exit 1; }
  SAMPLE=$(basename "$R1" | sed -E 's/_R1\.fastq\.gz$//' | grep -oE 'OBNC|OB[0-9]+|[A-Za-z0-9._-]+')
  echo "== Mapping $SAMPLE"

  bowtie2 --very-sensitive -p "$THREADS" \
    --no-unal --no-mixed --no-discordant -k 1 \
    -x "$IDX" -1 "$R1" -2 "$R2" \
    2> "$OUT/logs/${SAMPLE}_bowtie2.log" \
  | samtools view -b -q 30 -F 4 -F 256 -F 2048 \
  | samtools sort -@ "$THREADS" -o "$OUT/bam/${SAMPLE}.q30.primary.bam"
  samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"

  samtools calmd -bAr "$OUT/bam/${SAMPLE}.q30.primary.bam" "$REF" > "$OUT/bam/${SAMPLE}.tmp.bam"
  mv "$OUT/bam/${SAMPLE}.tmp.bam" "$OUT/bam/${SAMPLE}.q30.primary.bam"
  samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"

  samtools idxstats "$OUT/bam/${SAMPLE}.q30.primary.bam" > "$OUT/counts/${SAMPLE}_idxstats.tsv"
done
```

# E) Compute per-MAG metrics with CoverM

```bash
coverm genome \
  --bam-files "$OUT"/bam/*.q30.primary.bam \
  --genome-fasta-directory "$(dirname "$REN_DIR")" \
  --methods relative_abundance tpm rpkm coverage breadth \
  --min-read-percent-identity 95 \
  --min-read-aligned-percent 75 \
  --contig-end-exclusion 75 \
  --exclude-supplementary-alignments \
  --threads "$THREADS" \
  --output-file "$OUT/coverm/coverm_metrics.tsv"

head -n 5 "$OUT/coverm/coverm_metrics.tsv"
```

# F) Call “active” per MAG × depth (breadth + TPM)

```bash
awk -v OFS="\t" -v B=0.15 -v T=1 '
BEGIN{print "sample","genome","breadth","coverage","TPM","relative_abundance","active"}
NR>1{
  s=$1; g=$2; ra=$3; tpm=$4; cov=$5; br=$6;
  print s,g,br,cov,tpm,ra,(br>=B && tpm>=T?"yes":"no")
}' "$OUT/coverm/coverm_metrics.tsv" > "$OUT/coverm/active_calls.tsv"

head -n 10 "$OUT/coverm/active_calls.tsv"
```

---

### Quick checkpoints (paste results if anything looks off)

* Bowtie2 log: overall alignment rate not crazy low/high (a few–tens of % is common).
* `samtools view -h … | grep MD:` shows MD/NM present after `calmd`.
* CoverM file has sensible **breadth** (0–1) and **TPM** (>0 where breadth is decent).

When you finish step C for one sample, tell me the alignment rate line + the first few rows of `coverm_metrics.tsv`, and I’ll sanity-tune thresholds if needed.































Here’s a **2-week, no-nonsense sprint plan** to get a strong “Dark Microbial Matter” poster out the door using exactly what you’re already doing (decontam → map to MAGs). It’s optimized to **only deep-analyze the MAGs that are actually active**, so you don’t sink time into all 400.

# Strategy (TL;DR)

1. **Finish mapping** trimmed+decontaminated reads → dereplicated MAGs (competitive, strict filters).
2. **Call active MAGs by depth** using **breadth ≥0.15 & TPM ≥1** (tweak later).
3. **Focus all downstream work on just those active MAGs** (usually tens, not hundreds).
4. Do **fast taxonomy** on the active set only: GTDB-Tk → if unclassified, CAT/BAT + AAI to nearest refs.
5. Make 4–5 clean figures and 6–8 bullet points that tell the story.

---

# Day-by-day (what to actually do)

### Days 1–3: Mapping & Active calls

* Run the mapping pipeline we set up (you’re already there): Bowtie2 (`-k1`, MAPQ≥30) → `samtools calmd` → **CoverM** (`breadth, TPM, coverage`, `--min-read-percent-identity 95`, `--min-read-aligned-percent 75`).
* Produce `coverm_metrics.tsv` and `active_calls.tsv`.

**Quick selector (keeps only what you’ll analyze further):**

```bash
# Active per (sample, MAG)
awk 'NR>1 && $7=="yes"{print $1,$2}' OFS='\t' active_calls.tsv \
 | sort -u > active_pairs.tsv

# Unique active MAG list
cut -f2 active_pairs.tsv | sort -u > active_MAGs.txt
```

### Days 3–5: Minimal taxonomy on **active MAGs only**

1. **GTDB-Tk v2** on the active set:

```bash
mkdir active_bins && xargs -I{} cp MAGS/{}.fa active_bins/ < active_MAGs.txt
gtdbtk classify_wf --genome_dir active_bins --out_dir gtdbtk_out --cpus 16
```

2. For MAGs still “unclassified” or low confidence:

   * **CAT/BAT** (bin LCA across all proteins) — great fallback.
   * **AAI** to nearest refs (CompareM AAI or MMseqs2 ezAAI) to assign **family/order level** when genus/species isn’t possible.

**AAI quickie (pairwise to GTDB reps list):**

```bash
comparem aai_wf -x fa active_bins aai_out -T 16
# summarize each MAG to its top reference AAI
```

> Poster text: report the **highest supported rank** across GTDB-Tk / CAT-LCA / AAI bands; call the rest “novel at genus/family/order”.

### Day 5: QC (fast, only on active MAGs)

* **CheckM2** for completeness/contam; **GUNC** for chimeras.
* Drop obvious chimeras; keep notes for methods box.

### Days 6–7: Community context (reads level)

* Run **Kraken2/Bracken** (already underway) on the same decontaminated libraries → get order/family profiles by depth.
* This gives you the “who’s there” panel to pair with “who’s active (MAG breadth)”.

### Days 7–10: Figures (lock these)

1. **Depth line**: % BONCAT-active cells vs depth (your FACS result).
2. **Stacked bars**: Bracken family (or order) composition by depth.
3. **Heatmap**: **MAG breadth** by depth (rows=active MAGs, columns=depths; binarize ≥0.15 for “active”).
4. **Novelty plot**: bar/pie showing fraction of active MAGs that are **novel at genus/family/order** (from GTDB-Tk + AAI bands + CAT).
5. **Mini QC table**: for top \~15 headline MAGs (length, completeness/contam, breadth max, TPM max, taxonomy call, AAI to nearest ref, GUNC pass/fail).

> Optional if time allows: tiny phylogenomic tree of the **active set only** with rings for depth of activity + novelty rank. (Skip if it risks the deadline.)

### Days 10–12: Story + methods boxes

* **Methods box** (bullets): BONCAT-FACS→MDA; competitive mapping (Bowtie2, `-k1`, MAPQ≥30), CoverM with **95% ID & 75% aligned**, presence by **breadth ≥0.15 + TPM ≥1**, GTDB-Tk → CAT/BAT → AAI; CheckM2 & GUNC; interpret abundance qualitatively due to MDA.
* **Results bullets** (the sound bites):

  * “Detected **N active MAGs** across **D depths**; activity peaks at the particulate-rich halocline.”
  * “**X%** of active MAGs are **novel at or above family level** (AAI bands + GTDB-Tk unclassified).”
  * “Active lineages include **\[named groups]** plus **previously unclassified clades** unique to Orca brine/interface.”
  * “Community profiles (Bracken) **agree** with MAG recruitment at family/order level.”

### Days 12–14: Polish

* Tighten fonts/legends; add **one schematic** of the workflow (BONCAT-FACS → MDA → mapping → active calls → taxonomy stack).
* Add a one-line **limitations** note (MDA bias; breadth-based detection).
* Add a **QR code** to a GitHub gist with methods/commands (optional but slick).

---

## Fast taxonomy calls (when time is tight)

Use this ladder for each active MAG; **stop** as soon as you get a confident rank:

1. **GTDB-Tk v2** → if genus/species assigned with reasonable placement, use it.
2. If not: **CAT/BAT** → take **family/order** LCA if consistent across proteins.
3. **AAI band** to nearest GTDB ref to label **novel genus/family/order** (e.g., “58% AAI to Family X → likely novel genus within X”).
4. If 16S exists (barrnap): corroborate with **SILVA SINA**; don’t over-weight it.

---

## Decision rules (so you don’t over-think)

* **Active call** per depth: `breadth ≥ 0.15` **and** `TPM ≥ 1` after 95% ID & 75% aligned filters.
* **Report taxonomy** at the **highest rank that’s consistent** across GTDB-Tk / CAT / AAI.
* **Novelty bins**:

  * **Novel genus**: AAI \~65–95% to nearest ref **and** not placed to named genus.
  * **Novel family**: AAI \~45–65%.
  * **Novel order**: AAI \~35–45%.
* **Drop** any active MAG that **fails GUNC** or has **contam >10%** (or flag with an asterisk).

---

## Minimal commands you’ll likely need next

**Select active set and copy FASTAs**

```bash
mkdir -p active_bins
while read m; do cp MAGS/"$m".fa active_bins/; done < active_MAGs.txt
```

**GTDB-Tk on active set**

```bash
gtdbtk classify_wf --genome_dir active_bins --out_dir gtdbtk_out --cpus 16
```

**CAT/BAT (fallback)**

```bash
CAT bins -b active_bins -d CAT_DB/ -t GTDB_tax/ -o cat_bins.tsv -n 16
```

**AAI summary**

```bash
comparem aai_wf -x fa active_bins aai_out -T 16
# then parse each MAG’s top reference and its AAI to assign genus/family/order band
```

---

## Poster map (what goes where)

* **Title**: “Dark microbial matter in a deep-sea brine: BONCAT-FACS read recruitment reveals novel active lineages in Orca Basin”
* **Left column**: schematic + methods bullets.
* **Center top**: % BONCAT-active vs depth.
* **Center**: MAG breadth heatmap (active only).
* **Right top**: Kraken/Bracken stacked bars by depth.
* **Right middle**: novelty bar (genus/family/order novel fractions).
* **Bottom**: mini table of headline MAGs + limitations note.

---

If you want, I can:

* give you a **tiny R snippet** that takes `coverm_metrics.tsv` → breadth heatmap + novelty bar, or
* a one-pager **methods text** you can paste straight onto the poster.























## MAP
```
mkdir -p /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/MAG_mapping_logs

for R1 in /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/clean_reads/*_clean_R1.fastq.gz; do
  R2=${R1/_R1/_R2}
  SAMPLE=$(basename "$R1" | grep -oE 'OB[0-9]+|OBNC')

  echo "Mapping $SAMPLE reads to MAGs..."

  bowtie2 --very-sensitive -p 12 \
    -x /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/renamed_mags/all_MAGs_index \
    -1 "$R1" -2 "$R2" \
    -S "/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/MAG_mapping_sam/${SAMPLE}_vs_MAGs.sam" \
    2> "/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/MAG_mapping_logs/${SAMPLE}_bowtie2.log"
done
```

## Conver SAM to BAM
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

## MAG abundance - coverM
```
coverm genome \
  --bam-files /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/MAG_mapping_sam/*.bam \
  --genome-fasta-directory /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/renamed_mags \
  --threads 12 \
  --methods relative_abundance \
  --min-read-percent-identity 95 \
  --min-read-aligned-percent 75 \
  --output-file coverm_abundance.tsv

coverm genome \
  --bam-files /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/MAG_mapping_sam/*.bam \
  --genome-fasta-directory /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/renamed_mags \
  --threads 12 \
  --methods relative_abundance \
  --output-file coverm_abundance.tsv
```





