
## Predict ORFs on both in situ/PS BONCAT-active MAGs
```
# A) predict genes (CDS) on each MAG

mkdir -p /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/active_mags_b20_insitu/ref/cds

screen
conda activate prodigal

for fa in *.fna; do
  base=$(basename "$fa" .fna)
  prodigal -i "$fa" \
           -a "ref/cds/${base}.faa" \
           -d "ref/cds/${base}.fna" \
           -p meta -q
done

mkdir -p /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/active_mags_b20_PS/ref/cds

for fa in *.fna; do
  base=$(basename "$fa" .fna)
  prodigal -i "$fa" \
           -a "ref/cds/${base}.faa" \
           -d "ref/cds/${base}.fna" \
           -p meta -q
done

# We make not making a nonredundant gene catalog because we are interested in reads mapped to each MAG.
```

## Make a combined CDS reference (keep MAG ID in headers)
```
cd /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/active_mags_b20_insitu

# Put the MAG base name in every CDS header (safe, idempotent)
mkdir -p ref/cds_renamed
for f in ref/cds/*.fna; do
  b=$(basename "$f" .fna)
  awk -v p="$b" 'BEGIN{OFS=""}
    /^>/{sub(/^>/,">"p"|"); print; next}
         {print}' "$f" > "ref/cds_renamed/${b}.fna"
done

# Concatenate to one reference
cat ref/cds_renamed/*.fna > ref/all_cds_MAGcontext.fna

cd /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/active_mags_b20_PS

# Put the MAG base name in every CDS header (safe, idempotent)
mkdir -p ref/cds_renamed
for f in ref/cds/*.fna; do
  b=$(basename "$f" .fna)
  awk -v p="$b" 'BEGIN{OFS=""}
    /^>/{sub(/^>/,">"p"|"); print; next}
         {print}' "$f" > "ref/cds_renamed/${b}.fna"
done

# Concatenate to one reference
cat ref/cds_renamed/*.fna > ref/all_cds_MAGcontext.fna
```

## Build a Salmon index on CDS
```
screen
conda activate salmon_env
cd /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/active_mags_b20_insitu

salmon index \
  -t ref/all_cds_MAGcontext.fna \
  -i ref/salmon_idx_cds \
  -k 31
```

