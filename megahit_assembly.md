## Install MEGAHIT
```
conda create --name megahit
conda activate megahit
conda install -c bioconda megahit
megahit --version
```
## Make directory
```
mkdir -p megahit_assemblies
```

## Make script
```
nano run_megahit_all.sh
```

## Script:
```
#!/usr/bin/env bash
set -euo pipefail

# where to write assemblies
OUTROOT="megahit"
mkdir -p "$OUTROOT"

# tweak these if you like
THREADS=20
MINLEN=1000       # min contig length to report
PRESET="meta-sensitive"

shopt -s nullglob
for R1 in *_paired_R1.fastq.gz; do
  base="${R1%_paired_R1.fastq.gz}"
  R2="${base}_paired_R2.fastq.gz"
  outdir="${OUTROOT}/${base}"
  log="${outdir}.log"

  # sanity: make sure R2 exists
  if [[ ! -f "$R2" ]]; then
    echo "[SKIP] $base – missing $R2" | tee -a "$log"
    continue
  fi

  # resume if assembly already done
  if [[ -s "${outdir}/final.contigs.fa" ]]; then
    echo "[DONE] $base – ${outdir}/final.contigs.fa exists" | tee -a "$log"
    continue
  fi

  echo "[RUN ] Assembling $base → $outdir" | tee "$log"
  megahit \
    -1 "$R1" \
    -2 "$R2" \
    -o "$outdir" \
    --min-contig-len "$MINLEN" \
    --presets "$PRESET" \
    --num-cpu-threads "$THREADS" \
    >>"$log" 2>&1

  # quick post-check
  if [[ -s "${outdir}/final.contigs.fa" ]]; then
    echo "[OK  ] $base – contigs: $(grep -c '^>' "${outdir}/final.contigs.fa")" | tee -a "$log"
  else
    echo "[FAIL] $base – contigs file missing/empty" | tee -a "$log"
  fi
done
```
## Make it executable
```
chmod +x run_megahit_all.sh
```
## Run it
```
./run_megahit_all.sh
```
## Check assembly quality with quast
```
conda create --name quast_env
conda activate quast_env
conda install -c bioconda quast
quast --version
```

```
quast -o QC *_contigs.fa
scp -r migdesmarais@fram.ucsd.edu:/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads/megahit_assemblies/contigs/QC /Users/migueldesmarais/Downloads
```



