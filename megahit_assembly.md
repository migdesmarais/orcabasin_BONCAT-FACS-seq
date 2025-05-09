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
for r1 in *_paired_R1.fastq.gz; do
  base=$(echo "$r1" | sed 's/_paired_R1.fastq.gz//')
  r2="${base}_paired_R2.fastq.gz"
  outdir="megahit_assemblies/${base}"

  echo "Assembling $base..."
  megahit \
    -1 "$r1" \
    -2 "$r2" \
    -o "$outdir" \
    --min-contig-len 1000 \
    --presets meta-sensitive \
    --num-cpu-threads 20
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



