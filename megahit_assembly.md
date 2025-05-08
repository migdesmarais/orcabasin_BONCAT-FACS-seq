## Install MEGAHIT


## Make directory
mkdir -p megahit_assemblies

## Make script
nano run_megahit_all.sh

## Script:

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

## Make it executable
chmod +x run_megahit_all.sh

## Run it
./run_megahit_all.sh
