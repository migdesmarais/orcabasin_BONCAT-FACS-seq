## Merging files from two sequencing runs
### For R2 files, change R*

```
cat /data_store/seq_data/igm_20250307_miguel_cat/OB161_S45_L003_R2_001.fastq.gz /data_store/seq_data/igm_20250502/OB161_S99_L001_R2_001.fastq.gz > /data_store/seq_data/igm_20250307_miguel_cat/merged_OB/OB161_S45_L003_R2_001.fastq.gz
```

## Check on number of reads per file
```
for f in OB*.fastq.gz; do
  reads=$(zcat "$f" | wc -l | awk '{printf "%.0f", $1/4}')
  echo "$f: $reads reads"
done
```

## Create derivative files to my /scratch directory
```
ln -s /data_store/seq_data/igm_20250307_miguel_cat/merged_OB/* /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/raw_data/
```

## Install graftM on conda environment
```
git clone https://github.com/geronimp/graftM
cd graftM
conda env create -n graftM -f graftm.yml
conda activate graftM
cd bin
export PATH=$PWD:$PATH
graftM -h
```

## Install trimmomatic
```
conda create --name trimmomatic_env 
conda activate trimmomatic_env
conda install -c bioconda trimmomatic
trimmomatic -version
```

## Install fastqc
```
conda create --name fastqc
conda activate fastqc
conda install -c bioconda fastqc
fastqc -version
```

## Run pre-QC fastqc
```
fastqc *fastq.gz
```

##	Download and visualize QC files
```
scp -r mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/raw_reads/pre_fastqc ~/Downloads/
```

##	Trimmomatic

```
conda activate trimmomatic_env
```

```
conda activate trimmomatic_env

# Make sure output directory exists
mkdir -p /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads

# Run Trimmomcatic
for file in $(cat samples.txt); do
  trimmomatic PE -phred33 -threads 12 \
    ${file}_R1_001.fastq.gz ${file}_R2_001.fastq.gz \
    /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads/${file}_paired_R1.fastq.gz \
    /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads/${file}_unpaired_R1.fastq.gz \
    /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads/${file}_paired_R2.fastq.gz \
    /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads/${file}_unpaired_R2.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```

## Run post-QC fastqc
```
fastqc *fastq.gz
```

##	Download and visualize QC files
```
scp -r mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads/post_fastqc ~/Downloads/
```

## Use NC to decontaminate reads DID NOT WORK
```
conda create -n bbtools_env -c agbiome bbtools -y
conda activate bbtools_env

cat OBNC_S54_L003_*.fastq.gz > OBNC_combined.fastq.gz

# check sizes
ls -lh OBNC_*.fastq.gz

for f in OBNC_S54_L003_*.fastq.gz; do
    echo -n "$f: "
    echo $(( $(zcat "$f" | wc -l) / 4 ))
done

# create kmers
kmercountexact.sh in=OBNC_combined.fastq.gz out=OBNC_kmers.fa k=31

#use kmers to decontam

mkdir -p ../clean_reads
mkdir -p ../clean_reads/stats

for r1 in *_paired_R1.fastq.gz; do
  r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
  base=$(basename "$r1" _paired_R1.fastq.gz)

  bbduk.sh \
    in1="$r1" in2="$r2" \
    out1="../clean_reads/${base}_clean_R1.fastq.gz" \
    out2="../clean_reads/${base}_clean_R2.fastq.gz" \
    ref=OBNC_kmers.fa \
    k=31 hdist=1 mink=11 ktrim=r ftm=5 \
    stats="../clean_reads/stats/${base}_bbduk_stats.txt" \
    threads=8
done
```

##	Calculate how many reads are leftover in each paired fastq files
```
echo -e "Filename\tReads" > paired_read_counts.tsv
for f in *_paired_*.fastq.gz; do
  count=$(zcat "$f" | wc -l | awk '{print $1/4}')
  echo -e "$f\t$count" >> paired_read_counts.tsv
done
```

## Remove contaminants found in OBNC from the other reads and place in clean_reads folder
## Map samples reads to OBNC contigs
```
#need to clean header and produce index
bowtie2-build OBNC_contigs.fa OBNC_clean_index

for R1 in *_paired_R1.fastq.gz; do
  R2=${R1/_R1/_R2}
  SAMPLE=$(echo "$R1" | grep -oE 'OB[0-9]+|OBNC')

  echo "Mapping sample: $SAMPLE"

  bowtie2 --very-sensitive -p 12 \
    -x /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/megahit_assemblies/contigs/OBNC_clean_index \
    -1 "$R1" -2 "$R2" \
    -S "${SAMPLE}_OBNC.sam" \
    2> "logs/${SAMPLE}_bowtie2.log"
done
```


## Convert SAM to BAM files
```
for SAM in *_OBNC.sam; do
  SAMPLE=$(echo "$SAM" | grep -oE 'OB[0-9]+|OBNC')
  echo "Converting $SAMPLE to sorted BAM..."

  samtools view -bS "$SAM" | samtools sort -o "${SAMPLE}_OBNC.sorted.bam"
  samtools index "${SAMPLE}_OBNC.sorted.bam"
done
```

## Extract mapped (contaminants) and unmapped (clean) reads
```
for BAM in *_OBNC.sorted.bam; do
  SAMPLE=$(echo "$BAM" | grep -oE 'OB[0-9]+|OBNC')
  echo "Extracting reads for $SAMPLE..."

  # Clean reads: both mates unmapped
  samtools view -b -f 12 -F 256 "$BAM" > "${SAMPLE}_unmapped.bam"

  # Contaminant reads: both mates mapped
  samtools view -b -f 3 "$BAM" > "${SAMPLE}_mapped_contaminants.bam"
done
```


