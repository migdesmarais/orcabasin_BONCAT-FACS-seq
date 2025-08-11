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

## Check taxonomy of reads
```
mkdir -p kraken
shopt -s nullglob

# get unique sample IDs without the lane/_paired bits
samples=$(ls *_paired_R1.fastq.gz | sed -E 's/_L[0-9]{3}_paired_R1\.fastq\.gz$//; s/_paired_R1\.fastq\.gz$//' | sort -u)

for s in $samples; do
  echo ">> Combining $s"

  # gather all R1/R2 lane files for this sample (e.g., L001, L002, L003…)
  r1=(${s}_L*_paired_R1.fastq.gz ${s}_paired_R1.fastq.gz)
  r2=(${s}_L*_paired_R2.fastq.gz ${s}_paired_R2.fastq.gz)

  # combine (concatenate gzip streams safely)
  cat "${r1[@]}" > "kraken/${s}_R1.combined.fastq.gz"
  cat "${r2[@]}" > "kraken/${s}_R2.combined.fastq.gz"
done





conda activate kraken_env

for fq in *_combined.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)

    kraken2 --db /data_store/kraken_database \
            --gzip-compressed \
            --output - \
            --use-names \
            "$fq" | sed "s/^/${sample}\t/" >> all_samples.kraken
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


## Convert SAM to BAM files and index manually
```
samtools view -@ 12 -bS OBNC_OBNC.sam | samtools sort -@ 12 -o OBNC_OBNC.sorted.bam
samtools index OBNC_OBNC.sorted.bam
```


## Extract mapped (contaminants) and unmapped (clean) reads + summary
```
echo -e "Sample\tTotal_Reads\tClean_Reads\tContaminant_Reads" > /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/clean_reads/decontamination_summary.tsv

for bam in /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads/*_OBNC.sorted.bam; do
    sample=$(basename "$bam" | cut -d'_' -f1)

    R1=$(ls /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads/${sample}_S*_L003_paired_R1.fastq.gz)
    R2=$(ls /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads/${sample}_S*_L003_paired_R2.fastq.gz)

    echo "Processing $sample"

    samtools view -F 4 "$bam" | cut -f 1 | sort | uniq > /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/clean_reads/${sample}_contaminant_ids.txt

    filterbyname.sh in="$R1" in2="$R2" \
        out=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/clean_reads/${sample}_contam_R1.fastq.gz \
        out2=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/clean_reads/${sample}_contam_R2.fastq.gz \
        names=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/clean_reads/${sample}_contaminant_ids.txt include=t overwrite=t

    filterbyname.sh in="$R1" in2="$R2" \
        out=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/clean_reads/${sample}_clean_R1.fastq.gz \
        out2=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/clean_reads/${sample}_clean_R2.fastq.gz \
        names=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/clean_reads/${sample}_contaminant_ids.txt include=f overwrite=t

    total=$(zcat "$R1" | echo $((`wc -l` / 4)))
    clean=$(zcat /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/clean_reads/${sample}_clean_R1.fastq.gz | echo $((`wc -l` / 4)))
    removed=$((total - clean))

    echo -e "${sample}\t${total}\t${clean}\t${removed}" >> /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/clean_reads/decontamination_summary.tsv

    echo "  Total: $total | Clean: $clean | Removed: $removed"
done
```

## Check taxonomy of reads/contaminants
```
mkdir -p kraken

for idfile in *_contaminant_ids.txt; do
    sample=$(basename "$idfile" _contaminant_ids.txt)

    # Combine clean R1 and R2
    cat ${sample}_clean_R1.fastq.gz ${sample}_clean_R2.fastq.gz > kraken/${sample}_clean_combined.fastq.gz

    # Combine contaminant R1 and R2
    cat ${sample}_contam_R1.fastq.gz ${sample}_contam_R2.fastq.gz > kraken/${sample}_contam_combined.fastq.gz
done


conda activate kraken_env

for fq in *_combined.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)

    kraken2 --db /data_store/kraken_database \
            --gzip-compressed \
            --output - \
            --use-names \
            "$fq" | sed "s/^/${sample}\t/" >> all_samples.kraken
done

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


