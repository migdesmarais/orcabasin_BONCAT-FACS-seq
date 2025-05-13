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

##	Calculate how many reads are leftover in each paried fastq files
```
echo -e "Filename\tReads" > paired_read_counts.tsv
for f in *_paired_*.fastq.gz; do
  count=$(zcat "$f" | wc -l | awk '{print $1/4}')
  echo -e "$f\t$count" >> paired_read_counts.tsv
done
```
