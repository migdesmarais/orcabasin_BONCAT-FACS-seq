## Install SingleM
```
curl -fsSL https://pixi.sh/install.sh | sh
git clone https://github.com/wwood/singlem
cd singlem
pixi shell
singlem -h
lyrebird -h

singlem data --output-directory singlem_db
singlem list-markers --refpkg singlem_db

export SINGLEM_METAPACKAGE_PATH=$/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/singlem/singlem_db/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb
```

## Assign taxonomy to BONCAT-FACS-seq data
```
singlem pipe -1 /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/megahit_assemblies/contigs/*.fa \
 --taxonomic-profile /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/megahit_assemblies/contigs/singlem_output/boncat_taxonomic_profile.tsv \
 --taxonomic-profile-krona /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/megahit_assemblies/contigs/singlem_output/boncat_taxonomic_profile.krona.tsv  \
 --otu-table /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/megahit_assemblies/contigs/singlem_output/boncat_otu_table.tsv  \
 --threads 16 \
 --metapackage /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/singlem/singlem_db/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb
```

```
singlem pipe -1 /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/megahit_assemblies/contigs/OBNC_contigs.fa \
 --taxonomic-profile /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/megahit_assemblies/contigs/singlem_output_OBNC/boncat_taxonomic_profile.tsv \
 --taxonomic-profile-krona /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/megahit_assemblies/contigs/singlem_output_OBNC/krona_OBNC.html  \
 --otu-table /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/megahit_assemblies/contigs/singlem_output_OBNC/boncat_otu_table.tsv  \
 --threads 16 \
 --metapackage /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/singlem/singlem_db/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb
```

## Assign taxonomy to dereplicated genomes
```
singlem pipe -1 /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/*.fa \
 --taxonomic-profile /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/singlem_output/boncat_taxonomic_profile.tsv \
 --taxonomic-profile-krona /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/singlem_output/boncat_taxonomic_profile.krona.tsv  \
 --otu-table /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/singlem_output/boncat_otu_table.tsv  \
 --threads 16 \
 --metapackage /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/singlem/singlem_db/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb
```

## Visualize krona from BONCAT seqs
```
conda create -n krona_env -c bioconda krona
conda activate krona_env



```


