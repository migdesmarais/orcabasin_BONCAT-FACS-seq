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
```

## Assign taxonomy to BONCAT-FACS-seq data
```
singlem search \
  --sequences /path/to/your/contigs_or_MAGs/ \
  --output_directory /path/to/output_folder/ \
  --threads 16 \
  --singlem_package /path/to/marker_package.smpkg
```

