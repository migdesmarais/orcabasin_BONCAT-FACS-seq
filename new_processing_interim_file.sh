# Set paths for the GTDB r226 build
REF=/data_store/gtdbtk_db/release226/skani/database
OLD_DB=/scratch/mdesmarais/kraken_gtdb_r226_128g
NEW_DB=/scratch/mdesmarais/kraken_gtdb_r226_128g_v2

mkdir -p "$NEW_DB/library"
printf 'Build log: %s/build.log\n' "$NEW_DB"

# Build the GTDB r226 Kraken2 and Bracken databases
(
set -eo pipefail
conda activate kraken2_env
kraken_bin=$(dirname "$(command -v kraken2)")

test -d "$NEW_DB/library"
test ! -e "$NEW_DB/hash.k2d"
mkdir -p "$NEW_DB/taxonomy"
cp "$OLD_DB"/taxonomy/{names,nodes}.dmp "$NEW_DB/taxonomy/"

python3 -u - "$REF" "$OLD_DB" "$NEW_DB" <<'PY'
import gzip
import sys
from pathlib import Path

reference, old_db, new_db = map(Path, sys.argv[1:])
taxids = {}

with (old_db / "taxonomy/gtdb.accession2taxid").open() as source:
    for line in source:
        _, accession, taxid, *_ = line.split()
        if accession in taxids and taxids[accession] != taxid:
            raise SystemExit(f"Conflicting tax IDs for {accession}")
        taxids[accession] = taxid

fastas = sorted(reference.rglob("*_genomic.fna.gz"))
matched = []
missing = []
suffix = "_genomic.fna.gz"

for fasta in fastas:
    accession = fasta.name[:-len(suffix)]
    if accession in taxids:
        matched.append((fasta, accession))
    else:
        missing.append(accession)

(new_db / "unmapped_assemblies.txt").write_text(
    "".join(f"{accession}\n" for accession in missing)
)
print(f"Staging {len(matched):,} assemblies. Skipping {len(missing):,} without an exact tax ID match.")
if not matched:
    raise SystemExit("No assemblies matched the GTDB tax ID map.")
if len({accession for _, accession in matched}) != len(matched):
    raise SystemExit("Duplicate assembly IDs in the FASTA library.")

with (new_db / "library/gtdb_r226.fna").open("wb") as library, \
     (new_db / "seqid2taxid.map").open("w") as seq_map:
    for i, (fasta, accession) in enumerate(matched, 1):
        taxid = taxids[accession]
        contig = 0
        with gzip.open(fasta, "rb") as reads:
            for line in reads:
                if line.startswith(b">"):
                    contig += 1
                    seqid = f"{accession}_{contig}"
                    library.write(f">{seqid}\n".encode())
                    seq_map.write(f"{seqid}\t{taxid}\n")
                else:
                    library.write(line.rstrip(b"\r\n") + b"\n")
        if i % 5000 == 0 or i == len(matched):
            print(f"Staged {i:,}/{len(matched):,} assemblies")
PY

kraken2-build --build --db "$NEW_DB" --threads 8 \
    --kmer-len 35 --max-db-size 137438953472

test -s "$NEW_DB/hash.k2d"
test -s "$NEW_DB/opts.k2d"
test -s "$NEW_DB/taxo.k2d"

conda activate bracken_env
bracken-build -d "$NEW_DB" -x "$kraken_bin" -t 8 -k 35 -l 150
test -s "$NEW_DB/database150mers.kmer_distrib"
echo "Kraken2 and Bracken builds finished"
) 2>&1 | tee "$NEW_DB/build.log"
