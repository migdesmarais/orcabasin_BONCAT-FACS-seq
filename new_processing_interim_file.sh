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

# rerun Kraken2 on all with new database
screen -S OB_kraken2_v2

conda activate kraken2_env
cd /scratch/mdesmarais/OB_BONCAT-FACS-SEQ

DB=/scratch/mdesmarais/kraken_gtdb_r226_128g_v2
OUT=kraken2_results_r226_v2
T0_IN=

for d in T0_reads/trimmed_20M T0_reads/trimmed; do
    if [[ -d "$d" ]] && compgen -G "$d/*_R1.fastq.gz" >/dev/null; then
        T0_IN="$d"
        break
    fi
done

if [[ -z "$T0_IN" ]]; then
    echo "Can't find trimmed T0 read pairs. Stopping before running Kraken2."
    exit 1
fi

for f in hash.k2d opts.k2d taxo.k2d database150mers.kmer_distrib; do
    if [[ ! -s "$DB/$f" ]]; then
        echo "Missing database file: $DB/$f"
        exit 1
    fi
done

mkdir -p "$OUT"/{T0,ALL,BON}/{reports,logs}
shopt -s nullglob

T0_R1=("$T0_IN"/*_R1.fastq.gz)
ALL_R1=(mapping_reads/*_ALL_R1.fastq.gz)
BON_R1=(mapping_reads/*_BON_R1.fastq.gz)

echo "T0 input: $T0_IN — ${#T0_R1[@]} pairs"
echo "ALL: ${#ALL_R1[@]} pairs"
echo "BON: ${#BON_R1[@]} pairs"
echo "Reports and logs: $OUT"

run_kraken () {
    local r1="$1"
    local group="$2"
    local r2="${r1%_R1.fastq.gz}_R2.fastq.gz"
    local sample="${r1##*/}"
    sample="${sample%_R1.fastq.gz}"

    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "Missing mate for $r1; stopping."
        return 1
    fi

    echo "Starting $group / $sample — $(date)"
    kraken2 --db "$DB" --threads 8 --paired "$r1" "$r2" \
        --report "$OUT/$group/reports/$sample.kreport" \
        --output /dev/null 2>&1 | tee "$OUT/$group/logs/$sample.log"
    echo "Finished $group / $sample — $(date)"
}

for r1 in "${T0_R1[@]}"; do run_kraken "$r1" T0; done
for r1 in "${ALL_R1[@]}"; do run_kraken "$r1" ALL; done
for r1 in "${BON_R1[@]}"; do run_kraken "$r1" BON; done

echo "All Kraken2 samples finished — $(date)"

cd /scratch/mdesmarais/OB_BONCAT-FACS-SEQ
conda activate bracken_env

DB=/scratch/mdesmarais/kraken_gtdb_r226_128g_v2
ROOT=/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/kraken2_results_r226_v2

test -s "$DB/database150mers.kmer_distrib" || { echo "Bracken distribution file missing"; exit 1; }

mkdir -p "$ROOT"/bracken/{genus,phylum,logs}

echo "Kraken reports by group:"
for group in T0 ALL BON; do
    count=$(find "$ROOT/$group/reports" -maxdepth 1 -type f -name '*.kreport' 2>/dev/null | wc -l)
    echo "$group: $count reports"
done

set -euo pipefail

for group in T0 ALL BON; do
    for report in "$ROOT/$group/reports/"*.kreport; do
        [[ -f "$report" ]] || continue
        sample=$(basename "$report" .kreport)

        echo "Bracken genus: $group / $sample"
        bracken -d "$DB" -i "$report" \
            -o "$ROOT/bracken/genus/${sample}.bracken.tsv" \
            -r 150 -l G -t 10 \
            2>&1 | tee "$ROOT/bracken/logs/${group}_${sample}_genus.log"

        echo "Bracken phylum: $group / $sample"
        bracken -d "$DB" -i "$report" \
            -o "$ROOT/bracken/phylum/${sample}.bracken.tsv" \
            -r 150 -l P -t 10 \
            2>&1 | tee "$ROOT/bracken/logs/${group}_${sample}_phylum.log"

        echo "Finished: $group / $sample"
    done
done

echo "All Bracken runs finished."
echo "Genus outputs:  $(find "$ROOT/bracken/genus" -type f -name '*.tsv' | wc -l)"
echo "Phylum outputs: $(find "$ROOT/bracken/phylum" -type f -name '*.tsv' | wc -l)"

