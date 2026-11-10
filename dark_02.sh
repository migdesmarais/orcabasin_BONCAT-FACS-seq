# run hhmer
conda create -n dark_oxygen_hmmer -c conda-forge -c bioconda hmmer -y
conda activate dark_oxygen_hmmer

# chcekc hmm
ls HCO/HMM/*.hmm | wc -l

for hmm in HCO/HMM/*.hmm; do
    echo "$(basename "$hmm"): $(grep -m 1 '^NAME ' "$hmm")"
done


# run hhmer
# ============================================================
# ORCA BASIN — DARK OXYGEN
# Screen proteins against Emeline's complete HMM collection
# ============================================================

mkdir -p hmm_results

for hmm in HCO/HMM/*.hmm; do

    model=$(basename "$hmm" .hmm)

    echo "========================================"
    echo "Searching model: $model"
    echo "========================================"

    hmmsearch \
        --cpu 8 \
        -E 10 \
        --domE 10 \
        --tblout "hmm_results/${model}.tblout" \
        --domtblout "hmm_results/${model}.domtblout" \
        -o "hmm_results/${model}.log" \
        "$hmm" \
        orca_MT_reference_proteins.faa

done

echo "ALL HMM SEARCHES COMPLETE"


# prepapre combined fasta
mkdir -p phylogeny

cat \
  HCO/Sequence_database/HCO-family-sequences_unique/NOD_sequences_unique.txt \
  HCO/Sequence_database/HCO-family-sequences_unique/QNOR_sequences_unique.txt \
  NOD_all_68_proteins.faa \
  > phylogeny/NOD_QNOR_ORCA_combined.faa

# then run iqtree, lost command



# NLD
#?
mkdir -p phylogeny

(
  set -euo pipefail
  for pdb in 3Q08 3NN1 3QPI 4WWS; do
    curl -fsSL "https://www.rcsb.org/fasta/entry/$pdb/display" \
      -o "phylogeny/${pdb}.fasta"
    awk -v label="REF_${pdb}" '
      /^>/ { if (++n > 1) exit; print ">" label; next }
      n == 1 { print }
    ' "phylogeny/${pdb}.fasta" > "phylogeny/${pdb}_one.faa"
    grep -q '^[A-Za-z]' "phylogeny/${pdb}_one.faa"
  done

  cat phylogeny/{3Q08,3NN1,3QPI,4WWS}_one.faa \
    CLD_all_hit_proteins.faa > phylogeny/CLD_ORCA_combined.faa

  grep -c '^>' phylogeny/CLD_ORCA_combined.faa
  mafft --auto --thread 8 phylogeny/CLD_ORCA_combined.faa \
    > phylogeny/CLD_ORCA_aligned.faa
)

cd /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dark_oxygen
conda activate dark_oxygen_hmmer
command -v iqtree2
conda install -n dark_oxygen_hmmer -c conda-forge -c bioconda iqtree -y &&
conda activate dark_oxygen_hmmer &&
cd /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dark_oxygen &&
IQTREE=$(command -v iqtree2 || command -v iqtree) &&
"$IQTREE" -s phylogeny/CLD_ORCA_aligned.faa \
  -m MFP -B 1000 --alrt 1000 -T 8 \
  --prefix phylogeny/CLD_ORCA
