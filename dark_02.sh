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

