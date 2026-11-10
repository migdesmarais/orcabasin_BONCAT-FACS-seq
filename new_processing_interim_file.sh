# Set paths for the GTDB r226 build
REF=/data_store/gtdbtk_db/release226/skani/database
OLD_DB=/scratch/mdesmarais/kraken_gtdb_r226_128g
NEW_DB=/scratch/mdesmarais/kraken_gtdb_r226_128g_v2

mkdir -p "$NEW_DB/library"
printf 'Build log: %s/build.log\n' "$NEW_DB"


