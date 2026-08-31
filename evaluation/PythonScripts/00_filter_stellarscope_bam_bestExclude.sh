#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 -p PREFIX -d WORKDIR [-s SUFFIX]"
    echo "  -p PREFIX   Sample prefix, e.g. 'old' (required)"
    echo "  -d WORKDIR  Directory containing input BAM and where outputs are written (required)"
    echo "  -s SUFFIX   Input BAM suffix to strip for basename (default: '-updated.bam')"
    exit 1
}

SUFFIX="-updated.bam"

while getopts "p:d:s:h" opt; do
    case $opt in
        p) PREFIX="$OPTARG" ;;
        d) WORKDIR="$OPTARG" ;;
        s) SUFFIX="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

: "${PREFIX:?Must specify -p PREFIX}"
: "${WORKDIR:?Must specify -d WORKDIR}"

# Make sure we're using gawk, since match(str, re, arr) is a gawk extension
if ! awk --version 2>/dev/null | grep -qi gnu; then
    echo "ERROR: This script requires gawk (GNU awk) for the 3-argument match() function." >&2
    echo "Your default 'awk' does not appear to be gawk. Install gawk or adjust PATH." >&2
    exit 1
fi

cd "$WORKDIR"

INPUT_BAM="${PREFIX}_pseudobulk${SUFFIX}"
BASENAME=$(basename "$INPUT_BAM" "$SUFFIX")

if [[ ! -f "$INPUT_BAM" ]]; then
    echo "ERROR: Input BAM not found: $WORKDIR/$INPUT_BAM" >&2
    exit 1
fi

# Derived filenames
PRIMARY_BAM="${BASENAME}_primary_updated_vScript.bam"
NONZERO_BAM="${BASENAME}_nonZeroXP_primary_updated_vScript.bam"
BESTEXCLUDE_BAM="${BASENAME}_bestExclude_vScript.bam"
MAX_XP_FILE="${BASENAME}_max_xp_per_read.txt"
TIE_COUNT_FILE="${BASENAME}_tie_count_per_read.txt"

# Clean up intermediate text files on exit (success or failure)
trap 'rm -f "$MAX_XP_FILE" "$TIE_COUNT_FILE"' EXIT

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "Working dir: $WORKDIR"
log "Input BAM: $INPUT_BAM"
log "Basename: $BASENAME"

# Step 1: Keep only PRI alignments
log "Step 1: Keeping PRI alignments..."
samtools view -@ 8 -h "$INPUT_BAM" | \
awk 'BEGIN{OFS="\t"} /^@/ {print; next} $0 ~ /ZT:Z:PRI/ {print}' | \
samtools view -@ 8 -b -o "$PRIMARY_BAM"
log "Step 1 complete: $(samtools view -c "$PRIMARY_BAM") PRI reads"

# Step 2: Filter out zero posteriors
log "Step 2: Filtering out zero posteriors..."
samtools view -@ 8 -h "$PRIMARY_BAM" | \
awk 'BEGIN{OFS="\t"} /^@/ {print; next} $0 !~ /XP:f:0(\.0*)?(\t|$)/ {print}' | \
samtools view -@ 8 -b -o "$NONZERO_BAM"
log "Step 2 complete: $(samtools view -c "$NONZERO_BAM") reads with non-zero posterior"

# Step 3: find max XP per read AND tie counts in a single pass
log "Step 3: Finding max XP and tie counts per read..."
samtools view -@ 8 "$NONZERO_BAM" | \
awk '{
    match($0, /XP:f:([0-9.eE+\-]+)/, arr)
    xp = arr[1] + 0
    qname = $1

    count[qname, xp]++            # tally how many alignments have this exact XP for this read
    if (!(qname in max_xp) || xp > max_xp[qname]) {
        max_xp[qname] = xp
    }
    if (!(qname in seen)) {
        seen[qname] = 1
        qnames[++nq] = qname       # remember read names in order of first appearance
    }
}
END {
    for (i = 1; i <= nq; i++) {
        q = qnames[i]
        print q, max_xp[q] > "'"$MAX_XP_FILE"'"
        print q, count[q, max_xp[q]] > "'"$TIE_COUNT_FILE"'"
    }
}'
log "Step 3 complete: $(wc -l < "$MAX_XP_FILE") reads processed, $(awk '$2 > 1' "$TIE_COUNT_FILE" | wc -l) tied reads found"

# Step 4: Filter to bestExclude
log "Step 4: Writing bestExclude BAM..."
samtools view -@ 8 -h "$NONZERO_BAM" | \
awk -v maxfile="$MAX_XP_FILE" -v tiefile="$TIE_COUNT_FILE" 'BEGIN{
    while ((getline line < maxfile) > 0) {
        split(line, a, " ")
        max_xp[a[1]] = a[2]
    }
    while ((getline line < tiefile) > 0) {
        split(line, a, " ")
        tie_count[a[1]] = a[2]
    }
}
/^@/ { print; next }
{
    match($0, /XP:f:([0-9.eE+\-]+)/, arr)
    xp = arr[1] + 0
    if (xp == max_xp[$1] && tie_count[$1] == 1) print
}' | \
samtools view -@ 8 -b -o "$BESTEXCLUDE_BAM"
log "Step 4 complete: $(samtools view -c "$BESTEXCLUDE_BAM") reads in final output"

log "Pipeline complete!"