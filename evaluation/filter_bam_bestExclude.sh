bamfile="young_pseudobulk_nonZeroXP_primary_updated.bam"

# Pass 1: find max XP per read
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pass 1: Finding max XP per read..."
samtools view -@ 8  $bamfile | \
awk '{
    match($0, /XP:f:([0-9.eE+\-]+)/, arr)
    xp = arr[1] + 0
    if (!(($1) in max_xp) || xp > max_xp[$1])
        max_xp[$1] = xp
}
END {
    for (qname in max_xp)
        print qname, max_xp[qname]
}' > max_xp_per_read.txt
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pass 1 complete: $(wc -l < max_xp_per_read.txt) reads processed"

# Pass 2: count ties
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pass 2: Counting ties..."
samtools view -@ 8 $bamfile | \
awk 'BEGIN{
    while ((getline line < "max_xp_per_read.txt") > 0) {
        split(line, a, " ")
        max_xp[a[1]] = a[2]
    }
}
{
    match($0, /XP:f:([0-9.eE+\-]+)/, arr)
    xp = arr[1] + 0
    if (xp == max_xp[$1])
        tie_count[$1]++
}
END {
    for (qname in tie_count)
        print qname, tie_count[qname]
}' > tie_count_per_read.txt
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pass 2 complete: $(awk '$2 > 1' tie_count_per_read.txt | wc -l) tied reads found"

# Pass 3: keep only reads with exactly one max XP (no ties)
samtools view -@ 8 -h $bamfile | \
awk 'BEGIN{
    while ((getline line < "max_xp_per_read.txt") > 0) {
        split(line, a, " ")
        max_xp[a[1]] = a[2]
    }
    while ((getline line < "tie_count_per_read.txt") > 0) {
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
samtools view -@ 8 -b -o young_pseudobulk_bestExclude.bam

