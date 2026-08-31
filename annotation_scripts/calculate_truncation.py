import pandas as pd
import numpy as np

def annotate_rmsk(input_file, output_file, threshold=0.9, flank_threshold=0.1):
    """
    Annotate RepeatMasker (UCSC rmsk format) entries as full-length or truncated.

    Parameters
    ----------
    input_file      : path to UCSC rmsk table (plain or .gz)
    output_file     : path for output TSV
    threshold       : min coverage to call Full_Length (default 0.9)
    flank_threshold : min fractional gap to call a truncation at one end (default 0.1)

    """

    cols = [
        'bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns',
        'genoName', 'genoStart', 'genoEnd', 'genoLeft', 'strand',
        'repName', 'repClass', 'repFamily', 'repStart', 'repEnd', 'repLeft', 'id'
    ]

    df = pd.read_csv(input_file, sep='\t', names=cols, compression='infer')

    # Cast to int — repStart and repLeft are negative in UCSC rmsk, preserve sign
    rep_start = df['repStart'].astype(int)
    rep_end   = df['repEnd'].astype(int)
    rep_left  = df['repLeft'].astype(int)

    plus = df['strand'] == '+'

    # --- Consensus length ---
    cons_len = np.where(plus,
                        rep_end + abs(rep_left),   # plus:  repEnd + abs(repLeft)
                        rep_end + abs(rep_start))  # minus: repEnd + abs(repStart)

    # --- Aligned length on consensus ---
    aln_len = np.where(plus,
                   rep_end - rep_start,      # plus strand
                   rep_end - abs(rep_left))  # minus strand (simpler/more direct)

    # --- Coverage ---
    df['cons_len'] = cons_len
    df['aln_len']  = aln_len
    df['coverage'] = np.where(cons_len > 0, aln_len / cons_len, 0.0)

    # --- Full-length call ---
    is_full = df['coverage'] >= threshold
    df['status'] = np.where(is_full, 'Full_Length', 'Truncated')

    # --- Gap sizes (always positive, strand-aware) ---
    # Plus:  5' gap = repStart,       3' gap = abs(repLeft)
    # Minus: 5' gap = abs(repLeft),   3' gap = abs(repStart)
    gap_5p = np.where(plus, rep_start,       abs(rep_left))
    gap_3p = np.where(plus, abs(rep_left),   abs(rep_start))

    # --- Fractional gaps ---
    flank = cons_len * flank_threshold
    trunc_5p = ~is_full & (gap_5p > flank)
    trunc_3p = ~is_full & (gap_3p > flank)

    # --- Annotation (np.select evaluates in order; first match wins) ---
    conditions = [
        is_full,
        trunc_5p & trunc_3p,
        trunc_5p,
        trunc_3p,
    ]
    choices = [
        'Full_Length',
        'Internal_Fragment',
        '5prime_Truncated',
        '3prime_Truncated',
    ]
    # default catches the rare edge case: coverage < threshold but gaps < flank_threshold
    df['annotation'] = np.select(conditions, choices, default='Truncated_Other')

    df.to_csv(output_file, sep='\t', index=False)

    print(f"Annotated {len(df):,} entries → {output_file}")
    print(df['annotation'].value_counts().to_string())


annotate_rmsk('rmsk_mm10.txt.gz', 'annotated_tes_mm10.tsv')