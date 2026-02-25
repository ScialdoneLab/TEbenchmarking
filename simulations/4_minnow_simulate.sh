#!/bin/bash

echo "Old"
src/minnow simulate --splatter-mode \
    --g2t data_mouse/mm10_GENCODE_rmsk_TE.annotation.tx2gene.tsv \
    --inputdir splatter_old_mm \
    --PCR 4 \
    -r minnow_index_mm10_TE/ref_k101_fixed.fa \
    -e 0.01 -p 24 \
    -o minnow_simulate_mm_RA_TE_old \
    --dbg --gfa minnow_index_mm10_TE/dbg.gfa \
    -w ../data/737K-august-2016.txt \
    --countProb minnow_estimate_mm_RA_TE/countProb.txt --custom --gencode
echo "Done with old"


echo "Young"
src/minnow simulate --splatter-mode \
    --g2t data_mouse/mm10_GENCODE_rmsk_TE.annotation.tx2gene.tsv \
    --inputdir splatter_young_mm \
    --PCR 4 \
    -r minnow_index_mm10_RA_TE/ref_k101_fixed.fa \
    -e 0.01 -p 24 \
    -o minnow_simulate_mm_RA_TE_young \
    --dbg --gfa minnow_index_mm10_TE/dbg.gfa \
    -w ../data/737K-august-2016.txt \
    --countProb minnow_estimate_mm_RA_TE/countProb.txt --custom --gencode
echo "Done with young"


echo "All"
src/minnow simulate --splatter-mode \
    --g2t data_mouse/mm10_GENCODE_rmsk_TE.annotation.tx2gene.tsv \
    --inputdir splatter_all_mm \
    --PCR 4 \
    -r minnow_index_mm10_TE/ref_k101_fixed.fa \
    -e 0.01 -p 24 \
    -o minnow_simulate_mm_RA_TE_all \
    --dbg --gfa minnow_index_mm10_TE/dbg.gfa \
    -w ../data/737K-august-2016.txt \
    --countProb minnow_estimate_mm_RA_TE/countProb.txt --custom --gencode
echo "Done with all"


echo "TEs + Genes"
src/minnow simulate --splatter-mode \
    --g2t data_mouse/mm10_GENCODE_rmsk_TE_wGene.annotation.tx2gene.tsv \
    --inputdir splatter_all_wGenes_mm \
    --PCR 4 \
    -r minnow_index_mm10_RA_TEsNgenes/ref_k101_fixed.fa \
    -e 0.01 -p 24 \
    -o minnow_simulate_mm_RA_TE_all_TEsNgenes \
    --dbg --gfa minnow_index_mm10_RA_TEsNgenes/dbg.gfa \
    -w ../data/737K-august-2016.txt \
    --countProb minnow_estimate_mm10_RA_TEsNgenes/countProb.txt --custom --gencode
    
echo "Done with TEs + Genes"