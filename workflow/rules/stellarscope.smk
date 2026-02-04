
# Run stellarscope 

## Set up

rule setup_stellarscope:
    """
    setup stellarscope:
        - sort the bam file by barcodes
    """
    output:
        RESULTS_PATH_STAR + "/results/STARoutdir/{dataset}/{sample_id}/score_Aligned.sortedByCB.bam"
    input:
        bam=RESULTS_PATH_STAR+"/results/STARoutdir/{dataset}/{sample_id}/score_Aligned.sortedByCoord.out.bam",
        barcodes=RESULTS_PATH_STAR+"/results/STARoutdir/{dataset}/{sample_id}/score_Solo.out/Gene/raw/barcodes.tsv"
    log:
        "logs/stellarscope/{dataset}_{sample_id}_setup.log"
    params:
        genome_id=get_genome,
        results_path=RESULTS_PATH_STAR
    threads: 12
    resources:
        mem_gb=48
    benchmark:
        "benchmarks/setup_stellarscope_{dataset}_{sample_id}.txt"
    conda: 
        "../envs/stellarscope.yml"
    shell:
        """
        mkdir -p {params.results_path}/results/stellarscope/
        mkdir -p {params.results_path}/results/stellarscope/tmp/

        stellarscope cellsort \
            --nproc {threads} \
            --tempdir {params.results_path}/results/stellarscope/tmp/ \
            --outfile {params.results_path}/results/STARoutdir/{wildcards.dataset}/{wildcards.sample_id}/score_Aligned.sortedByCB.bam \
        	{input.bam} {input.barcodes} 2> {log}
        """

# rule make_minnow_unstranded:
#     """
#     transform minnow-derived, stellarscope-sorted bam files to be unstranded
#     """
#     output:
#         bam=RESULTS_PATH_STAR + "/results/STARoutdir/{dataset}/{sample_id}/score_Aligned.changedStrand.bam"
#     input:
#         bam=RESULTS_PATH_STAR + "/results/STARoutdir/{dataset}/{sample_id}/score_Aligned.sortedByCB.bam",
#         gtf=get_TEannotation_dataset
#     log:
#         "logs/stellarscope/{dataset}_{sample_id}_changeStrands.log"
#     params:
#         genome_id=get_genome,
#         results_path=RESULTS_PATH_STAR
#     threads: 12
#     resources:
#         mem_gb=48
#     benchmark:
#         "benchmarks/stellarscope_{dataset}_{sample_id}_changeStrands.txt"
#     conda: 
#         "../envs/pyBamEnv.yml"
#     shell:
#         "python workflow/scripts/changeStrands.py {input.bam} {input.gtf} {output.bam}"

rule unzip_TEannotation:
    """
    Unzip TE rmsk gtf file
    """
    output:
        TEgtf_unzipped="tmp/TEannotation_{genome_id}.gtf"
    input:
        TEannotation=get_TEannotation_genome
    log:
        "logs/unzip_TEannotation_{genome_id}.txt"
    threads: 3
    resources:
        mem_gb=16
    benchmark:
        "benchmark/unzip_TEannotation_{genome_id}.txt"
    conda:
        "../envs/stellarscope.yml"
    shell:
        "gunzip -c {input.TEannotation} > {output.TEgtf_unzipped}"



rule run_stellarscope_pseudobulk_onestep:
    """
    Run stellarscope in one command
    """
    output:
        outdir=directory(RESULTS_PATH_Stellarscope + "/results/stellarscope_out/{dataset}_pseudobulk_onestep/{sample_id}/"),
        matrix=RESULTS_PATH_Stellarscope + "/results/stellarscope_out/{dataset}_pseudobulk_onestep/{sample_id}/{sample_id}_pseudobulk-TE_counts.mtx",
        updated_bam=RESULTS_PATH_Stellarscope + "/results/stellarscope_out/{dataset}_pseudobulk_onestep/{sample_id}/{sample_id}_pseudobulk-updated.bam"
    input:
        bam=RESULTS_PATH_STAR + "/results/STARoutdir/{dataset}/{sample_id}/score_Aligned.sortedByCB.bam",
        barcodes=get_whitelist,
        TEannotation=get_TEannotation_dataset_unzipped
    log:
        "logs/stellarscope/stellarscope_stload_onestep_{dataset}_{sample_id}.txt"
    params:
        results_path=RESULTS_PATH_Stellarscope
    threads: 24
    resources:
        mem_gb=250
    benchmark:
        "benchmarks/run_stellarscope_onestep_{dataset}_{sample_id}.txt"
    conda: 
        "../envs/stellarscope.yml"
    shell:
        """
        # create a directory 
        mkdir -p "{params.results_path}/results/stellarscope_out"
        mkdir -p "{params.results_path}/results/stellarscope_out/{wildcards.dataset}_pseudobulk_onestep"
        mkdir -p "{params.results_path}/results/stellarscope_out/{wildcards.dataset}_pseudobulk_onestep/{wildcards.sample_id}/"
        
        outfolder={output.outdir}

        # stellarscope pseudobulk onestep
        stellarscope assign \
            --exp_tag {wildcards.sample_id}_pseudobulk \
            --whitelist <(gunzip -c {input.barcodes}) \
            --attribute transcript_id \
            --outdir $outfolder \
            --nproc {threads} \
            --pooling_mode pseudobulk \
            --reassign_mode best_exclude \
            --logfile logs/stellarscope/{wildcards.dataset}_{wildcards.sample_id}_onestep_pseudobulk.txt \
            --max_iter 500 \
            --updated_sam \
            {input.bam}  \
            {input.TEannotation}
        """



rule run_stellarscope_load:
    """
    Load the alignments intersecting TE annotation and remove PCR duplicates
    """
    output:
        directory(RESULTS_PATH_Stellarscope + "/results/stellarscope_out/{dataset}/{sample_id}/stload/"),
        RESULTS_PATH_Stellarscope + "/results/stellarscope_out/{dataset}/{sample_id}/stload/{sample_id}_stload-checkpoint.dedup_umi.pickle"

    input:
        #bam=get_stellarscope_bam,
        bam=RESULTS_PATH_STAR + "/results/STARoutdir/{dataset}/{sample_id}/score_Aligned.sortedByCB.bam",
        #barcodes=RESULTS_PATH_STAR + "/results/STARoutdir/{dataset}/{sample_id}/score_Solo.out/Gene/filtered/barcodes.tsv",
        barcodes=get_whitelist,
        TEannotation=get_TEannotation_dataset_unzipped
    log:
        "logs/stellarscope/stellarscope_stload_{dataset}_{sample_id}.txt"
    params:
        results_path=RESULTS_PATH_Stellarscope,
        strandedness=get_stranded_mode_stellarscope
    threads: 12
    resources:
        mem_gb=96
    benchmark:
        "benchmarks/run_stellarscope_load_{dataset}_{sample_id}.txt"
    conda: 
        "../envs/stellarscope.yml"
    shell:
        """
        # create a directory 
        mkdir -p "{params.results_path}/results/stellarscope_out"
        mkdir -p "{params.results_path}/results/stellarscope_out/{wildcards.dataset}"
        mkdir -p "{params.results_path}/results/stellarscope_out/{wildcards.dataset}/{wildcards.sample_id}/"
        mkdir -p "{params.results_path}/results/stellarscope_out/{wildcards.dataset}/{wildcards.sample_id}/stload/"
        
        outfolder="{params.results_path}/results/stellarscope_out/{wildcards.dataset}/{wildcards.sample_id}/stload/"

        # stellarscope analysis (load)
        stellarscope assign \
        --exp_tag {wildcards.sample_id}_stload \
        --outdir $outfolder \
        --nproc {threads} \
        {params.strandedness} \
        --whitelist <(gunzip -c {input.barcodes}) \
        --attribute transcript_id \
        --skip_em \
        --logfile logs/stellarscope/stellarscope_stload_{wildcards.dataset}_{wildcards.sample_id}.txt \
        --updated_sam \
        {input.bam} \
        {input.TEannotation}
        """


# then remove the tmp folders


rule run_stellarscope_pseudobulk:
    """
    Run stellarscope with pseudobulk pooling mode
    """
    output:
        directory(RESULTS_PATH_Stellarscope + "/results/stellarscope_out/{dataset}/{sample_id}/pseudobulk/")
    input:
        pickle=RESULTS_PATH_Stellarscope + "/results/stellarscope_out/{dataset}/{sample_id}/stload/{sample_id}_stload-checkpoint.dedup_umi.pickle"
    log:
        "logs/stellarscope/{dataset}_{sample_id}_resume_pseudobulk.txt"
    params:
        results_path=RESULTS_PATH_Stellarscope
    threads: 12
    resources:
        mem_gb=96
    benchmark:
        "benchmarks/run_stellarscope_pseudobulk_{dataset}_{sample_id}.txt"
    conda: 
        "../envs/stellarscope.yml"
    shell:
        """
        # create a directory for stellarscope pseudobulk results
        mkdir -p "{params.results_path}/results/stellarscope_out"
        mkdir -p "{params.results_path}/results/stellarscope_out/{wildcards.dataset}"
        mkdir -p "{params.results_path}/results/stellarscope_out/{wildcards.dataset}/{wildcards.sample_id}/"
        mkdir -p "{params.results_path}/results/stellarscope_out/{wildcards.dataset}/{wildcards.sample_id}/pseudobulk/"
        outfolder="{params.results_path}/results/stellarscope_out/{wildcards.dataset}/{wildcards.sample_id}/pseudobulk/"
        
        # use stellarscope resume to continue from the deduplication checkpoint
        stellarscope resume \
            --exp_tag {wildcards.sample_id}_pseudobulk \
            --outdir $outfolder \
            --nproc {threads} \
            --pooling_mode pseudobulk \
            --reassign_mode best_exclude \
            --logfile logs/stellarscope/{wildcards.dataset}_{wildcards.sample_id}_resume_pseudobulk.txt \
            --max_iter 500 \
            --updated_sam \
            {params.results_path}/results/stellarscope_out/{wildcards.dataset}/{wildcards.sample_id}/stload/{wildcards.sample_id}_stload-checkpoint.dedup_umi.pickle
        #  --use_every_reassign_mode 

        """




rule run_stellarscope_celltype:
    """
    Run stellarscope with celltype pooling mode
    """
    output:
        directory(RESULTS_PATH_Stellarscope + "/results/stellarscope_out/{dataset}/{sample_id}/celltype/")
    input:
        pickle=RESULTS_PATH_Stellarscope + "/results/stellarscope_out/{dataset}/{sample_id}/stload/{sample_id}_stload-checkpoint.dedup_umi.pickle",
        celltypes=DATA_PATH + "/data/{dataset}/annotation/{sample_id}/celltype_tsv.tsv"
    log:
        "logs/stellarscope/{dataset}_{sample_id}_run_celltype.log"
    params:
        results_path=RESULTS_PATH_Stellarscope
    threads: 12
    resources:
        mem_gb=96
    benchmark:
        "benchmarks/run_stellarscope_celltype_{dataset}_{sample_id}.txt"
    conda: 
        "../envs/stellarscope.yml"
    log:
        "logs/stellarscope/{dataset}_{sample_id}_resume_celltype.log"
    shell:
        """
        # create a directory for stellarscope pseudobulk results
        mkdir -p "{params.results_path}/results/stellarscope_out"
        mkdir -p "{params.results_path}/results/stellarscope_out/{wildcards.dataset}"
        mkdir -p "{params.results_path}/results/stellarscope_out/{wildcards.dataset}/{wildcards.sample_id}/"
        mkdir -p "{params.results_path}/results/stellarscope_out/{wildcards.dataset}/{wildcards.sample_id}/celltype/"
        outfolder="{params.results_path}/results/stellarscope_out/{wildcards.dataset}/{wildcards.sample_id}/celltype/"
        
        # use stellarscope resume to continue from the deduplication checkpoint
        stellarscope resume \
            --exp_tag {wildcards.sample_id}_celltype \
            --outdir $outfolder \
            --nproc {threads} \
            --pooling_mode celltype \
            --reassign_mode best_exclude \
            --logfile logs/stellarscope/stellarscope_resume_{wildcards.dataset}_{wildcards.sample_id}_celltype.txt \
            --max_iter 500 \
            --updated_sam \
            --celltype_tsv {input.celltypes} \
            {params.results_path}/results/stellarscope_out/{wildcards.dataset}/{wildcards.sample_id}/stload/{wildcards.sample_id}_stload-checkpoint.dedup_umi.pickle \
            2> {log}
        """