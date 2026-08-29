# SoloTE

rule setup_SoloTE:
    """
    Prepare TE annotation for SoloTE
    """
    output:
        "data/TEannotation/{genome_id}_rmsk.bed",
        "data/TEannotation/{genome_id}_rmsk.fa.out.gz"
    log:
        "logs/SoloTE/setup_{genome_id}.log"
    conda: 
        "../envs/SoloTE.yml"
    benchmark:
        "benchmarks/setup_SoloTE_{genome_id}.txt"
    shell:
        """
        python code/SoloTE/SoloTE_RepeatMasker_to_BED.py -g {wildcards.genome_id} 2> {log}

        # move the TE annotation files to the wanted folders
        mv {wildcards.genome_id}.fa.out.gz data/TEannotation
        mv {wildcards.genome_id}_rmsk.bed data/TEannotation
        # gzip data/TEannotation/{wildcards.genome_id}_rmsk.bed

        """


rule run_SoloTE:
    """
    Run SoloTE
    """
    output:
        directory(RESULTS_PATH_SoloTE + "/results/SoloTEout/{dataset}/{sample_id}/")
    input:
        bam=RESULTS_PATH_STAR+"/results/STARoutdir/{dataset}/{sample_id}/best_Aligned.sortedByCoord.out.bam",
        annotation=get_annotation_SoloTE
    log:
        "logs/SoloTE/{dataset}_{sample_id}.log"
    params:
        genome_id=get_genome,
        results_path=RESULTS_PATH_SoloTE
    threads: 12
    resources:
        mem_gb=64
    benchmark:
        "benchmarks/run_SoloTE_{dataset}_{sample_id}.txt"
    conda: 
        "../envs/SoloTE.yml"
    shell:
        """
        python code/SoloTE/SoloTE_pipeline.py \
            --threads {threads} \
        	--bam {input.bam} \
	        --teannotation {input.annotation} \
		    --outputprefix {wildcards.sample_id} \
		    --outputdir {params.results_path}/results/SoloTEout/{wildcards.dataset}/{wildcards.sample_id}/ 2> {log}
        mv {wildcards.sample_id}_SoloTE.stats {params.results_path}/results/SoloTEout/{wildcards.dataset}/{wildcards.sample_id}/
        """

rule run_SoloTE_MAPQthr:
    """
    Run SoloTE with the selected MAPQ threshold 
    """
    output:
        directory(RESULTS_PATH_SoloTE + "/results/SoloTEout_thr{MAPQthr}/{dataset}/{sample_id}/")
    input:
        bam=RESULTS_PATH_STAR+"/results/STARoutdir/{dataset}/{sample_id}/best_Aligned.sortedByCoord.out.bam",
        annotation=get_annotation_SoloTE
    log:
        "logs/SoloTE_thr{MAPQthr}/{dataset}_{sample_id}.log"
    params:
        genome_id=get_genome,
        results_path=RESULTS_PATH_SoloTE
    threads: 12
    resources:
        mem_gb=64
    benchmark:
        "benchmarks/run_SoloTE_thr{MAPQthr}_{dataset}_{sample_id}.txt"
    conda: 
        "../envs/SoloTE.yml"
    shell:
        """
        python code/SoloTE_selectMAPQthr/SoloTE_pipeline.py \
            --threads {threads} \
        	--bam {input.bam} \
	        --teannotation {input.annotation} \
		    --outputprefix {wildcards.sample_id} \
		    --outputdir {params.results_path}/results/SoloTEout_thr{wildcards.MAPQthr}/{wildcards.dataset}/{wildcards.sample_id}/ \
            --locusMAPQthr {wildcards.MAPQthr} 2> {log}
        mv {wildcards.sample_id}_SoloTE.stats {params.results_path}/results/SoloTEout_thr{wildcards.MAPQthr}/{wildcards.dataset}/{wildcards.sample_id}/
        """
