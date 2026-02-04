rule cellrangerTE_setup:
    """
    Run CellRangerTE script to generate annotation
    """
    output:
        directory(RESULTS_PATH_Cellranger + "/results/CellRangerTE/{genome_id}_CRTE_build/"),
        directory(RESULTS_PATH_Cellranger + "/results/CellRangerTE/{genome_id}_TE/")
    log:
        "../logs/CellRangerTE/CellRangerTE_{genome_id}.log"
    params:
        results_path=RESULTS_PATH_Cellranger,
        genome_cr_id=get_cr_genome_id,
        genome_cr_version=get_cr_genome_version
    threads: 10
    resources:
        mem_gb=200
    benchmark:
        "../benchmarks/CellRangerTE_{genome_id}.txt"
    singularity: "docker://litd/docker-cellranger:v9.0.0"
    shell:
        """
        sh code/CellRangerTE/CellRanger-TE_database_generation.sh -g {params.genome_cr_id} -r 35

        mv -r {params.genome_cr_id}_{params.genome_cr_version}_CRTE_build/ {params.results_path}/results/CellRangerTE/{wildcards.genome_id}_CRTE_build/"
        mv -r {params.genome_cr_id}_{params.genome_cr_version}_TE/ {params.results_path}/results/CellRangerTE/{wildcards.genome_id}_TE/"
        """
    

# cellranger count --transcriptome=GRCh38_GCv35_TE ...

