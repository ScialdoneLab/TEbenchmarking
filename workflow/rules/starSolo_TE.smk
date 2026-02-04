rule create_TE_STARsolo_annotation:
    """
    Create TE gtf file for STARsolo (put locus ID in the gene_id attribute)
    """
    output:
        "data/TEannotation/{genome_id}_GENCODE_rmsk_TElocus.gtf.gz", 
        "tmp/TEannotation_locus_{genome_id}.gtf" # add temp()
    input:
        fasta=get_fasta, 
        annotation=get_TEannotation_genome
    log:
        "logs/unzip_TEannotation/{genome_id}.log"
    threads: 6
    resources:
        mem_gb=8
    benchmark:
        "../benchmarks/unzip_TEannotation_{genome_id}.txt"
    conda:
        "../envs/align.yml"
    shell:
        """
        zcat {input.annotation} | sed 's/gene_name "[^"]*";\\s*//g' |
        awk '{{
            match($0, /gene_id "([^"]+)"/, g);
            match($0, /transcript_id "([^"]+)"/, t);
            if (t[1] != "") {{
                sub(/gene_id "[^"]+"/, "gene_id \\"" t[1] "\\"");
            }}
            print
        }}' |  gzip > data/TEannotation/{wildcards.genome_id}_GENCODE_rmsk_TElocus.gtf.gz
        
        gunzip -c data/TEannotation/{wildcards.genome_id}_GENCODE_rmsk_TElocus.gtf.gz > \
        tmp/TEannotation_locus_{wildcards.genome_id}.gtf
        """


# STARsolo EM for TE quantification

rule index_genome_TE:
    """
    Index a genome using star.
    """
    output:
        directory(RESULTS_PATH_STAR+"/indexes/TE_genomeIndexes_{genome_id}_{read_length}")
    input:
        fasta="tmp/genome_fasta_{genome_id}", 
        annotation="tmp/TEannotation_locus_{genome_id}.gtf"
    log:
        "logs/index_genome_TE/{genome_id}_{read_length}.log"
    threads: 12
    resources:
        mem_gb=48
    benchmark:
        "../benchmarks/index_genome_TE_{genome_id}_{read_length}.txt"
    conda:
        "../envs/align.yml"
    shell:
        """
        # compute overhang parameter based on read length
        overhang=$(({wildcards.read_length}-1))

        # index the genome with star
        STAR --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --runThreadN {threads} \
            --sjdbOverhang $overhang \
            --sjdbGTFfile {input.annotation} \
            > {log}
        """


rule run_STARsolo_EM:
    """
    STARsolo to quantify TE expression with EM algorithm
    """
    output:
        directory(RESULTS_PATH_STARsoloTE + "/results/STARsolo_EM_TE/{dataset}/{sample_id}/TE_Solo.out")
    input:
        fastaFolder=DATA_PATH+"/data/{dataset}/fastqs/{sample_id}/",
        refDir=get_index_TE_dir,
        TEannotation=get_TEannotation_locus_dataset_unzipped
    log:
        "logs/runSTARsolo_EM_TE/{dataset}_{sample_id}.log"
    params:
        whitelist=get_whitelist,
        UMIlen=get_UMIlength,
        results_path=RESULTS_PATH_STARsoloTE,
        strand=get_strandedness
    threads: 12
    resources:
        mem_gb=96
    benchmark:
        "benchmarks/STARsolo_EM_TE_align_{dataset}_{sample_id}.txt"
    conda:
        "../envs/align.yml"
    shell:
        """
        # create comma separated lists of R1s and R2s
        R1s=`ls {input.fastaFolder} | grep R1 | sed "s|^|{input.fastaFolder}/|" | paste -sd,`
        R2s=`ls {input.fastaFolder} | grep R2 | sed "s|^|{input.fastaFolder}/|" | paste -sd,`
        echo $R1s 
        echo $R2s

        # run genome alignment with STAR
        STAR --genomeDir {input.refDir} \
        	--readFilesIn $R2s $R1s \
            --readFilesCommand gunzip -c \
	        --runThreadN {threads} \
            --outFilterMismatchNmax 18 \
            --winAnchorMultimapNmax 100 \
            --outFilterMultimapNmax 100 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --alignMatesGapMax 1000000 \
            --outFilterScoreMin 30 \
            --outFileNamePrefix {params.results_path}/results/STARsolo_EM_TE/{wildcards.dataset}/{wildcards.sample_id}/TE_ \
            --outSAMtype BAM SortedByCoordinate \
	        --soloType CB_UMI_Simple \
	        --soloCBwhitelist  <(gunzip -c {params.whitelist}) \
	        --soloUMIlen {params.UMIlen} \
	        --outSAMattributes NH HI AS nM NM MD jM jI XS MC ch cN CR CY UR UY GX GN CB UB sM sS sQ \
            --sjdbGTFfile {input.TEannotation} \
            --soloMultiMappers EM \
            --soloFeatures Gene \
            --soloStrand {params.strand} \
            > {log}
        """



