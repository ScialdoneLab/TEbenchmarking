# STAR Alignment

## Indexing

rule unzip_genomes:
    """
    Unzip genome and annotation
    """
    output:
        "tmp/genome_fasta_{genome_id}", #add temp() to these
        "tmp/genome_annotation_{genome_id}"
    input:
        fasta=get_fasta, 
        annotation=get_annotation
    log:
        "logs/unzip_genome/{genome_id}.log"
    threads: 6
    resources:
        mem_gb=8
    benchmark:
        "../benchmarks/unzip_genomes_{genome_id}.txt"
    conda:
        "../envs/align.yml"
    shell:
        """
        # unzip fasta file
        gunzip -c {input.fasta} > tmp/genome_fasta_{wildcards.genome_id} 2> {log}
        gunzip -c {input.annotation} > tmp/genome_annotation_{wildcards.genome_id}
        """


rule index_genome:
    """
    Index a genome using star.
    """
    output:
        directory(RESULTS_PATH_STAR+"/indexes/genomeIndexes_{genome_id}_{read_length}")
    input:
        fasta="tmp/genome_fasta_{genome_id}", 
        annotation="tmp/genome_annotation_{genome_id}"
    log:
        out = 'logs/index_genome/{genome_id}_{read_length}_stdout.log',
        err = 'logs/index_genome/{genome_id}_{read_length}_stderr.err'
    threads: 12
    resources:
        mem_gb=48
    benchmark:
        "../benchmarks/index_genome_{genome_id}_{read_length}.txt"
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
            --sjdbGTFfile {input.annotation} 2> {log.err} 1> {log.out}
        """


## Alignment



rule STARsolo_align:
    """
    STARsolo with parameters for TE analysis
    """
    output:
        RESULTS_PATH_STAR + "/results/STARoutdir/{dataset}/{sample_id}/{conf}_Aligned.sortedByCoord.out.bam",
        RESULTS_PATH_STAR + "/results/STARoutdir/{dataset}/{sample_id}/{conf}_Solo.out/Gene/filtered/barcodes.tsv",
        RESULTS_PATH_STAR + "/results/STARoutdir/{dataset}/{sample_id}/{conf}_Solo.out/Gene/raw/barcodes.tsv"
    input:
        fastaFolder=DATA_PATH+"/data/{dataset}/fastqs/{sample_id}/",
        refDir=get_index_dir,
        annotation=get_annotation_dataset
    log:
        "logs/align_to_genome/{dataset}_{sample_id}_{conf}.log"
    params:
        whitelist=get_whitelist,
        UMIlen=get_UMIlength,
        extra_params=get_extra_params,
        results_path=RESULTS_PATH_STAR,
        strand=get_strandedness
    threads: 12
    resources:
        mem_gb=96
    benchmark:
        "../benchmarks/STARsolo_align_{dataset}_{sample_id}_{conf}.txt"
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
            --outFilterMultimapScoreRange 5 \
            --outFileNamePrefix {params.results_path}/results/STARoutdir/{wildcards.dataset}/{wildcards.sample_id}/{wildcards.conf}_ \
            --outSAMtype BAM SortedByCoordinate \
	        --soloType CB_UMI_Simple \
	        --soloCBwhitelist  <(gunzip -c {params.whitelist}) \
	        --soloUMIlen {params.UMIlen} \
	        --outSAMattributes NH HI AS nM NM MD jM jI XS MC ch cN CR CY UR UY GX GN CB UB sM sS sQ \
            --sjdbGTFfile {input.annotation} \
            --clipAdapterType CellRanger4 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
            --limitOutSJcollapsed 5000000 \
            --soloStrand {params.strand} \
            {params.extra_params} \
            &> {log}
        """
 