configfile: "config/config.yml"

# For indexing

def get_genome(wildcards):
    return config["datasets"][wildcards.dataset]["genome"]

def get_genome_fasta(wildcards):
    return config["genomes"][wildcards.genome]["fasta"]

def get_read_length(wildcards):
    return config["datasets"][wildcards.dataset]["read_length"]

def get_fasta(wildcards):
    return DATA_PATH + "/" + config["genomes"][wildcards.genome_id]["fasta"]


def get_annotation(wildcards):
    return DATA_PATH + "/" + config["genomes"][wildcards.genome_id]["gtf"]

# get TE annotation rmsk gtf, gz, from genome wildcard
def get_TEannotation_genome(wildcards):
    return DATA_PATH + "/" + config["genomes"][wildcards.genome_id]["rmsk_gtf"]
# get TE annotation rmsk gtf, gz, from dataset wildcard
def get_TEannotation_dataset(wildcards):
    genomeid=config["datasets"][wildcards.dataset]["genome"]
    return DATA_PATH + "/" + config["genomes"][genomeid]["rmsk_gtf"]
# get TE annotation rmsk gtf, unzipped, from genome wildcard
def get_TEannotation_dataset_unzipped(wildcards):
    genomeid=config["datasets"][wildcards.dataset]["genome"]
    return f"tmp/TEannotation_{genomeid}.gtf"
# get TE annotation where the gene_id is the locus instead of the family, unzipped, from dataset wildcard
def get_TEannotation_locus_dataset_unzipped(wildcards):
    genomeid=config["datasets"][wildcards.dataset]["genome"]
    return f"tmp/TEannotation_locus_{genomeid}.gtf"



# For alignment

def get_annotation_dataset(wildcards):
    genomeid=config["datasets"][wildcards.dataset]["genome"]
    return f"tmp/genome_annotation_{genomeid}"

def get_whitelist(wildcards):
    return config["datasets"][wildcards.dataset]["whitelist"]

def get_nmax(wildcards):
    return config["conf"][wildcards.conf]["nmax"]
def get_UMIlength(wildcards):
    return config["datasets"][wildcards.dataset]["UMI_length"]
def get_extra_params(wildcards):
    return config["conf"][wildcards.conf]["extra_params"]
    
def get_index_dir(wildcards):
    genomeid=config["datasets"][wildcards.dataset]["genome"]
    readlength=config["datasets"][wildcards.dataset]["read_length"]
    return directory(f"{RESULTS_PATH_STAR}/indexes/genomeIndexes_{genomeid}_{readlength}")

def get_index_TE_dir(wildcards):
    genomeid=config["datasets"][wildcards.dataset]["genome"]
    readlength=config["datasets"][wildcards.dataset]["read_length"]
    return directory(f"{RESULTS_PATH_STAR}/indexes/TE_genomeIndexes_{genomeid}_{readlength}")




# for STARsolo with TEs
def get_GeneTEannotation_gtf(wildcards):
    genomeid=config["datasets"][wildcards.dataset]["genome"]
    return DATA_PATH + "/" + config["genomes"][genomeid]["GeneTE_gtf"]


# For SoloTE
# I could also use TEannotation/{wildcards.genome_id}_rmsk.bed
def get_annotation_SoloTE(wildcards):
    genomeid=config["datasets"][wildcards.dataset]["genome"]
    return f"data/TEannotation/{genomeid}_rmsk.bed"


# For CellrangerTE
def get_cr_genome_id(wildcards):
    return config["genomes"][wildcards.genome_id]["cr_id"]
def get_cr_genome_version(wildcards):
    return config["genomes"][wildcards.genome_id]["cr_version"]


# get all sample names for a given dataset (needed for MATES)
def get_sample_names(wildcards):
    return config["datasets"][wildcards.dataset]["sample_ids"]

def get_bam_indexes(wildcards):
    dataset=wildcards.dataset
    sample_ids=list(config["datasets"][wildcards.dataset]["sample_ids"])
    return expand(RESULTS_PATH_STAR + f"/results/STARoutdir/{dataset}/{{sample_ids}}/score_Aligned.sortedByCoord.out.bam.bai", sample_ids=sample_ids)

# translate genome id into species name (needed for MATES)
def get_species(wildcards):
    if wildcards.genome_id=="mm10":
        return "Mouse"
    elif wildcards.genome_id=="hg38":
        return "Human"

# get path to bed file of TE annotation for MATES
def get_TE_nooverlap_bed(wildcards):
    genomeid=config["datasets"][wildcards.dataset]["genome"]
    return RESULTS_PATH_MATES + "/{genomeid}/TE_nooverlap.bed"

def get_strandedness(wildcards):
    strandedness = config["datasets"][wildcards.dataset]["strandedness"]
    return strandedness

def get_stranded_mode_stellarscope(wildcards):
    strandedness = config["datasets"][wildcards.dataset]["strandedness"]
    
    if strandedness == 'Unstranded':
        return ""
    elif strandedness == "Forward":
        return "--stranded_mode F"
    elif strandedness == "Reverse":
        return "--stranded_mode R"
         

def get_stellarscope_bam(wildcards):
    if config["datasets"][wildcards.dataset]["simulated"]==True:
        bam=RESULTS_PATH_STAR + f"/results/STARoutdir/{wildcards.dataset}/{wildcards.sample_id}/score_Aligned.changedStrand.bam",
    else:
        bam=RESULTS_PATH_STAR + f"/results/STARoutdir/{wildcards.dataset}/{wildcards.sample_id}/score_Aligned.sortedByCB.bam"
    return bam

def get_strandedness_stellarscope(wildcards):
    if config["datasets"][wildcards.dataset]["strandedness"]=="Reverse":
        strandedness="R"
    else:
        strandedness="F" # forward stranded or unstranded, which will be first made forward stranded
    return strandedness