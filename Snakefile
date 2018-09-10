# Snakemake file for ChIP-Seq analysis

###############
# Libraries
###############
import os
import pandas as pd
from snakemake.utils import validate, min_version
#############################################
# Configuration and sample sheets
#############################################
configfile: "configs/config_tomato_sub.yaml"

#FASTQ_DIR = config["fastq_dir"]        # where to find the fastq files
WORKING_DIR = config["working_dir"]    # where you want to store your intermediate files (this directory will be cleaned up at the end)
RESULT_DIR = config["result_dir"]      # what you want to keep
SCRIPT_DIR = config["script_dir"]

GENOME_FASTA_URL = config["refs"]["genome_url"]
GENOME_FASTA_FILE = os.path.basename(config["refs"]["genome_url"])
TOTALCORES = 16                         #check this via 'grep -c processor /proc/cpuinfo'

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
SAMPLES = list(set(samples.index.values))

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
UNITS = units.index.get_level_values('unit').unique().tolist()

CASES = ['ChIP1']                       #complete list according to the sample.tsv file
CONTROLS = ['ChIP2']                    #complete list according to the sample.tsv file
###############
# Helper Functions
###############
def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

#def get_forward_fastq(wildcards):
    #return units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna()         #unused

#def get_reverse_fastq(wildcards):
    #return units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna()         #unused

def get_treatment(wildcards):
    return samples.loc[(wildcards.sample), ["condition"] == 'treatment'].dropna()

def get_control(wildcards):
    return samples.loc[(wildcards.sample), ["condition"] == 'control'].dropna()

##############
# Wildcards
##############
wildcard_constraints:
    sample = "[A-Za-z0-9]+"

wildcard_constraints:
    unit = "L[0-9]+"

##############
# Desired output
##############
FASTQC_REPORTS = expand(RESULT_DIR + "fastqc/{sample}_{unit}_{pair}_fastqc.zip", sample=SAMPLES,unit=UNITS, pair={"forward", "reverse"})
BAM_INDEX = expand(RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam.bai", sample=SAMPLES,unit=UNITS)
BAM_RMDUP = expand(RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam", sample=SAMPLES,unit=UNITS)
BEDGRAPH = expand(RESULT_DIR + "bedgraph/{sample}_{unit}.sorted.rmdup.bedgraph", sample=SAMPLES,unit=UNITS)
BIGWIG = expand(RESULT_DIR + "bigwig/{sample}_{unit}.bw", sample=SAMPLES,unit=UNITS)
BAM_COMPARE = expand(RESULT_DIR + "bamcompare/log2_{treatment}_{control}_{unit}.bamcompare.bw", zip, treatment = CASES, control = CONTROLS,unit=UNITS) #add zip function in the expand to compare respective treatment and control
BED_NARROW = expand(RESULT_DIR + "bed/{treatment}_vs_{control}_{unit}_peaks.narrowPeak", zip, treatment = CASES, control = CONTROLS,unit=UNITS)
BED_BROAD = expand(RESULT_DIR + "bed/{treatment}_vs_{control}_{unit}_peaks.broadPeak", zip, treatment = CASES, control = CONTROLS,unit=UNITS)
################
# Final output
################
rule all:
    input:
        BAM_INDEX,
        BAM_RMDUP,
        FASTQC_REPORTS,
        BEDGRAPH,
        BIGWIG,
        BAM_COMPARE,
        BED_NARROW,
        BED_BROAD
    message: "ChIP-seq pipeline succesfully run."		#finger crossed to see this message!

    shell:"#rm -rf {WORKING_DIR}"

###############
# Rules
###############
rule get_genome_fasta:
    output:
        WORKING_DIR + "genome.fasta"
    message:"downloading {GENOME_FASTA_FILE} genomic fasta file"
    shell: "wget -O {output} {GENOME_FASTA_URL}"

rule trimmomatic:
    input:
        reads = get_fastq,
        adapters = config["adapters"]
    output:
        forward_reads = WORKING_DIR + "trimmed/{sample}_{unit}_forward.fastq.gz",
        reverse_reads = WORKING_DIR + "trimmed/{sample}_{unit}_reverse.fastq.gz",
        forwardUnpaired = temp(WORKING_DIR + "trimmed/{sample}_{unit}_forward_unpaired.fastq.gz"),
        reverseUnpaired = temp(WORKING_DIR + "trimmed/{sample}_{unit}_reverse_unpaired.fastq.gz")
    message: "trimming {wildcards.sample} reads"
    log:
        RESULT_DIR + "logs/trimmomatic/{sample}_{unit}.log"
    params :
        seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold =      str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['trimmomatic']['minReadLength']),
        phred = 		            str(config["trimmomatic"]["phred"])
    threads: 10
    shell:
        "trimmomatic PE {params.phred} -threads {threads} "
        "{input.reads} "
        "{output.forward_reads} "
        "{output.forwardUnpaired} "
        "{output.reverse_reads} "
        "{output.reverseUnpaired} "
        "ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} 2>{log}"

rule fastqc:
    input:
        fwd=WORKING_DIR + "trimmed/{sample}_{unit}_forward.fastq.gz",
        rev=WORKING_DIR + "trimmed/{sample}_{unit}_reverse.fastq.gz"
    output:
        fwd=RESULT_DIR + "fastqc/{sample}_{unit}_forward_fastqc.zip",
        rev=RESULT_DIR + "fastqc/{sample}_{unit}_reverse_fastqc.zip"
    log:
        RESULT_DIR + "logs/fastqc/{sample}_{unit}.fastqc.log"
    params:
        RESULT_DIR + "fastqc/"
    message:
        "---Quality check of trimmed {wildcards.sample} sample with FASTQC" 		#removed, it was not working
    shell:
        "fastqc --outdir={params} {input.fwd} {input.rev} 2>{log}"


rule index:
    input:
        WORKING_DIR + "genome.fasta"
    output:
        [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)],
        WORKING_DIR + "genome.rev.1.bt2",
        WORKING_DIR + "genome.rev.2.bt2"
    message:"indexing genome"
    params:
        WORKING_DIR + "genome"
    threads: 10
    shell:"bowtie2-build --threads {threads} {input} {params}"

rule align:
    input:
        forward = WORKING_DIR + "trimmed/{sample}_{unit}_forward.fastq.gz",
        reverse = WORKING_DIR + "trimmed/{sample}_{unit}_reverse.fastq.gz",
        forwardUnpaired = WORKING_DIR + "trimmed/{sample}_{unit}_forward_unpaired.fastq.gz",
        reverseUnpaired = WORKING_DIR + "trimmed/{sample}_{unit}_reverse_unpaired.fastq.gz",
        index = [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)]
    output:
        temp(WORKING_DIR + "mapped/{sample}_{unit}.sam")
    message: "Mapping files"
    params:
        bowtie = " ".join(config["bowtie2"]["params"].values()), #take argument separated as a list separated with a space
        index = WORKING_DIR + "genome"
    threads: 10
    shell:
        "bowtie2 {params.bowtie} "
        "--threads {threads} "
        "-x {params.index} "
        "-1 {input.forward} -2 {input.reverse} "
        "-U {input.forwardUnpaired},{input.reverseUnpaired} "   # also takes the reads unpaired due to trimming
        "-S {output}"                                           #changed to match the pipeline of cambridge

rule header_sam:
    input:
        WORKING_DIR + "mapped/{sample}_{unit}.sam"
    output:
        WORKING_DIR + "mapped/{sample}_{unit}_header.sam"
    message : "Creation of the header for {wildcards.sample}"
    shell:
        "sed -n 1,15p {input} > {output}"

rule unique_align:
    input:
        WORKING_DIR + "mapped/{sample}_{unit}.sam"
    output:
        WORKING_DIR + "mapped/{sample}_{unit}_unique.txt"
    shell:
        "samtools view -S -f0x02 {input} "
        "| grep -v "XS:i:" > {output} "                         # grep -v : look for "XS:i:" which represent multi align read and -v exclude them from the text file


rule multi_align:
    input:
        WORKING_DIR + "mapped/{sample}_{unit}.sam"
    output:
        WORKING_DIR + "mapped/{sample}_{unit}_multi.txt"
    shell:
        "samtools view -S -f0x02 {input} "
        "| grep  "XS:i:" > {output} "

rule python:
    input:
        unique = WORKING_DIR + "mapped/{sample}_{unit}_unique.txt"
        multi  = WORKING_DIR + "mapped/{sample}_{unit}_multi.txt"
    output:
        unique = WORKING_DIR + "mapped/{sample}_{unit}_unique.py.txt"
        multi  = WORKING_DIR + "mapped/{sample}_{unit}_multi.py.txt"
    script:
        "scripts/foo.py"
rule word_count:
    input:
        WORKING_DIR + "mapped/{sample}_{unit}_multi.py.txt"
    output:
        WORKING_DIR + "mapped/{sample}_{unit}_multi.py.stats"
    shell:
        "wc -l {input} >> {output}"

rule multi_to_unique:
    input:
        txt   = WORKING_DIR + "mapped/{sample}_{unit}_multi.py.txt"
        stats = WORKING_DIR + "mapped/{sample}_{unit}_multi.py.stats"
    output:
        txt   = WORKING_DIR + "mapped/{sample}_{unit}_multi.to.unique.txt"
        RDATA = WORKING_DIR + "mapped/{sample}_{unit}_MU.RData"
    script:
        "scripts/multi_unique_extract.r {input.txt} {input.stats} {output.RDATA} {output.txt}  15 20000 21 16"


rule txt_to_sam:
    input:
    output:
    shell:

rule sam_to_bam:
    input:
    output:
    shell:

rule merge_bam:
    input:
    output:
    shell:

rule sort:
    input:
        WORKING_DIR + "mapped/{sample}_{unit}.bam"
    output:
        RESULT_DIR + "mapped/{sample}_{unit}.sorted.bam"
    message:"sorting {wildcards.sample} bam file"
    threads: 10
    shell:"samtools sort -@ {threads} -o {output} {input}"

rule rmdup:
    input:
        RESULT_DIR + "mapped/{sample}_{unit}.sorted.bam"
    output:
        RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam"
    message: "Removing duplicate from file {input}"
    shell:
        "samtools rmdup {input} {output}"                       #samtools manual says "This command is obsolete. Use markdup instead."

rule bam_index:
    input:
        RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam.bai"
    message: "Indexing {wildcards.sample} for rapid access"
    shell:
        "samtools index {input}"

rule bedgraph:
    input:
        RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bedgraph/{sample}_{unit}.sorted.rmdup.bedgraph"
    params:
        genome = WORKING_DIR + "genome"
    message:
        "Creation of {input} bedgraph file"
    shell:
        "bedtools genomecov -bg -ibam {input} -g {params.genome} > {output}"
        # require a sorted bam file as input
        # -ibam the input file is in BAM format
        # -bga  Report Depth in BedGraph format, regions with zero coverage are also reported. Extract those regions with "grep -w 0$"
        # -pc Calculate coverage of pair-end fragments. Works for BAM files only.

rule bigwig:
    input:
        RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bigwig/{sample}_{unit}.bw"
    message:
        "Converting {input} bam into bigwig file"
    log:
        RESULT_DIR + "logs/deeptools/{sample}_{unit}_bamtobigwig.log"
    params:
        EFFECTIVEGENOMESIZE = str(config["bamCoverage"]["params"]["EFFECTIVEGENOMESIZE"]) #take argument separated as a list separated with a space
    shell:
        "bamCoverage --bam {input} -o {output} --effectiveGenomeSize {params.EFFECTIVEGENOMESIZE} 2>{log}"

rule bamcompare:
    input:
        treatment = RESULT_DIR + "mapped/{treatment}_{unit}.sorted.rmdup.bam",              #input requires an indexed bam file
        control = RESULT_DIR + "mapped/{control}_{unit}.sorted.rmdup.bam"                   #input requires an indexed bam file
    output:
        bigwig = RESULT_DIR + "bamcompare/log2_{treatment}_{control}_{unit}.bamcompare.bw"
    message:
        "Running bamCompare"
    shell:
        "bamCompare -b1 {input.treatment} -b2 {input.control} -o {output.bigwig}"

rule call_narrow_peaks:
    input:
        treatment = RESULT_DIR + "mapped/{treatment}_{unit}.sorted.rmdup.bam",
        control = RESULT_DIR + "mapped/{control}_{unit}.sorted.rmdup.bam"
    output:
        bed = RESULT_DIR + "bed/{treatment}_vs_{control}_{unit}_peaks.narrowPeak"
    message:
        "Calling narrowPeak"
    params:
        name = "{treatment}_vs_{control}_{unit}",        #this option will give the output name, has to be similar to the output
        format = str(config['macs2']['format']),
        genomesize = str(config['macs2']['genomesize']),
        qvalue = str(config['macs2']['qvalue'])
    conda:
        "envs/mac2_env.yaml"
    shell:
        """
        source activate macs2
        macs2 callpeak -t {input.treatment} -c {input.control} {params.format} {params.genomesize} --name {params.name} --nomodel --bdg -q {params.qvalue} --outdir results/bed/
        """
# -g define the mappable genome size, for human change 'mm' to 'hs'
# --name will be used to create output files like NAME_peaks.xls', 'NAME_negative_peaks.xls', 'NAME_peaks.bed' , 'NAME_summits.bed', 'NAME_model.r'
# --bdg provides the files for the calculation of the FDR
# -q define the minimum FDR to call significant region, default is 0.05
rule call_broad_peaks:
    input:
        treatment = RESULT_DIR + "mapped/{treatment}_{unit}.sorted.rmdup.bam",
        control = RESULT_DIR + "mapped/{control}_{unit}.sorted.rmdup.bam"
    output:
        bed = RESULT_DIR + "bed/{treatment}_vs_{control}_{unit}_peaks.broadPeak"
    message:
        "Calling broadPeak"
    params:
        name = "{treatment}_vs_{control}_{unit}",
        format = str(config['macs2']['format']),
        genomesize = str(config['macs2']['genomesize']),
        qvalue = str(config['macs2']['qvalue'])
    conda:
        "envs/mac2_env.yaml"
    shell:
        """
        source activate macs2
        macs2 callpeak -t {input.treatment} -c {input.control} {params.format} --broad --broad-cutoff 0.1 {params.genomesize} --name {params.name} --nomodel --bdg -q {params.qvalue} --outdir results/bed/
        """
