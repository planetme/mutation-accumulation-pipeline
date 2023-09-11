import pathlib as pathlib
import glob as glob

configfile: "config.yaml"

#use config.yaml to get the outputdirs instead for portability
REF_FNA_FILE = config['ref_fna_file']
REF_FNA_DICT = config['ref_fna_dict']
REF_FNA_FAI = config['ref_fna_fai']
REF_FNA_INDEXING_LOG = config['ref_fna_indexing_log']
READS_DIR = config['reads_dir']
FASTQ_SCREEN_REFERENCES = config['fastq_screen_references']
ADAPTERS = config['adapters']
READ_GROUP_SCRIPT = config['read_group_script']

#FASTQC_REPORTS_DIR = config['fastqc_reports_dir']
PREQC_DIR = config['preqc_dir']
ADAPTRIM_DIR = config['adaptrim_dir']
ADAPTRIM_QC_DIR = config['adaptrim_qc_dir']
ADAPTRIM_QUALTRIM_DIR = config['adaptrim_qualtrim_dir']
ADAPTRIM_QUALTRIM_QC_DIR = config['adaptrim_qualtrim_qc_dir']
SCREENING_DIR = config['screening_dir']
FASTQ_SCREEN_CONF = config['fastq_screen_conf']
FASTQ_SCREEN_LOGFILE_EXTENSION = config['fastq_screen_logfile_extension']
SNAKEMAKE_LOGS = config['snakemake_logs']

ALIGNED_BAM_DIR = config['aligned_bam_dir']
SAMTOOLS_REPORTS_DIR = config['samtools_reports_dir']
VCFS_DIR = config['vcfs_dir']
GVCFS_DIR = config['gvcfs_dir']
PICARD_DIR = config['picard_dir']

READ1_TAG = config['read1_tag']
READ2_TAG = config['read2_tag']
ADAPTRIM1_TAG = config['adaptrim1_tag']
ADAPTRIM2_TAG = config['adaptrim2_tag']
QUALTRIM_P_TAG = config['qualtrim_p_tag']
QUALTRIM_U_TAG = config['qualtrim_u_tag']
#ADAPTRIM_QUALTRIM_1P_TAG = config['adaptrim_qualtrim_1p_tag']
#ADAPTRIM_QUALTRIM_2P_TAG = config['adaptrim_qualtrim_2p_tag']

READS_FILE_EXTENSION = config['reads_file_extension']
PREQC_FILE_EXTENSION = config['preqc_file_extension']
ADAPTRIM_FILE_EXTENSION = config['adaptrim_file_extension']
ADAPTRIM_QC_FILE_EXTENSION = config['adaptrim_qc_file_extension']
ADAPTRIM_QUALTRIM_FILE_EXTENSION = config['adaptrim_qualtrim_file_extension']
ADAPTRIM_QUALTRIM_QC_FILE_EXTENSION = config['adaptrim_qualtrim_qc_file_extension']
ADAPTRIM_QUALTRIM_QC_HTML_FILE_EXTENSION = config['adaptrim_qualtrim_qc_html_file_extension']
ALIGNED_BAM_FILE_EXTENSION = config['aligned_bam_file_extension']
ALIGNED_BAM_FLAGSTAT_FILE_EXTENSION = config['aligned_bam_flagstat_file_extension']
INDEXED_ALIGNED_BAM_FILE_EXTENSION = config['indexed_aligned_bam_file_extension']
GVCFS_FILE_EXTENSION = config['gvcfs_file_extension']
VCFS_FILE_EXTENSION = config['vcfs_file_extension']

# Single file patterns: use Python string formatting to build
# The safest way to mix global variables and wildcards in a formatted string is to remember the following:

# Global variables are surrounded in single curly braces (e.g. {INPUT_DIR}).
# Wildcards are surrounded with double curly braces (e.g. {{reads_file_prefix}}).
# Use upper-case for globals and lower-case for wildcards.
#
# Now define all the wildcard patterns that either depend on
# directory and file configuration, or are used more than once.

# Note the use of single curly braces for global variables
# and double curly braces for snakemake wildcards
# Calbicans-1_S29_L006_R1_001.fastq.gz
# reads_file_prefix = Calbicans-1_S29_L006_
# reads_file_suffix = _001
SAMPLE = f'{{reads_file_prefix}}'

READ_ONE_FILE = f'{READS_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{READS_FILE_EXTENSION}'
READ_TWO_FILE = f'{READS_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{READS_FILE_EXTENSION}'
#TEST_FILE=f'{BAM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{SAM_FILE_EXTENSION}'

# Calbicans-1_S29_L006_R1_001_fastqc.zip
PREQC_ONE_FILE = f'{PREQC_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{PREQC_FILE_EXTENSION}'
PREQC_TWO_FILE = f'{PREQC_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{PREQC_FILE_EXTENSION}'

# Calbicans-1_S29_L006_R1_001.adaptrim.1P.fq.gz
ADAPTRIM_ONE_FILE = f'{ADAPTRIM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{ADAPTRIM1_TAG}{ADAPTRIM_FILE_EXTENSION}'
ADAPTRIM_TWO_FILE = f'{ADAPTRIM_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{ADAPTRIM2_TAG}{ADAPTRIM_FILE_EXTENSION}'

# Calbicans-1_S29_L006_R1_001.adaptrim.1P_fastqc.zip
ADAPTRIM_QC_ONE_FILE = f'{ADAPTRIM_QC_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{ADAPTRIM1_TAG}{ADAPTRIM_QC_FILE_EXTENSION}'
ADAPTRIM_QC_TWO_FILE = f'{ADAPTRIM_QC_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{ADAPTRIM2_TAG}{ADAPTRIM_QC_FILE_EXTENSION}'

# Calbicans-1_S29_L006_R1_001_P.adaptrim.qualtrim.fq.gz
ADAPTRIM_QUALTRIM_ONE_P_FILE = f'{ADAPTRIM_QUALTRIM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{QUALTRIM_P_TAG}{ADAPTRIM_QUALTRIM_FILE_EXTENSION}'
ADAPTRIM_QUALTRIM_ONE_U_FILE = f'{ADAPTRIM_QUALTRIM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{QUALTRIM_U_TAG}{ADAPTRIM_QUALTRIM_FILE_EXTENSION}'
ADAPTRIM_QUALTRIM_TWO_P_FILE = f'{ADAPTRIM_QUALTRIM_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{QUALTRIM_P_TAG}{ADAPTRIM_QUALTRIM_FILE_EXTENSION}'
ADAPTRIM_QUALTRIM_TWO_U_FILE = f'{ADAPTRIM_QUALTRIM_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{QUALTRIM_U_TAG}{ADAPTRIM_QUALTRIM_FILE_EXTENSION}'
#(config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim.fq.gz")
ADAPTRIM_QUALTRIM_BASEOUT_FILE = f'{ADAPTRIM_QUALTRIM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{ADAPTRIM_QUALTRIM_FILE_EXTENSION}'

# Calbicans-1_S29_L006_R1_001.adaptrim.qualtrim_1P_fastqc.zip
#fastqc does not support output file name
#ADAPTRIM_QUALTRIM_QC_ONE_FILE = f'{ADAPTRIM_QUALTRIM_QC_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{ADAPTRIM_QUALTRIM_1P_TAG}{ADAPTRIM_QUALTRIM_QC_FILE_EXTENSION}'
#ADAPTRIM_QUALTRIM_QC_TWO_FILE = f'{ADAPTRIM_QUALTRIM_QC_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{ADAPTRIM_QUALTRIM_2P_TAG}{ADAPTRIM_QUALTRIM_QC_FILE_EXTENSION}'

# Calbicans-1_S29_L006_R1_001_P.adaptrim.qualtrim_fastqc.zip
ADAPTRIM_QUALTRIM_QC_ONE_FILE = f'{ADAPTRIM_QUALTRIM_QC_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{QUALTRIM_P_TAG}{ADAPTRIM_QUALTRIM_QC_FILE_EXTENSION}'
ADAPTRIM_QUALTRIM_QC_TWO_FILE = f'{ADAPTRIM_QUALTRIM_QC_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{QUALTRIM_P_TAG}{ADAPTRIM_QUALTRIM_QC_FILE_EXTENSION}'

# Calbicans-1_S29_L006_R1_001.adaptrim.qualtrim.html
# ADAPTRIM_QUALTRIM_QC_ONE_SCREEN_FILE = f'{ADAPTRIM_QUALTRIM_QC_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{QUALTRIM_P_TAG}{ADAPTRIM_QUALTRIM_QC_HTML_FILE_EXTENSION}'
# ADAPTRIM_QUALTRIM_QC_TWO_SCREEN_FILE = f'{ADAPTRIM_QUALTRIM_QC_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{QUALTRIM_P_TAG}{ADAPTRIM_QUALTRIM_QC_HTML_FILE_EXTENSION}'

FASTSCREEN_LOGFILE = f'{SCREENING_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{FASTQ_SCREEN_LOGFILE_EXTENSION}'

# Calbicans-1_S29_L006_001_sorted.bam
#ALIGNED_BAM_ONE_FILE = f'{BAM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{ALIGNED_BAM_FILE_EXTENSION}'
#ALIGNED_BAM_TWO_FILE = f'{BAM_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{ALIGNED_BAM_FILE_EXTENSION}'
ALIGNED_BAM_FILE = f'{ALIGNED_BAM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{ALIGNED_BAM_FILE_EXTENSION}'

# Calbicans-1_S29_L006_001_sorted.flagstat.txt
ALIGNED_BAM_FLAGSTAT_FILE = f'{SAMTOOLS_REPORTS_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{ALIGNED_BAM_FLAGSTAT_FILE_EXTENSION}'

#MERGED_BAM_ONE_FILE = f'{BAM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{MERGED_BAM_FILE_EXTENSION}'
#MERGED_BAM_TWO_FILE = f'{BAM_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{MERGED_BAM_FILE_EXTENSION}'

# Calbicans-1_S29_L006_001_sorted.bam.bai
INDEXED_ALIGNED_BAM_FILE = f'{ALIGNED_BAM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{INDEXED_ALIGNED_BAM_FILE_EXTENSION}'

# Calbicans-1_S29_L006_001_sorted.g.vcf.gz
GVCFS_FILE = f'{GVCFS_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{GVCFS_FILE_EXTENSION}'

# Calbicans-1_S29_L006_001_sorted.vcf.gz
VCFS_FILE = f'{VCFS_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{VCFS_FILE_EXTENSION}'

#------------------------------------------------------------
# File lists
#
# Now we can use the single file patterns in conjunction with
# glob_wildcards and expand to build the lists of all expected files.

# the list of Read 1 and 2 file names 
READ_ONE_PREFIX,READ_ONE_SUFFIX = glob_wildcards(READ_ONE_FILE)
READ_TWO_PREFIX,READ_TWO_SUFFIX = glob_wildcards(READ_TWO_FILE)

SAMPLE_LIST = expand(SAMPLE, reads_file_prefix=READ_ONE_PREFIX)
READ_ONE_FILES = expand(READ_ONE_FILE, reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
READ_TWO_FILES = expand(READ_TWO_FILE, reads_file_prefix=READ_TWO_PREFIX, reads_file_suffix=READ_TWO_SUFFIX)

#READONE_FILES = glob_wildcards(READONE_FILE).reads_file_prefix+{READ1_TAG}+glob_wildcards(READONE_FILE).reads_file_suffix+{READS_FILE_EXTENTION}
#READTWO_FILES = glob_wildcards(READTWO_FILE).reads_file_prefix #+{READ1_TAG}+glob_wildcards(READTWO_FILE).reads_file_suffix+{READS_FILE_EXTENTION}
#TEST_FILES=expand(TEST_FILE,reads_file_prefix=READONE_PREFIX, reads_file_suffix=READONE_SUFFIX)

PREQC_ONE_FILES=expand(PREQC_ONE_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
PREQC_TWO_FILES=expand(PREQC_TWO_FILE,reads_file_prefix=READ_TWO_PREFIX, reads_file_suffix=READ_TWO_SUFFIX)

ADAPTRIM_ONE_FILES=expand(ADAPTRIM_ONE_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
ADAPTRIM_TWO_FILES=expand(ADAPTRIM_TWO_FILE,reads_file_prefix=READ_TWO_PREFIX, reads_file_suffix=READ_TWO_SUFFIX)

ADAPTRIM_QC_ONE_FILES=expand(ADAPTRIM_QC_ONE_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
ADAPTRIM_QC_TWO_FILES=expand(ADAPTRIM_QC_TWO_FILE,reads_file_prefix=READ_TWO_PREFIX, reads_file_suffix=READ_TWO_SUFFIX)

ADAPTRIM_QUALTRIM_ONE_P_FILES=expand(ADAPTRIM_QUALTRIM_ONE_P_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
ADAPTRIM_QUALTRIM_ONE_U_FILES=expand(ADAPTRIM_QUALTRIM_ONE_U_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
ADAPTRIM_QUALTRIM_TWO_P_FILES=expand(ADAPTRIM_QUALTRIM_TWO_P_FILE,reads_file_prefix=READ_TWO_PREFIX, reads_file_suffix=READ_TWO_SUFFIX)
ADAPTRIM_QUALTRIM_TWO_U_FILES=expand(ADAPTRIM_QUALTRIM_TWO_U_FILE,reads_file_prefix=READ_TWO_PREFIX, reads_file_suffix=READ_TWO_SUFFIX)

ADAPTRIM_QUALTRIM_BASEOUT_FILES=expand(ADAPTRIM_QUALTRIM_BASEOUT_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)

ADAPTRIM_QUALTRIM_QC_ONE_FILES=expand(ADAPTRIM_QUALTRIM_QC_ONE_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
ADAPTRIM_QUALTRIM_QC_TWO_FILES=expand(ADAPTRIM_QUALTRIM_QC_TWO_FILE,reads_file_prefix=READ_TWO_PREFIX, reads_file_suffix=READ_TWO_SUFFIX)

#ADAPTRIM_QUALTRIM_QC_ONE_SCREEN_FILES=expand(ADAPTRIM_QUALTRIM_QC_ONE_SCREEN_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
#ADAPTRIM_QUALTRIM_QC_TWO_SCREEN_FILES=expand(ADAPTRIM_QUALTRIM_QC_TWO_SCREEN_FILE,reads_file_prefix=READ_TWO_PREFIX, reads_file_suffix=READ_TWO_SUFFIX)

FASTSCREEN_LOGFILES=expand(FASTSCREEN_LOGFILE,reads_file_prefix=READ_TWO_PREFIX, reads_file_suffix=READ_TWO_SUFFIX)

#ALIGNED_BAM_ONE_FILES=expand(ALIGNED_BAM_ONE_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
#ALIGNED_BAM_TWO_FILES=expand(ALIGNED_BAM_TWO_FILE,reads_file_prefix=READ_TWO_PREFIX, reads_file_suffix=READ_TWO_SUFFIX)
ALIGNED_BAM_FILES=expand(ALIGNED_BAM_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)

ALIGNED_BAM_FLAGSTAT_FILES=expand(ALIGNED_BAM_FLAGSTAT_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)

#MERGED_BAM_ONE_FILES=expand(ALIGNED_BAM_ONE_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
#MERGED_BAM_TWO_FILES=expand(ALIGNED_BAM_TWO_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)

INDEXED_ALIGNED_BAM_FILES=expand(INDEXED_ALIGNED_BAM_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)

GVCFS_FILES=expand(GVCFS_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
VCFS_FILES=expand(VCFS_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)

# The list of all sam files
# BOOK_FILE = f'{INPUT_DIR}{{book}}.txt'
# DAT_FILE = f'{DAT_DIR}{{book}}.dat'
# BOOK_NAMES = glob_wildcards(BOOK_FILE).book
# ALL_DATS = expand(DAT_FILE, book=BOOK_NAMES)
#READONE_SAM_FILES = expand(READONE_SAM_FILE, sam_file_name=READONE_FILES)
#READTWO_SAM_FILES = expand(READTWO_SAM_FILE, sam_file_name=READONE_FILES)

#------------------------------------------------------------
# Rules
#
# Note that when using this pattern, it is rare for a rule to
# define filename patterns directly. Nearly all inputs and outputs
# can be specified using the existing global variables.

# pseudo-rule that tries to build everything.
# Just add all the final outputs that you want built.
rule all:
    input: PREQC_ONE_FILES, PREQC_TWO_FILES, \
    ADAPTRIM_QUALTRIM_QC_ONE_FILES, ADAPTRIM_QUALTRIM_QC_TWO_FILES, \
    ADAPTRIM_QC_ONE_FILES, ADAPTRIM_QC_TWO_FILES, \
    FASTSCREEN_LOGFILES, \
    ALIGNED_BAM_FILES, ALIGNED_BAM_FLAGSTAT_FILES, INDEXED_ALIGNED_BAM_FILES, VCFS_FILES

rule clean:
    shell: f'rm -rf {ALIGNED_BAM_DIR}'

# Prequality control
rule preqc:
    input: READ_ONE_FILE, READ_TWO_FILE
    output: PREQC_ONE_FILE, PREQC_TWO_FILE
    shell: 'fastqc {input} --outdir={PREQC_DIR}'

# Trim adapter sequences using BBduk
rule adaptrim:
    input: 
        in1=READ_ONE_FILE,
        in2=READ_TWO_FILE
    output: 
        out1=ADAPTRIM_ONE_FILE,
        out2=ADAPTRIM_TWO_FILE
    shell: 'bbduk.sh -Xmx1g in1={input.in1} in2={input.in2} out1={output.out1} out2={output.out2} ref={ADAPTERS} ktrim=r k=23 mink=11 hdist=1 tpe tbo' 


# Generate fastqc report from trimmed sequences
rule adaptrim_qc:
    input:
        in1=ADAPTRIM_ONE_FILE,
        in2=ADAPTRIM_TWO_FILE
    output:
        out1=ADAPTRIM_QC_ONE_FILE,
        out2=ADAPTRIM_QC_TWO_FILE
    shell: 'fastqc {input.in1} {input.in2} --outdir={ADAPTRIM_QC_DIR}'

# Trim back low quality scores at the end of reads with Trimmomatic
rule adaptrim_qualtrim:
    input:
        in1=ADAPTRIM_ONE_FILE,
        in2=ADAPTRIM_TWO_FILE
    output:
        out1=ADAPTRIM_QUALTRIM_ONE_P_FILE,
        out2=ADAPTRIM_QUALTRIM_ONE_U_FILE,
        out3=ADAPTRIM_QUALTRIM_TWO_P_FILE,
        out4=ADAPTRIM_QUALTRIM_TWO_U_FILE,
        #baseout=ADAPTRIM_QUALTRIM_BASEOUT_FILE
    #log:
        #{SNAKEMAKE_LOGS}ADAPTRIM_QUALTRIM.log
    #shell: 'trimmomatic PE {input.in1} {input.in2} -baseout {output.baseout} SLIDINGWINDOW:4:15 MINLEN:36'
    shell: 'trimmomatic PE {input.in1} {input.in2} {output.out1} {output.out2} {output.out3} {output.out4} SLIDINGWINDOW:4:15 MINLEN:36'

# Generate fastqc report from trimmomatic outputs
rule adaptrim_qualtrim_qc:
    input:
        in1 = ADAPTRIM_QUALTRIM_ONE_P_FILE,
        in2 = ADAPTRIM_QUALTRIM_TWO_P_FILE
    output:
        out1 = ADAPTRIM_QUALTRIM_QC_ONE_FILE,
        out2 = ADAPTRIM_QUALTRIM_QC_TWO_FILE
    shell: 'fastqc {input.in1} {input.in2} --outdir={ADAPTRIM_QUALTRIM_QC_DIR}'

# Use fastq screen to determine the amount of library contamination
# By default the program looks for a configuration file named “fastq_screen.conf” in the folder
# where the FastQ Screen script it is located. If you wish to specify a different configuration file,
# which may be placed in different folder, then use the –conf option:
rule fastq_screen:
    input:
        in1 = ADAPTRIM_QUALTRIM_QC_ONE_FILE,
        in2 = ADAPTRIM_QUALTRIM_QC_TWO_FILE
    log: FASTSCREEN_LOGFILE
#   output:
#       out1 = ADAPTRIM_QUALTRIM_QC_ONE_SCREEN_FILE,
#       out2 = ADAPTRIM_QUALTRIM_QC_TWO_SCREEN_FILE
    #shell: 'fastq_screen --conf {FASTQ_SCREEN_CONF} --aligner bowtie2 --outdir {SCREENING_DIR} {input.in1}'
    shell: 'fastq_screen --conf {FASTQ_SCREEN_CONF} --aligner bwa --outdir {SCREENING_DIR} {input.in1} >{log}'

rule indexRef:
    input: REF_FNA_FILE
    log: REF_FNA_INDEXING_LOG
    shell: 'bwa-mem2 index {input}'

# Align paired reads with bwa-mem2
rule align_reads:
    input:
        in1 = ADAPTRIM_QUALTRIM_ONE_P_FILE,
        in2 = ADAPTRIM_QUALTRIM_TWO_P_FILE
    output: ALIGNED_BAM_FILE
    params: SAMPLE
    shell: 'bash {READ_GROUP_SCRIPT} {REF_FNA_FILE} {input.in1} {input.in2} {params} {output}'

# Calculate alignment statistics with samtools
rule align_stats:
    input: ALIGNED_BAM_FILE
    output: ALIGNED_BAM_FLAGSTAT_FILE
    shell: 'samtools flagstat {input} > {output}'

# if running alignment jobs in parallel, distributing horizontally across hundreds of nodes
# rather than trying to run a single job with dozens of cores, sorted BAM files from all nodes
# need to be merged together for further downstream analysis.

# Merge filtered bam files with samtools merge (done manually)

# Mark and remove duplicates with Picard

# Duplicate marking and removal
# this rule is only needed if alignment is done across multiple computer nodes

rule picard_markdup:
	input:
		bam = (config["merged_bams_path"] + "{sample}_sorted_merged.bam")
	output:
		(config["picard_report_path"] + "{sample}_sorted_merged_DuplicationMetrics.txt"),
		(config["merged_bams_path"] + "{sample}_sorted_merged_rmdups.bam")
	params:
		outfile = (config["merged_bams_path"] + "{sample}_sorted_merged_rmdups.bam"),
		metrics = (config["picard_report_path"] + "{sample}_sorted_merged_DuplicationMetrics.txt"),
		tmp_dir = config["merged_bams_path"],
		to_rm = (config["merged_bams_path"] + "{sample}_sorted_merged.bam")
	shell:
       """
		set -euxo pipefail
		gatk MarkDuplicates --REMOVE_DUPLICATES--TMP_DIR {params.tmp_dir} -I {input.bam} -O {params.outfile} -M {params.metrics}
		if [ -s {params.outfile} ]
		then
			rm {params.to_rm}
		fi		
       """

#---------------------------------------
# Variant calling with HaplotypeCaller
#---------------------------------------

# Index bam files with samtools
rule index_bam:
    input: ALIGNED_BAM_FILE
    output: INDEXED_ALIGNED_BAM_FILE
    shell: 'samtools index {input} {output}'

# Calculate genotype likelihoods with HaplotypeCaller
# The GATK uses two files to access and safety check access to the reference files: a .dict dictionary of the contig names and sizes
# and a .fai fasta index file to allow efficient random access to the reference bases.
# You have to generate these files in order to be able to use a Fasta file as reference.
# samtools dict GCF_000182965.3_ASM18296v3_genomic.fa -o GCF_000182965.3_ASM18296v3_genomic.dict
# samtools faidx GCF_000182965.3_ASM18296v3_genomic.fa -o GCF_000182965.3_ASM18296v3_genomic.fa.fai
rule create_ref_dict:
    input: REF_FNA_FILE
    output: REF_FNA_DICT
    shell: 'samtools dict {input} -o {output}'

rule create_ref_fai:
    input: REF_FNA_FILE
    output: REF_FNA_FAI
    shell: 'samtools faidx {input} -o {output}'

rule haplotypecaller:
    input: 
        bam = ALIGNED_BAM_FILE,
        idx = INDEXED_ALIGNED_BAM_FILE
    output: GVCFS_FILE
    shell: 'gatk HaplotypeCaller --TMP_DIR {GVCFS_DIR} -R {REF_FNA_FILE} -I {input.bam} -ERC GVCF -O {output}'

# Joint call genotypes for each sample
rule joint_genotype:
    input: GVCFS_FILE
    output: VCFS_FILE
    shell: 'gatk GenotypeGVCFs --TMP_DIR {VCFS_DIR} -R {REF_FNA_FILE} -V {input} -O {output}'

