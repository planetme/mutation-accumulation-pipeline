import pathlib as pathlib
import glob as glob

configfile: "config.yaml"

#use config.yaml to get the outputdirs instead for portability
REF_FNA_FILE = config['ref_fna_file']
REF_FNA_INDEXING_LOG = ['ref_fna_indexing_log']
READS_DIR = config['reads_dir']
FASTQC_REPORTS_DIR = config['fastqc_reports_dir']
SAM_DIR = config['sam_dir'] ###################
SAMTOOLS_REPORTS_DIR = config['samtools_reports_dir']
ALIGNED_BAM_DIR = config['aligned_bam_dir']
MERGED_BAM_DIR = config['merged_bam_dir']
VCFS_DIR = config['vcfs_dir']
GVCFS_DIR = config['gvcfs_dir']
PREQC_DIR = config['preqc_dir']
ADAPTRIM_DIR = config['adaptrim_dir']
ADAPTRIM_QC_DIR = config['adaptrim_qc_dir']
ADAPTRIM_QUALTRIM_DIR = config['adaptrim_qualtrim_dir']
ADAPTRIM_QUALTRIM_QC_DIR = config['adaptrim_qualtrim_qc_dir']
QUALTRIM_DIR = config['qualtrim_dir']
SCREENING_DIR = config['screening_dir']
FASTQ_SCREEN_REFERENCES = config['fastq_screen_references']
PICARD_DIR = config['picard_dir']

READ1_TAG = config['read1_tag']
READ2_TAG = config['read2_tag']
ADAPTRIM1_TAG = config['adaptrim1_tag']
ADAPTRIM2_TAG = config['adaptrim2_tag']
QUALTRIM_P_TAG = config['qualtrim_p_tag']
QUALTRIM_U_TAG = config['qualtrim_u_tag']
ADAPTRIM_QUALTRIM_1P_TAG = config['adaptrim_qualtrim_1p_tag']
ADAPTRIM_QUALTRIM_2P_TAG = config['adaptrim_qualtrim_2p_tag']
ADAPTERS = config['adapters']
READS_FILE_EXTENSION = config['reads_file_extension']
PREQC_FILE_EXTENSION = config['preqc_file_extension']
ADAPTRIM_FILE_EXTENSION = config['adaptrim_file_extension']
ADAPTRIM_QC_FILE_EXTENSION = config['adaptrim_qc_file_extension']
ADAPTRIM_QUALTRIM_FILE_EXTENSION = config['adaptrim_qualtrim_file_extension']
ADAPTRIM_QUALTRIM_QC_FILE_EXTENSION = config['adaptrim_qualtrim_qc_file_extension']
SAM_FILE_EXTENSION = config['sam_file_extension']


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
# a_1_S28_R1_001.fastq.gz
READ_ONE_FILE = f'{READS_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{READS_FILE_EXTENSION}'
READ_TWO_FILE = f'{READS_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{READS_FILE_EXTENSION}'
#TEST_FILE=f'{SAM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{SAM_FILE_EXTENSION}'

PREQC_ONE_FILE = f'{PREQC_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{PREQC_FILE_EXTENSION}'
PREQC_TWO_FILE = f'{PREQC_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{PREQC_FILE_EXTENSION}'

ADAPTRIM_ONE_FILE = f'{ADAPTRIM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{ADAPTRIM1_TAG}{ADAPTRIM_FILE_EXTENSION}'
ADAPTRIM_TWO_FILE = f'{ADAPTRIM_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{ADAPTRIM1_TAG}{ADAPTRIM_FILE_EXTENSION}'

ADAPTRIM_QC_ONE_FILE = f'{ADAPTRIM_QC_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{ADAPTRIM1_TAG}{ADAPTRIM_QC_FILE_EXTENSION}'
ADAPTRIM_QC_TWO_FILE = f'{ADAPTRIM_QC_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{ADAPTRIM1_TAG}{ADAPTRIM_QC_FILE_EXTENSION}'

ADAPTRIM_QUALTRIM_ONE_P_FILE = f'{ADAPTRIM_QUALTRIM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{QUALTRIM_P_TAG}{ADAPTRIM_QUALTRIM_FILE_EXTENSION}'
ADAPTRIM_QUALTRIM_ONE_U_FILE = f'{ADAPTRIM_QUALTRIM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{QUALTRIM_U_TAG}{ADAPTRIM_QUALTRIM_FILE_EXTENSION}'
ADAPTRIM_QUALTRIM_TWO_P_FILE = f'{ADAPTRIM_QUALTRIM_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{QUALTRIM_P_TAG}{ADAPTRIM_QUALTRIM_FILE_EXTENSION}'
ADAPTRIM_QUALTRIM_TWO_U_FILE = f'{ADAPTRIM_QUALTRIM_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{QUALTRIM_U_TAG}{ADAPTRIM_QUALTRIM_FILE_EXTENSION}'
#(config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim.fq.gz")
ADAPTRIM_QUALTRIM_BASEOUT_FILE = f'{ADAPTRIM_QUALTRIM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{ADAPTRIM_QUALTRIM_FILE_EXTENSION}'

ADAPTRIM_QUALTRIM_QC_ONE_FILE = f'{ADAPTRIM_QUALTRIM_QC_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{ADAPTRIM_QUALTRIM_1P_TAG}{ADAPTRIM_QUALTRIM_QC_FILE_EXTENSION}'
ADAPTRIM_QUALTRIM_QC_TWO_FILE = f'{ADAPTRIM_QUALTRIM_QC_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{ADAPTRIM_QUALTRIM_1P_TAG}{ADAPTRIM_QUALTRIM_QC_FILE_EXTENSION}'

ADAPTRIM_QUALTRIM_QC_ONE_SCREEN_FILE = f'{ADAPTRIM_QUALTRIM_QC_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{ADAPTRIM_QUALTRIM_1P_TAG}.html'
ADAPTRIM_QUALTRIM_QC_TWO_SCREEN_FILE = f'{ADAPTRIM_QUALTRIM_QC_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{ADAPTRIM_QUALTRIM_1P_TAG}.html'

SAM_ONE_FILE = f'{SAM_DIR}{{reads_file_prefix}}{READ1_TAG}{{reads_file_suffix}}{SAM_FILE_EXTENSION}'
SAM_TWO_FILE = f'{SAM_DIR}{{reads_file_prefix}}{READ2_TAG}{{reads_file_suffix}}{SAM_FILE_EXTENSION}'
#READONE_SAM_FILE = f'{SAM_DIR}{{sam_file_name}}{SAM_FILE_EXTENSION}'

#------------------------------------------------------------
# File lists
#
# Now we can use the single file patterns in conjunction with
# glob_wildcards and expand to build the lists of all expected files.

# the list of Read 1 and 2 file names
READ_ONE_PREFIX,READ_ONE_SUFFIX = glob_wildcards(READ_ONE_FILE)
READ_TWO_PREFIX,READ_TWO_SUFFIX = glob_wildcards(READ_TWO_FILE)

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
ADAPTRIM_QUALTRIM_TWO_P_FILES=expand(ADAPTRIM_QUALTRIM_TWO_P_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
ADAPTRIM_QUALTRIM_TWO_U_FILES=expand(ADAPTRIM_QUALTRIM_TWO_U_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)

ADAPTRIM_QUALTRIM_BASEOUT_FILES=expand(ADAPTRIM_QUALTRIM_BASEOUT_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)

ADAPTRIM_QUALTRIM_QC_ONE_FILES=expand(ADAPTRIM_QUALTRIM_QC_ONE_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
ADAPTRIM_QUALTRIM_QC_TWO_FILES=expand(ADAPTRIM_QUALTRIM_QC_TWO_FILE,reads_file_prefix=READ_TWO_PREFIX, reads_file_suffix=READ_TWO_SUFFIX)

ADAPTRIM_QUALTRIM_QC_ONE_SCREEN_FILES=expand(ADAPTRIM_QUALTRIM_QC_ONE_SCREEN_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
ADAPTRIM_QUALTRIM_QC_TWO_SCREEN_FILES=expand(ADAPTRIM_QUALTRIM_QC_TWO_SCREEN_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)

SAM_ONE_FILES=expand(SAM_ONE_FILE,reads_file_prefix=READ_ONE_PREFIX, reads_file_suffix=READ_ONE_SUFFIX)
SAM_TWO_FILES=expand(SAM_TWO_FILE,reads_file_prefix=READ_TWO_PREFIX, reads_file_suffix=READ_TWO_SUFFIX)



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
    #input: SAM_ONE_FILES, SAM_TWO_FILES, PREQC_ONE_FILES, PREQC_TWO_FILES, ADAPTRIM_ONE_FILES, ADAPTRIM_TWO_FILES, ADAPTRIM_QC_ONE_FILES, ADAPTRIM_QC_TWO_FILES
    input: SAM_ONE_FILES, SAM_TWO_FILES, PREQC_ONE_FILES, PREQC_TWO_FILES, \
    ADAPTRIM_QUALTRIM_QC_ONE_FILES, ADAPTRIM_QUALTRIM_QC_TWO_FILES, \
    ADAPTRIM_QUALTRIM_QC_ONE_SCREEN_FILES, ADAPTRIM_QUALTRIM_QC_TWO_SCREEN_FILES

rule clean:
    shell: f'rm -rf {SAM_DIR}'

"""
rule test:
    input: READ_ONE_FILE, READ_TWO_FILE
    output: SAM_ONE_FILE, SAM_TWO_FILE
    #output: READONE_SAM_FILE
    #shell: 'cp {input} {output}'
    shell: 'ls {input}; touch {output}'
"""

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
        in1=ADAPTRIM_QC_ONE_FILE,
        in2=ADAPTRIM_QC_TWO_FILE
    output:
        out1=ADAPTRIM_QUALTRIM_ONE_P_FILE,
        out2=ADAPTRIM_QUALTRIM_ONE_U_FILE,
        out3=ADAPTRIM_QUALTRIM_TWO_P_FILE,
        out4=ADAPTRIM_QUALTRIM_TWO_U_FILE,
        baseout=ADAPTRIM_QUALTRIM_BASEOUT_FILE
    shell: 'trimmomatic PE {input.in1} {input.in2} -baseout {output.baseout} SLIDINGWINDOW:4:15 MINLEN:36'

# Generate fastqc report from trimmomatic outputs
rule adaptrim_qualtrim_qc:
    input:
        in1 = ADAPTRIM_QUALTRIM_ONE_P_FILE,
        in2 = ADAPTRIM_QUALTRIM_TWO_P_FILE
    output:
        out1 = ADAPTRIM_QUALTRIM_QC_ONE_FILE,
        out2 = ADAPTRIM_QUALTRIM_QC_TWO_FILE
    shell: 'fastqc {input.in1} {input.in2} --outdir={ADAPTRIM_QUALTRIM_DIR}'

# Use fastq screen to determine the amount of library contamination
rule fastq_screen:
    input:
        #input1 = (config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_{read}P.fq.gz")
        in1 = ADAPTRIM_QUALTRIM_QC_ONE_FILE,
        in2 = ADAPTRIM_QUALTRIM_QC_TWO_FILE
    output:
        #(config["fastq_screen_report_path"] + "{sample}_{lane}.adaptrim.qualtrim_{read}P_screen.html")
        out1 = ADAPTRIM_QUALTRIM_QC_ONE_SCREEN_FILE,
        out2 = ADAPTRIM_QUALTRIM_QC_TWO_SCREEN_FILE
    shell: 'fastq_screen --conf {FASTQ_SCREEN_REFERENCES} --aligner bowtie2 --outdir {SCREENING_DIR} {input.in1}'

rule indexRef:
    input: REF_FNA_FILE
    log: REF_FNA_INDEXING_LOG
    shell: 'bwa-mem2 index {input}'

