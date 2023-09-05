import pathlib as pathlib
import glob as glob
import os as os

configfile: "config.json"

workingPath = pathlib.Path(os.path.dirname(__file__))
outputsPath = workingPath.joinpath('outputs')

# Specify file output directories
outputdirs = [
    "fastqc_reports",
    "samtools_reports",
    "aligned_bams",
    "merged_bams",
    "gvcfs"
    "vcfs"
]

# Specify quality control reports (to be placed under outputsPath/fastqc_reports)
qcdirs = [
    "preqc",
    "adaptrim",
    "qualtrim",
    "screening",
    "picard"
]

# Generate output folders
def generateOutputDirs():
    for outputdir in outputdirs:
        outputsPath.mkdir(outputdir)
    for qcdir in qcdirs:
        outputsPath.joinpath('fastqc_reports').mkdir(qcdir)

# Get inputs
def getSamples():
	samples = list()

	glob.glob((inputsPath + "/*S"))

def getLanes():

# This function should return a dictionary indexed by samples containing both reads. Further processing allows for arbitrary separation of each read regardless of file naming scheme into separate lists for smk rules.
def mapPairs(inputsPath, inputExtension, wildcard):
	pairsDict = dict()

	for sampleNumber in range(config["sample_starting_number"], config["sample_ending_number"]):
		reads = glob.glob((inputsPath + "/" + str(wildcard[0]) + "S" + str(sampleNumber) + str(wildcard[1]) + inputExtension))
		pairsDict[sampleNumber] = reads

	return pairsDict

# Aforementioned further processing
def returnInputs(inputsPath, inputExtension, wildcard):
	r1 = list()
	r2 = list()
	pairsDict = mapPairs(inputsPath, inputExtension, wildcard)

	for pair in pairsDict:
		r1.append(os.path.basename(pair[0]))
		r2.append(os.path.basename(pair[1]))
	
	return [r1, r2]

def getInitialInputsR1(wildcard):
	return [x for x in returnInputs(config["inputs_path"], config["input_extension"])[0], wildcard]
def getInitialInputsR2(wildcard):
	return [x for x in returnInputs(config["inputs_path"], config["input_extension"])[1], wildcard]

def getFastQCTargets():
	r1_targets = list()
	r2_targets = list()

	for s in config["samples"]:
		for l in config["lanes"]:
			r1_targets.append([str(config["fastqc_report_path"] + "preqc/" + os.path.basename(x).split(".")[0] + "_fastqc.zip") for x in glob.glob((config["fastq_path"] + str(s) + "_*_" + str(l) + "_1" + ".fq.gz"))])
			r2_targets.append([str(config["fastqc_report_path"] + "preqc/" + os.path.basename(x).split(".")[0] + "_fastqc.zip") for x in glob.glob((config["fastq_path"] + str(s) + "_*_" + str(l) + "_2" + ".fq.gz"))])
	return r1_targets,r2_targets

# This rule designates the desired outputs, snakemake will work to generate these outputs using necessary rules
rule all:
	input:
		# Pre-QC reports
		getFastQCTargets(),
		# Adapter-trimmed files from BBduk
		expand((config["adaptrim_path"] + "{sample}_{lane}.adaptrim.1P.fq.gz"), sample = config["samples"], lane = config["lanes"]),
		expand((config["adaptrim_path"] + "{sample}_{lane}.adaptrim.2P.fq.gz"), sample = config["samples"], lane = config["lanes"]),
		# QC reports on BBduk output
		expand((config["fastqc_report_path"] + "adaptrim/" + "{sample}_{lane}.adaptrim.1P_fastqc.zip"), sample = config["samples"], lane = config["lanes"]),
		expand((config["fastqc_report_path"] + "adaptrim/" + "{sample}_{lane}.adaptrim.2P_fastqc.zip"), sample = config["samples"], lane = config["lanes"]),
		# Adapter- and quality-trimmed files from Trimmomatic
		expand((config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_1P.fq.gz"), sample = config["samples"], lane = config["lanes"]),
		expand((config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_1U.fq.gz"), sample = config["samples"], lane = config["lanes"]),
		expand((config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_2P.fq.gz"), sample = config["samples"], lane = config["lanes"]),
		expand((config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_2U.fq.gz"), sample = config["samples"], lane = config["lanes"]),
		# QC reports on Trimmomatic output
		expand((config["fastqc_report_path"] + "adaptrim_qualtrim/" + "{sample}_{lane}.adaptrim.qualtrim_1P_fastqc.zip"), sample = config["samples"], lane = config["lanes"]),
		expand((config["fastqc_report_path"] + "adaptrim_qualtrim/" + "{sample}_{lane}.adaptrim.qualtrim_2P_fastqc.zip"), sample = config["samples"], lane = config["lanes"]),
		# FastQC-Screen outputs
		expand((config["fastq_screen_report_path"] + "{sample}_{lane}.adaptrim.qualtrim_1P_screen.html"), sample = config["samples"], lane = config["lanes"]),
		expand((config["fastq_screen_report_path"] + "{sample}_{lane}.adaptrim.qualtrim_2P_screen.html"), sample = config["samples"], lane = config["lanes"]),
		# BWA-MEM2 outputs
		expand((config["aligned_bams_path"] + "{sample}_{lane}_sorted.bam"), sample = config["samples"], lane = config["lanes"]),
		# Samtools flagstat outputs
		expand((config["samtools_report_path"] + "{sample}_{lane}_sorted.flagstat.txt"), sample = config["samples"], lane = config["lanes"]),
		# Samtools merge outputs
		#expand((config["merged_bams_path"] + "{sample}_sorted_merged.bam"), sample = config["samples"]),
		# Picard MarkDuplicates outputs
		expand((config["picard_report_path"] + "{sample}_sorted_merged_DuplicationMetrics.txt"), sample = config["samples"]),
		expand((config["merged_bams_path"] + "{sample}_sorted_merged_rmdups.bam"), sample = config["samples"]),
		# Indexed bams
		expand((config["merged_bams_path"] + "{sample}_sorted_merged_rmdups.bam.bai"), sample = config["samples"]),
		# HaplotypeCaller outputs
		expand((config["rd_gvcfs_path"] + "{sample}_sorted_merged_rmdups.g.vcf.gz"), sample = config["samples"]),
		# GenotypeGVCFs outputs
		expand((config["rd_vcfs_path"] + "{sample}_sorted_merged_rmdups.vcf.gz"), sample = config["samples"])

# Prequality control
rule preqc:
	input:
	output:
		(outputsPath.joinpath("fastqc_reports/preqc") + "{wildcard}" + "fastqc.zip"),
        (outputsPath.joinpath("fastqc_reports/preqc") + "{wildcard}" + "fastqc.zip")
	params:
		fastq_r1 = getInitialInputsR1
        fastq_r2 = getInitialInputsR2
		indir = config["inputs_path"],
		outdir = outputsPath.joinpath("fastqc_reports/preqc")
		wkdir = workingPath
	shell:
		"""
		cd {params.indir}
		fastqc {params.fastq_r1} {params.fastq_r2} --outdir={params.outdir}
		cd {params.wkdir}
		"""


# Trim adapter sequences using BBduk
rule adaptrim:
	input:
	output:
		(config["adaptrim_path"] + "{wildcard}_R1_001" + ".adaptrim.1P.fq.gz"),
		(config["adaptrim_path"] + "{wildcard}_R2_001" + ".adaptrim.2P.fq.gz")
	params:
		indir = config["inputs_path"],
		wkdir = workingPath
		input1 = getBBdukInputsR1,
		input2 = getBBdukInputsR2,
		output1 = (outputsPath.joinpath("fastqc/adaptrim") + "{sample}_{lane}.adaptrim.1P.fq.gz"),
		output2 = (outputsPath.joinpath("fastqc/adaptrim") + "{sample}_{lane}.adaptrim.2P.fq.gz"),
		adapters = config["adapters"]
	shell:
		"""
		cd {params.indir}
		bbduk.sh -Xmx1g in1={params.input1} in2={params.input2} out1={params.output1} out2={params.output2} ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tpe tbo 
		cd {params.wkdir}
		"""

# Generate fastqc report from trimmed sequences
rule adaptrim_qc:
	input:
		input1 = (config["adaptrim_path"] + "{sample}_{lane}.adaptrim.1P.fq.gz"),
		input2 = (config["adaptrim_path"] + "{sample}_{lane}.adaptrim.2P.fq.gz")
	output:
		(config["fastqc_report_path"] + "adaptrim/" + "{sample}_{lane}.adaptrim.1P_fastqc.zip"),
		(config["fastqc_report_path"] + "adaptrim/" + "{sample}_{lane}.adaptrim.2P_fastqc.zip")
	params:
		indir = config["adaptrim_path"],
		outdir = config["fastqc_report_path"] + "adaptrim/",
		wkdir = config["working_path"]
	shell:
		"""
		cd {params.indir}
		fastqc {input.input1} {input.input2} --outdir={params.outdir}
		cd {params.wkdir}
		"""

# Trim back low quality scores at the end of reads with Trimmomatic
rule adaptrim_qualtrim:
	input:
		input1 = (config["adaptrim_path"] + "{sample}_{lane}.adaptrim.1P.fq.gz"),
		input2 = (config["adaptrim_path"] + "{sample}_{lane}.adaptrim.2P.fq.gz")
	output:
		(config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_1P.fq.gz"),
		(config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_1U.fq.gz"),
		(config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_2P.fq.gz"),
		(config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_2U.fq.gz")
	params:
		output = (config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim.fq.gz")
	shell:
		"trimmomatic PE {input.input1} {input.input2} -baseout {params.output} SLIDINGWINDOW:4:15 MINLEN:36"

# Generate fastqc report from trimmomatic outputs
rule adaptrim_qualtrim_qc:
	input:
		input1 = (config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_1P.fq.gz"),
		input2 = (config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_2P.fq.gz")
	output:
		(config["fastqc_report_path"] + "adaptrim_qualtrim/" + "{sample}_{lane}.adaptrim.qualtrim_1P_fastqc.zip"),
		(config["fastqc_report_path"] + "adaptrim_qualtrim/" + "{sample}_{lane}.adaptrim.qualtrim_2P_fastqc.zip")
	params:
		indir = config["adaptrim_qualtrim_path"],
		outdir = config["fastqc_report_path"] + "adaptrim_qualtrim/",
		wkdir = config["working_path"]
	shell:
		"""
		cd {params.indir}
		fastqc {input.input1} {input.input2} --outdir={params.outdir}
		cd {params.wkdir}
		"""

# Use fastq screen to determine the amount of library contamination
rule fastq_screen:
	input:
		input1 = (config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_{read}P.fq.gz")
	output:
		(config["fastq_screen_report_path"] + "{sample}_{lane}.adaptrim.qualtrim_{read}P_screen.html")
	params:
		outdir = config["fastq_screen_report_path"],
		config = config["fastq_screen_references"]
	shell:
		"fastq_screen --conf {params.config} --aligner bowtie2 --outdir {params.outdir} {input.input1}"

rule indexRef:
    input:
        (config["alignment_reference"])
    log:
        (config["alignment_reference"].parent() + config["alignment_reference"].name + ".log")
    shell:
        "bwa-mem2 index {input}"

#-----------------------------------
# Alignment and quality control
#-----------------------------------

# Align paired reads with bwa-mem2
rule align_reads:
	input:
		input1 = (config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_1P.fq.gz"),
		input2 = (config["adaptrim_qualtrim_path"] + "{sample}_{lane}.adaptrim.qualtrim_2P.fq.gz")
	output:
		(config["aligned_bams_path"] + "{sample}_{lane}_sorted.bam")
	params:
		ref = config["bwa_ref"],
		sample = "{sample}",
		output = (config["aligned_bams_path"] + "{sample}_{lane}_sorted.bam"),
		script = config["read_group"]
	shell:
		"bash {params.script} {params.ref} {input.input1} {input.input2} {params.sample} {params.output}"

# Calculate alignment statistics with samtools
rule align_stats:
	input:
		bam = (config["aligned_bams_path"] + "{sample}_{lane}_sorted.bam")
	output:
		(config["samtools_report_path"] + "{sample}_{lane}_sorted.flagstat.txt")
	params:
		output = (config["samtools_report_path"] + "{sample}_{lane}_sorted.flagstat.txt")
	shell:
		"samtools flagstat {input.bam} > {params.output}"

# Merge filtered bam files with samtools merge (done manually)

# Mark and remove duplicates with Picard

# Duplicate marking and removal
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
		gatk MarkDuplicates --TMP_DIR {params.tmp_dir} -I {input.bam} -O {params.outfile} -M {params.metrics}
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
	input:
		bam = (config["merged_bams_path"] + "{sample}_sorted_merged_rmdups.bam")
	output: 
		idx = (config["merged_bams_path"] + "{sample}_sorted_merged_rmdups.bam.bai")
	shell:
		"samtools index {input.bam}"

# Calculate genotype likelihoods with HaplotypeCaller
rule haplotypecaller:
	input:
		bam = (config["merged_bams_path"] + "{sample}_sorted_merged_rmdups.bam"),
		idx = (config["merged_bams_path"] + "{sample}_sorted_merged_rmdups.bam.bai")
	output: 
		(config["rd_gvcfs_path"] + "{sample}_sorted_merged_rmdups.g.vcf.gz")
	params:
		outfile = (config["rd_gvcfs_path"] + "{sample}_sorted_merged_rmdups.g.vcf.gz"),
		ref = config["reference"],
		tmp_dir = config["rd_gvcfs_path"]
	shell:
		"gatk HaplotypeCaller --tmp-dir {params.tmp_dir} -R {params.ref} -I {input.bam} -ERC GVCF -O {params.outfile}"

# Joint call genotypes for each sample
rule joint_genotype:
	input:
		(config["rd_gvcfs_path"] + "{sample}_sorted_merged_rmdups.g.vcf.gz")
	output:
		(config["rd_vcfs_path"] + "{sample}_sorted_merged_rmdups.vcf.gz")
	params:
		outfile = config["rd_vcfs_path"] + "{sample}_sorted_merged_rmdups.vcf.gz",
		ref = config["reference"],
		tmp_dir = config["rd_vcfs_path"]
	shell:
		"gatk GenotypeGVCFs --tmp-dir {params.tmp_dir} -R {params.ref} -V {input} -O {params.outfile}" 


