import pathlib as pathlib

configfile: "config.json"

workingPath = pathlib.Path(os.path.abspath()).parent()
outputsPath = workingPath.joinpath('outputs')
pairsDict = dict()

# Specify file output directories
outputdirs = [
    "fastqc_reports",
    "samtools_reports",
    "aligned_bams",
    "merged_bams",
    "gvcfs"
    "vcfs"
]

# Specify quality controll reports (to be placed under outputsPath/fastqc_reports)
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
def mapPairs()


def getFastQCInputsR1(wildcard):
	return [os.path.basename(x) for x in glob.glob((config["inputs_path"] + str(wildcard) + "_R1_001" + config["input_extension"]))]
def getFastQCInputsR2(wildcard):
	return [os.path.basename(x) for x in glob.glob((config["inputs_path"] + str(wildcard) + "_R2_001" + config["input_extension"]))]

def getBBdukInputsR1(wildcard):
	return [os.path.basename(x) for x in glob.glob((config["fastq_path"] + str(wildcard[0]) + "_*_" + str(wildcard[1]) + "_1.fq.gz"))]
def getBBdukInputsR2(wildcard):
	return [os.path.basename(x) for x in glob.glob((config["fastq_path"] + str(wildcard[0]) + "_*_" + str(wildcard[1]) + "_2.fq.gz"))]

# Prequality control
rule preqc:
	input:
	output:
		(outputsPath.joinpath("fastqc_reports/preqc") + "{wildcard}_R1_001" + "fastqc.zip"),
        (outputsPath.joinpath("fastqc_reports/preqc") + "{wildcard}_R2_001" + "fastqc.zip")
	params:
		fastq_r1 = getFastQCInputsR1,
        fastq_r2 = getFastQCInputsR2
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
		(config["adaptrim_path"] + "{wildcard}_R1_001.adaptrim.1P.fq.gz"),
		(config["adaptrim_path"] + "{wildcard}_R2_001.adaptrim.2P.fq.gz")
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


