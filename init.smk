import pathlib as pathlib

configfile: "config.json"

# Inputs: 

# Specify file output directories
outputdirs = [
    "fastqc_reports",
    "samtools_reports",
    "aligned_bams",
    "merged_bams",
    "gvcfs"
    "vcfs"
]

qcdirs = [
    "preqc",
    "adaptrim",
    "qualtrim",
    "screening",
    "picard"
]

# Generate output folders
def generateOutputDirs()
    workingPath = pathlib.Path(os.path.abspath()).parent()
    outputsPath = workingPath.joinpath('outputs')

    for outputdir in outputdirs:
        outputsPath.mkdir(outputdir)
    for qcdir in qcdirs:
        outputsPath.joinpath('fastqc_reports').mkdir(qcdir)
shell:
    """
    cd config["working_directory"]
    mkdir outputs
    cd outputs
    """

for outputdir in outputdirs:
    shell:
        "mkdir outputdir"

shell:
    "cd fastqc_reports"

for qcdir in qcdirs
    shell:
        "mkdir qcdir"

# Get inputs


# Prequality control
rule preqc:
	input:
	output:
		(config["fastqc_report_path"] + "preqc/" + "{wildcard}" + "fastqc.zip"),
	params:
		fastq_r1 = getFastQCInputsR1,
		indir = config["fastq_path"],
		outdir = config["fastqc_report_path"] + "preqc/",
		wkdir = config["working_path"]
	shell:
		"""
		cd {params.indir}
		fastqc {params.fastq_r1} {params.fastq_r2} --outdir={params.outdir}
		cd {params.wkdir}
		"""




