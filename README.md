# mutation-accumulation-pipeline
A pipeline specifically to process small sequences generated from mutation accumulation experiments.
Credit to Nathaniel Sharp and Jacob Roman-Fredette for the general snakemake rule layout and read group shell script.

This pipeline currently does not support multiple-lane sequencing runs for now.

## Dependencies
fastqc, BBMap, Trimmomatic, fastq_screen, bwa-mem2, samtools, picard, gatk4

## Instructions
1. Install the dependencies using conda.
2. In your working directory, create a `pipelines` directory and put the pipeline config, readgroup shell script, your fastq_screen config, and the pipeline config in there.
3. Create a directory called `sequences` and move your raw sequencing outputs into a new directory underneath named `reads`.
4. Create another directory under `sequences` named `adapter`, and move the adapters that you used under there.
5. Next, download your species' reference genome under a new directory titled `alignment_reference`.
6. If you want to screen against suspected contaminants, download reference genomes into a directory named `screens`.
7. Run the pipeline on -np mode to check if the expected files will be generated.
8. Run the pipeline with the following command: snakemake -s pipeline.smk --cores <number>.
