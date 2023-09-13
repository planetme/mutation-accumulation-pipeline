# mutation-accumulation-pipeline
A pipeline specifically for the purpose of processing small sequences generated from mutation accumulation experiments.
Credit to Nathaniel Sharp and Jacob Roman-Fredette for the general snakemake rule layout and read group shell script.

This pipeline currently does not support multiple lane sequencing runs for now.

## Instructions
1. Install the following packages using conda:
    
  fastqc
  BBMap
  Trimmomatic
  fastq_screen
  bwa-mem2
  samtools
  picard
  gatk4

2. In your working directory, create a `pipelines` directory and put the config, readgroup shellscript, fastq_screen config, and pipeline config in there.
3. Create a directory called `sequences` and move your raw sequencing outputs into a new directory underneath named `reads`.
4. Create another directory under `sequences` named `adapter`, and move the adapters that you used under there.
5. Next, download your species' reference genome under a new directory titled `alignment_reference`.
6. If you want to screen against suspected contaminants, download reference genomes into a directory named `screens`.
7. Run the pipeline.
