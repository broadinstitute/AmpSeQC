# README
An amplicon sequencing pipeline for *Plasmodium falciparum*

**Input**: directory of demuxed fastq files, tsv file of amplicon panel or gff3 file of genes, optional parameter file

**Output**: tsv file of read counts per amplicon/gene per sample

Contact [Tim Straub](mailto:tim.straub@broadinstitute.org) for any questions.

## Description
This set of scripts will run sequencing quality control, alignments (with bowtie2), and generate read count per amplicon (or gene) per sample. It will also run FastQC on pre-QC, post-QC, and aligned reads. It uses an UGER task array to parallelize QC and alignment per sample.

### Steps in pipeline
1. Search for and validate paired-end fastq files in `DEMUX_DIRECTORY`
2. Run FastQC on raw reads.
3. Run Trim Galore, removing adapter content and low quality sequences.
4. Run polyG tail trimming on trimmed reads, which is most applicable to iSeq data. (TODO: only run if iSeq is true)
5. Run FastQC on post-QC reads.
6. Align to 3D7 genome using Bowtie2 global aligner. Retain only properly paired reads.
7. Run FastQC on bam file of aligned reads.
8. Once all samples finished aligning, generate read counts using `bedtools multicov`.

## Quick start guide
**Please copy or symlink the main script `ampseq_pipeline.sh` to your working directory, then edit (only if you copy!) if you need.**

**IMPORTANT:** The pipeline expects your reads to be paired gzipped fastq files with the following naming convention: `{SAMPLE}_R[12]*.fastq.gz`, where `{SAMPLE}` is a unique sample name. If they are in a different naming scheme, please copy and edit the `ampseq_pipeline.sh` script. If they are not gzipped, please gzip them. Save space!

You can submit it as an UGER job script (e.g. `use UGER; qsub ampseq_pipeline.sh`), or run the script directly (e.g. `sh ampseq_pipeline.sh`). The pipeline relies upon a conda environment that contains all the necessary tools, so no need to worry about dependencies.

You can edit the script directly: there is a small section of parameters in the beginning. You can also specify a parameter file (see `default_params.txt` or the [section below](#pipeline-parameters)) for an explanation of each parameter. The script accepts the parameter file as an argument (e.g., `sh ampseq_pipeline.sh yourParams.txt`), or it will look for `params.txt` or `default_params.txt` file. Otherwise, it will use the values specified in the script itself.

## Pipeline parameters
### Options you probably should change
- `DEMUX_DIRECTORY`
  - top-level demultiplexed directory of gzipped-fastq files (either from sequencer or demux pipeline)
- `BED_FILE`
  - output read count tsv file

### Options specific to quality control
- `ISEQ`
  - change to `true` if you're analyzing iSeq/NextSeq 2-color chemistry default_params
- `MIN_LENGTH`
  - minimum read length to retain (default 70 bp)
- `MIN_BQ`
  - minimum base quality (PHRED score) to keep (default 2)
- `MAX_N`
  - maximum number of "N" bases allowed in read (default 1)

### Options related to amplicon/gff3 file
- `ANNOTATION_FILE`
  - tsv file with amplicons, plus their chrom, start, and end, or gff3 file of genes
- `ANNOTATION_HEADER`
  - header for output file (adjust for whether you use a tsv or gff3 file)

### Options you can likely leave alone
- `PARASITE_GENOME`
  - Bowtie2 indexed fasta file containing parasite genome (default `PlasmoDB-46_Pfalciparum3D7_Genome.fasta`). If using your own, please run bowtie2-build on the fasta file first.
- `MAX_INSERT_SIZE`
  - max insert size for bowtie2 (default 1000 bp)
- `TASK_MEM`
  - task memory (default 16g, which should be plenty)
