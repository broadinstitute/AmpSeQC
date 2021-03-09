# README
An amplicon sequencing pipeline for *Plasmodium falciparum*

**Input**: directory of demuxed fastq files, gff3 file of amplicon panel or genes, optional parameter file

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
8. Run MultiQC to collate all FastQC results.
9. Once all samples finished aligning, generate read counts using `bedtools multicov`.

## Quick start guide
**Please copy the two shell scripts `AmpSeQC.sh` and `AmpSeQC_task.sh` to your working directory, then you can edit the main script if you need to.**

**IMPORTANT:** The pipeline expects your reads to be paired gzipped fastq files with the following naming convention: `{SAMPLE}_R[12]*.fastq.gz`, where `{SAMPLE}` is a unique sample name. If they are in a different naming scheme, please copy and edit the `AmpSeQC.sh` script. If they are not gzipped, please gzip them. Save space!

You can run the script directly (e.g. `./AmpSeQC.sh`).

You can edit the script directly: there is a small section of parameters in the beginning. You can also specify a parameter file (see `AmpSeQC_params.txt` or the [section below](#pipeline-parameters)) for an explanation of each parameter. The script accepts the parameter file as an argument (e.g., `./AmpSeQC.sh yourParams.txt`), or it will look for `params.txt` or `AmpSeQC_params.txt` file. Otherwise, it will use the values specified in the script itself.

## Pipeline parameters
### Options you probably should change
- `DEMUX_DIRECTORY`
  - top-level demultiplexed directory of gzipped-fastq files (either from sequencer or demux pipeline)
- `COUNT_FILE`
  - output read count tsv file

### Options specific to quality control
- `ISEQ`
  - change to `true` if you're analyzing iSeq/NextSeq 2-color chemistry
- `MIN_LENGTH`
  - minimum read length to retain (default 70 bp)
- `MIN_BQ`
  - minimum base quality (PHRED score) to keep (default 2)
- `MAX_N`
  - maximum number of "N" bases allowed in read (default 1)
- `MAX_SOFTCLIP`
  - maximum soft clipping allowed for bwa mem alignment (default 20 bp)

### Options related to amplicon/gff3 file
- `ANNOTATION_FILE`
  - gff3 file of amplicons or genes (default `amplicon.gff`)

### Options you can likely leave alone
- `REF_GENOME`
  - BWA/Bowtie2 indexed fasta file containing reference genome (default `reference.fasta`). If using your own, please run bowtie2-build on the fasta file first.
- `MAX_INSERT_SIZE`
  - max insert size for bwa mem or bowtie2 (default 500 bp)
