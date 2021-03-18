AmpSeQC: Amplicon Sequencing Quality Control
============================================

A general multiplexed amplicon sequencing and quality control pipeline specifically built
for _Plasmodium falciparum_ data

**Input**: paired-end fastq files, reference genome, gff3 annotaiton file of amplicon panel or genes

**Output**: tsv file of read counts per amplicon/gene per sample, MultiQC reports

### Steps in pipeline
1. Search for and validate specified paired-end fastq files
2. Run FastQC on raw reads.
3. Run Trim Galore, removing adapter content and low quality sequences.
4. Run FastQC on post-QC reads.
5. Align to reference genome using BWA-Mem or Bowtie2.
6. Filter for soft-clipping (BWA-Mem only) and retain only properly paired reads (both aligners).
7. Run FastQC on bam file of filtered aligned reads.
9. Once all samples finished QC and aligning, generate table of read counts per amplicon per sample
10. Run MultiQC to collate all FastQC results into reports.

Dependencies
------------

* Python>=3.7
* Bowtie2
* BWA
* FastQC
* HTSeq
* MultiQC
* Samclip
* Samtools
* Trim-Galore

Installation
------------

1. Install Anaconda or miniconda (if not already installed)
2. Clone Github repository:
    
    `git clone https://github.com/broadinstitute/AmpSeQC.git`

3. Move into the directory:
    
    `cd AmpSeQC`

4. Create a new conda environment:
    
    `conda env create -f environment.yml`

5. Activate the environment:

    `source activate ampseqc`

Quick start guide
-----------------

**Make sure you have activated the ampseqc conda environment `source activate ampseqc`!**

You can run the script like so.

`python3 /path/to/AmpSeQC/AmpSeQC.py -c output_read_counts.tsv -r /path/to/reference.fasta -a /path/to/annotations.gff fastq_R1.fq.gz fastq_R2.fq.gz ...`

where
* `output_read_counts.tsv` is the file you want your read counts to go to
* `/path/to/reference.fasta` is a samtools and bwa or bowtie2 indexed fasta file of your reference
  - default is reference.fasta in script directory
* `/path/to/annotations.gff` is a gff3 file of your amplicons/genes you want read counts for
  - default is amplicons.gff in script directory
* `fastq_R1.fq.gz` and `fastq_R2.fq.gz` are paired-end reads in fastq format (you can list as many pairs of files as you want)

There's a lot more to the script, but those are the basics. See help file for more details.

More details on parameters
--------------------------
```
usage: AmpSeQC.py [options] -c [read_counts] -r [ref.fasta] -a [annot.gff] fastq [fastq ...]

Amplicon sequencing quality control pipeline

positional arguments:
  fastq                 Paired-end fastq files (can be gzipped). Expects pairs to be in separate files.

optional arguments:
  -h, --help            show this help message and exit
  -c COUNTS, --counts COUNTS
                        Read count tsv file (default: read_counts.tsv)
  -r REF, --ref REF     Indexed reference fasta file to align to (default: reference.fasta)
  -a ANNOT, --annot ANNOT
                        Amplicon/gene gff3 file to generate read counts of (default: amplicons.gff)
  --2color              Run with 2 color chemistry (e.g. iSeq) QC parameters (default: False
  -l MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum read length (in bp) to retain after trimming (default: 70)
  -q MIN_BQ, --min_bq MIN_BQ
                        Minimum base quality (PHRED score) to retain after trimming (default: 20)
  -N MAX_N, --max_N MAX_N
                        Maximum number of 'N' bases allowed in read (default: 1)
  -I MAX_INSERT_SIZE, --max_insert_size MAX_INSERT_SIZE
                        Maximum insert size (in bp) expected for aligner (default: 500)
  -S SOFT_CLIP, --soft_clip SOFT_CLIP
                        Maximum soft clipping (in bp) allowed for BWA-Mem (default: 5)
  --bowtie2             Align with Bowtie2 instead of BWA-Mem (default: False)
  --min_amplicon_count MIN_AMPLICON_COUNT
                        Minimum total read count to retain an amplicon (default: 0)
  --min_sample_count MIN_SAMPLE_COUNT
                        Minimum total read count to retain a sample (default: 0)
  --no_fastqc           Do not run FastQC or MultiQC
  -p PROCS, --procs PROCS
                        Number of processors to use (default: 1)
```

### Input/output
* `-c COUNTS, --counts COUNTS`
  - This is your output file where you want read counts per amplicon per sample to go to.
* `-r REF, --ref REF`
  - This is a samtools faidx and bwa and/or bowtie2 indexed (depending on which aligner you used) reference genome.
  - The default is reference.fasta, the _Plasmodium falciparum_ 3D7 genome.
* `-a ANNOT, --annot ANNOT`
  - This is the annotation file with your amplicons or genes. Needs to be gff3 format!
  - The default is amplicons.gff, the Neafsey lab's own _Plasmodium_ amplicon panel.

### QC parameters
* `--2color`
  - This runs Trim-Galore in 2-color mode, useful for iSeq or NovaSeq 2-color chemistry data. (default: disabled)
* `-l MIN_LENGTH, --min_length MIN_LENGTH`
  - This specifies the minimum read length to retain after trimming with Trim-Galore. (default 70 bp)
* `-q MIN_BQ, --min_bq MIN_BQ`
  - This specifies the minimum base quality to retain during trimming. (default: Q20)
* `-N MAX_N, --max_N MAX_N`
  - This specifies the number of "N" bases allowed in a given read during trimming. (default: 1)

### Alignment parameters
* `-I MAX_INSERT_SIZE, --max_insert_size MAX_INSERT_SIZE`
  - This specifies the expected maximum insert size of your library. Set this to a bit bigger than your largest amplicon. (default: 500 bp)
* `-S SOFT_CLIP, --soft_clip SOFT_CLIP`
  - This specifies the amount of soft clipping allowed in alignment. Only applies when using BWA Mem. (default: 5 bp)
* `--bowtie2`
  - Run Bowtie2 global aligner instead of BWA Mem local aligner. No soft clipping is allowed. (default: disabled)

### Read count output filtering
* `--min_amplicon_count MIN_AMPLICON_COUNT`
  - Filter read count output to remove amplicons with fewer than this many reads across all samples (default: 0)
* `--min_sample_count MIN_SAMPLE_COUNT`
  - Filter samples with fewer than this many total counts across all amplicons (default: 0)

### Other
* `--no_fastqc`
  - Don't run FastQC or MultiQC. Speeds up analysis if only interested in read counts. (default: disabled)
* `-p PROCS, --procs PROCS`
  - Enabled parallel processing of QC and alignment of samples. Specify up to the number of processors your system has. (default: 1)

Citation
--------
Coming soon...