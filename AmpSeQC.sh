#!/bin/bash

##########################################
#### PIPELINE FOR AMPLICON SEQUENCING ####
####  CONTACT TIM STRAUB FOR DETAILS  ####
####   tim.straub@broadinstitute.org  ####
####                                  ####
####      QC --> ALIGN --> COUNT      ####
####        UPDATED MARCH 2021        ####
##########################################

#################################
#### STEP 0: Read parameters ####
#################################

# DEFAULT VALUES IF YOU DON'T SUPPLY A PARAMETER FILE
ISEQ=false
BWA=true
DEMUX_DIRECTORY="./"
BED_FILE="read_counts.tsv"
REF_GENOME="reference.fasta"
ANNOTATION_FILE="amplicons.gff"
MIN_LENGTH=70
MIN_BQ=2
MAX_INSERT_SIZE=500
MAX_SOFTCLIP=20
MAX_N=1
NTHREADS=1

# check for passing parameter file in arguments
# if not, then check for params.txt and default_params.txt files
# if not, then use default parameters
if [[ -s $1 ]]; then
    PARAM_FILE="$1"
elif [[ -s params.txt ]]; then
    PARAM_FILE="params.txt"
elif [[ -s AmpSeQC_params.txt ]]; then
    PARAM_FILE="AmpSeQC_params.txt"
else
    echo "WARNING: no parameter file specified or found. Using default parameters and checking for data in current directory."
fi

# read in parameter file
if [[ -s $PARAM_FILE ]]; then
    echo "INFO: reading in parameters from $PARAM_FILE"
    while read line; do
        if [[ ${line::1} == "#" ]]; then
            true # commented line
        elif [[ $line =~ ^\s*$ ]]; then
            true # blank line
        else
            KEY="$(echo "$line" | awk '{print $1}')"
            VALUE="$(echo "$line" | awk '{print $2}')"
            case "$KEY" in
                "ISEQ")
                    if [[ $VALUE =~ ^[Tt] ]]; then
                        ISEQ=true
                    else
                        ISEQ=false
                    fi
                    ;;
                "BWA")
                    if [[ $VALUE =~ ^[Tt] ]]; then
                        BWA=true
                    else
                        ISEQ=false
                    fi
                    ;;
                "DEMUX_DIRECTORY")
                    if [[  -d $VALUE ]]; then
                        DEMUX_DIRECTORY="$VALUE"
                    else
                        echo "WARNING: DEMUX_DIRECTORY '$VALUE' does not exist!"
                    fi
                    ;;
                "BED_FILE")
                    if [[ -s $VALUE ]]; then
                        echo "WARNING: $VALUE already exists... will overwrite."
                    fi
                    BED_FILE="$VALUE"
                    ;;
                "REF_GENOME")
                    if [[ -s $VALUE ]]; then
                        REF_GENOME="$VALUE"
                    else
                        echo "WARNING: REF_GENOME '$VALUE' is invalid!"
                    fi
                    ;;
                "ANNOTATION_FILE")
                    if [[ -s $VALUE ]]; then
                        ANNOTATION_FILE="$VALUE"
                    else
                        echo "WARNING: ANNOTATION_FILE '$VALUE' does not exist!"
                    fi
                    ;;
                "MIN_LENGTH")
                    if [[ $VALUE =~ ^[0-9]+$ ]]; then
                        MIN_LENGTH="$VALUE"
                    else
                        echo "WARNING: MIN_LENGTH '$VALUE' is invalid!"
                    fi
                    ;;
                "MIN_BQ")
                    if [[ $VALUE =~ ^[0-9]+$ ]]; then
                        MIN_BQ="$VALUE"
                    else
                        echo "WARNING: MIN_BQ '$VALUE' is invalid!"
                    fi
                    ;;
                "MAX_N")
                    if [[ $VALUE =~ ^[0-9]+$ ]]; then
                        MAX_N="$VALUE"
                    else
                        echo "WARNING: MAX_N '$VALUE' is invalid!"
                    fi
                    ;;
                "MAX_INSERT_SIZE")
                    if [[ $VALUE =~ ^[0-9]+$ ]]; then
                        MAX_INSERT_SIZE="$VALUE"
                    else
                        echo "WARNING: MAX_INSERT_SIZE '$VALUE' is invalid!"
                    fi
                    ;;
                "MAX_SOFTCLIP")
                    if [[ $VALUE =~ ^[0-9]+$ ]]; then
                        MAX_SOFTCLIP="$VALUE"
                    else
                        echo "WARNING: MAX_SOFTCLIP '$VALUE' is invalid!"
                    fi                    
                    ;;
                "NTHREADS")
                    if [[ $VALUE =~ ^[0-9]+$ ]]; then
                        NTHREADS="$VALUE"
                    else
                        echo "WARNING: NTHREADS '$VALUE' is invalid!"
                    fi
                    ;;
                *) echo "Warning: unknown parameter found: $KEY";;
            esac
        fi
    done < "$PARAM_FILE"
fi

echo "INFO: Setting DEMUX_DIRECTORY to $DEMUX_DIRECTORY"
echo "INFO: Setting BED_FILE to $BED_FILE"
echo "INFO: Setting REF_GENOME to $REF_GENOME"
echo "INFO: Setting ANNOTATION_FILE to $ANNOTATION_FILE"
echo "INFO: Setting MIN_LENGTH to $MIN_LENGTH"
echo "INFO: Setting MIN_BQ to $MIN_BQ"
echo "INFO: Setting MAX_N to $MAX_N"
echo "INFO: Setting MAX_INSERT_SIZE to $MAX_INSERT_SIZE"
echo "INFO: Setting MAX_SOFTCLIP to $MAX_SOFTCLIP"
echo "INFO: Setiting NTHREADS to $NTHREADS"
if $ISEQ; then echo "INFO: Setting ISEQ to true"; fi
if $BWA; then
    echo "INFO: Using BWA MEM aligner"
else
    echo "INFO: Using Bowtie2 aligner"
fi

# these need to be exported to be used in task script
export ISEQ BWA MIN_LENGTH MIN_BQ MAX_N MAX_INSERT_SIZE MAX_SOFTCLIP REF_GENOME DEMUX_DIRECTORY


#######################################################
#### STEP 1: Set up environment for pipeline tools ####
#######################################################

# find fastq files and dump to text files
find "$DEMUX_DIRECTORY" -xtype f -name "*_R1*.fastq.gz" | sort -V > fwd_reads.txt
find "$DEMUX_DIRECTORY" -xtype f -name "*_R2*.fastq.gz" | sort -V > rvs_reads.txt

# check for forward reads
if [ $(cat fwd_reads.txt | wc -l) == "0" ]; then
    echo "ERROR: No reads found in $DEMUX_DIRECTORY"
    echo "TIP: Maybe the fastq files are not in the *_R[12]*.fastq.gz format?"
    echo "TIP: If not, change lines 170 and 171 to the format they are in."
    exit 1
fi

# get number of samples
N_SAMPLES=$(cat fwd_reads.txt | wc -l )

# check reverse reads == forward reads
if [ "$N_SAMPLES" != $(cat rvs_reads.txt | wc -l) ]; then
    echo "ERROR: Forward and reverse reads do not have same number of files. Check fwd_reads.txt and rvs_reads.txt for details."
    exit 1
fi

echo "INFO: Found $N_SAMPLES samples in $DEMUX_DIRECTORY"

# test to make sure paired end reads line up correctly and contain reads!
# then for samples that have reads, get sample names
# for those that do not, put into bad_demux.txt

rm -f samples.txt bad_demux.txt

for i in `seq $N_SAMPLES`; do
    READ1=$(awk NR==$i fwd_reads.txt)
    READ2=$(awk NR==$i rvs_reads.txt)
    if [ "$READ1" != "${READ2/_R2/_R1}" ]; then
        echo "ERROR: Unmatched paired read files: $READ1 vs. $READ2";
        # exit 1
    else
        echo "$READ1" | sed "s|.*/||" | sed "s/_R1.*//" >> samples.txt
    fi
    if [ ! -s "$READ1" ] || [ ! -n "$(zcat $READ1 | head -n 1)" ] || [ ! -s "$READ2" ] || [ ! -n "$(zcat $READ2 | head -n 1)" ]; then
        echo "$READ1" | sed "s|.*/||" | sed "s/_R1.*//" >> bad_demux.txt
    fi
    
done

N_SAMPLES=$(cat samples.txt | wc -l)
if [ "$N_SAMPLES" -eq 0 ]; then
    echo "ERROR: No samples contained actual sequence data."
    exit 1
fi

bad_demux=$(cat bad_demux.txt 2> /dev/null | sed "s/$/-bad_demux/")
echo "WARNING: $(cat bad_demux.txt 2> /dev/null | wc -l) samples failed at demux."

# check sample names are in fact unique
if [ $(sort samples.txt | uniq | wc -l) != $(cat samples.txt |wc -l) ]; then
    echo "ERROR: Duplicate sample names detected."
    exit 1
fi

# clean/make all the required directories
rm -fr qc alignments bad_qc bad_alignment good_bams logs fastqc_preqc fastqc_postqc fastqc_aligned
mkdir -p qc alignments bad_qc bad_alignment good_bams logs fastqc_preqc fastqc_postqc fastqc_aligned

##################################################
#### STEP 2: QC and align reads (in parallel) ####
##################################################

if (( NTHREADS > 1 )); then
    seq $N_SAMPLE | parallel -j $NTHREADS ./AmpSeQC_task.sh 
else
    for TASK_ID in `seq N_SAMPLES`; do
        ./AmpSeQC_task.sh $TASK_ID
    done
fi

# run multiqc
echo "INFO: Generating MultiQC reports. Please wait..."
multiqc --force --dirs --outdir multiqc_report_full --export --interactive .
multiqc --force --dirs --outdir multiqc_report_preqc --export --interactive fastqc_preqc/
multiqc --force --dirs --outdir multiqc_report_postqc --export --interactive fastqc_postqc/
multiqc --force --dirs --outdir multiqc_report_aligned --export --interactive fastqc_aligned/

# get samples with bad qc
BAD_QC=$(find bad_qc/ -type f)
if [ $(echo "$BAD_QC" | wc -l) -gt 0 ]; then
    echo "WARNING: $BAD_QC samples failed QC."
    echo "WARNING: See bad_qc directory for list of samples."
fi

# get samples with bad alignments
BAD_ALIGNMENT=$(find bad_alignment/ -type f)
if [ $(echo "$BAD_ALIGNMENT" | wc -l) -gt 0 ]; then
    echo "WARNING: $BAD_ALIGNMENT samples failed alignment."
    echo "WARNING: See bad_alignment directory for list of samples and their bam files."
fi

# check to make sure there are good bam files
GOOD_BAMS=$(find good_bams/ -type f | wc -l)
if [ "$GOOD_BAMS" -gt 0 ]; then
    cat good_bams/* | sort -V > bams.txt
    echo "INFO: $GOOD_BAMS samples had good alignments."
else
    echo "ERROR: no samples had good alignments! Exiting..."
    exit 1
fi


#######################################################
#### STEP 3: Generate read count per amplicon/gene ####
#######################################################

echo "INFO: Generating read counts. Please wait..."

# add header with gff columns + sample names
echo -ne "amplicon\t" > "$BED_FILE"
cat bams.txt | sed "s|alignments/||" | sed "s|.proper_pairs.bam||" | tr '\n' '\t' >> "$BED_FILE"
FAILED_SAMPLES=$(echo -e "${bad_demux}${BAD_QC}${BAD_ALIGNMENT}")
echo "$FAILED_SAMPLES" | sed "s|\(bad_[a-z]\+\)/\(.*\)|\2-\1|" | tr '\n' '\t' >> "$BED_FILE"
sed -i "s/\t$/\n/" "$BED_FILE"

N_FAILED=$(echo "$FAILED_SAMPLES" | wc -l)
FAILED_COUNTS=$(printf '\t0%.0s' $(seq 0 $N_FAILED))

# get read counts per amplicon with htseq-count
htseq-count -r pos -s no -t "amplicon" -i "ID" -n "$NTHREADS" `cat bams.txt` "$ANNOTATION_FILE" | sed "s/$/$FAILED_COUNTS/" >> "$BED_FILE" && \
echo "INFO: Read counts per gene/amplicon for all samples can be found in $BED_FILE (amplicons as rows)"

# Check to make sure final output is produced
if [ ! -s "$BED_FILE" ]; then
    echo "ERROR: Could not create read count file $BED_FILE"
    exit 1
fi

echo "INFO: Completed successfully"
#### ALL DONE ####
