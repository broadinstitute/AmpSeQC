#!/bin/bash

TASK_ID="$@"

SAMPLE="$(awk NR==$TASK_ID samples.txt)"
FWD="$(awk NR==$TASK_ID fwd_reads.txt)"
RVS="$(awk NR==$TASK_ID rvs_reads.txt)"

rm -f logs/${SAMPLE}.*

exec 1<&-
exec 2<&-
exec 1<>logs/${SAMPLE}.out
exec 2<>logs/${SAMPLE}.err

source activate $CONDA_ENV

READ1="qc/${SAMPLE}_R1.fastq.gz"
READ2="qc/${SAMPLE}_R2.fastq.gz"

if [ -s bad_demux.txt ]&& grep -q "^${SAMPLE}$" bad_demux.txt; then
        echo "WARNING: Sample failed demux. Skipping QC and alignment."
        exit 0
elif [ -s "$READ1" ] && [ -s "$READ2" ]; then
    echo "WARNING: Already QC-ed $SAMPLE. Skipping QC step."
else
    echo "INFO: Running FastQC on pre-QC reads. please wait..."
    fastqc -o fastqc_preqc --noextract -f fastq "$FWD" "$RVS"

    echo "INFO: QC-ing ${SAMPLE}. Please wait..."
    if $ISEQ; then
        echo "INFO: Trimming with iSeq quality parameters!"
        trim_galore --paired --length "$MIN_LENGTH" --2colour "$MIN_BQ" --max_n "$MAX_N" --dont_gzip -o qc "$FWD" "$RVS"
    else
        trim_galore --paired --length "$MIN_LENGTH" -q "$MIN_BQ" --max_n "$MAX_N" --dont_gzip -o qc "$FWD" "$RVS"
    fi
    if [ $? -eq 0 ]; then
        cat qc/"${SAMPLE}"*_R1_*.fq > "$READ1" && \
        cat qc/"${SAMPLE}"*_R2_*.fq > "$READ2" && \
        rm -f qc/"${SAMPLE}"*_R1_*.fq qc/"${SAMPLE}"*_R2_*.fq
        echo "INFO: QC of ${SAMPLE} complete."
    else
        echo "WARNING: TrimGalore failed!"
    fi
fi

rm -f "bad_qc/${SAMPLE}" "good_bams/${SAMPLE}" "bad_alignment/${SAMPLE}"

if [ ! -s "$READ1" ] || [ ! -s "$READ2" ]; then
    echo "WARNING: $SAMPLE did not have QC-ed reads to align to."
    touch "bad_qc/${SAMPLE}"
    exit 1
fi

echo "INFO: Running FastQC on post-QC reads. Please wait..."
fastqc -o fastqc_postqc --noextract -f fastq "$READ1" "$READ2"
if $BWA; then
    echo "INFO: Aligning $SAMPLE with BWA-MEM to $PARASITE_GENOME. Please wait..."
    bwa mem -I "200,100,$MAX_INSERT_SIZE,50" "$PARASITE_GENOME" "$READ1" "$READ2" > "alignments/${SAMPLE}.sam"
    if [ ! -s "alignments/${SAMPLE}.sam" ] || [ $(samtools flagstat "alignments/${SAMPLE}.sam" | head -n 1 | awk '{ print $1 }') == "0" ]; then
            echo "WARNING: No alignment data produced for $SAMPLE!"
            echo "alignments/${SAMPLE}.sam" > "bad_alignment/$SAMPLE"
            exit 1
    fi
    samclip --ref "$PARASITE_GENOME" --max "$MAX_SOFTCLIP" "alignments/${SAMPLE}.sam" | samtools fixmate - - | samtools view -bf 3 | samtools sort -T "alignments/${SAMPLE}" > "alignments/${SAMPLE}.proper_pairs.bam" && \
    samtools index "alignments/${SAMPLE}.proper_pairs.bam" && \
    samtools sort -T "alignments/${SAMPLE}" "alignments/${SAMPLE}.sam" > "alignments/${SAMPLE}.bam" && \
    rm "alignments/${SAMPLE}.sam"
else
    echo "INFO: Aligning $SAMPLE with Bowtie2 to $PARASITE_GENOME. Please wait..."
    bowtie2 -X "$MAX_INSERT_SIZE" --very-sensitive -x "$PARASITE_GENOME" -1 "$READ1" -2 "$READ2" | samtools sort -T "alignments/${SAMPLE}" > "alignments/${SAMPLE}.bam" && \
    samtools index "alignments/${SAMPLE}.bam"
    if [ ! -s "alignments/${SAMPLE}.bam" ] || [ $(samtools flagstat "alignments/${SAMPLE}.bam" | head -n 1 | awk '{ print $1 }') == "0" ]; then
        echo "WARNING: No alignment data produced for $SAMPLE!"
        echo "alignments/${SAMPLE}.bam" > "bad_alignment/$SAMPLE"
        exit 1
    fi
    samtools view -bf 3 "alignments/${SAMPLE}.bam" > "alignments/${SAMPLE}.proper_pairs.bam" && \
    samtools index "alignments/${SAMPLE}.proper_pairs.bam"
fi

if [ ! -s "alignments/${SAMPLE}.proper_pairs.bam" ] || [ $(samtools flagstat "alignments/${SAMPLE}.proper_pairs.bam" | head -n 1 | awk '{ print $1 }') == "0" ]; then
    echo "WARNING: No good alignment data produced for $SAMPLE!"
    echo "alignments/${SAMPLE}.bam" > "bad_alignment/$SAMPLE"
    exit 1
else
    echo "INFO: ${SAMPLE} finished alignment with at least some good data."
    echo "alignments/${SAMPLE}.proper_pairs.bam" > "good_bams/${SAMPLE}"
    echo "INFO: Running FastQC on aligned reads. Please wait..."
    fastqc -o fastqc_aligned --noextract -f bam "alignments/${SAMPLE}.proper_pairs.bam"
fi

exit 0