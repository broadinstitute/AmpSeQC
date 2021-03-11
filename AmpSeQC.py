#!/usr/bin/env python
"""AmpSeQC pipeline

An amplicon sequencing quality control pipeline
"""

import os
import sys
import gzip
import shutil
import subprocess

from multiprocessing import Pool

from argparse import ArgumentParser

def check_commands(no_fastqc=False, bowtie2=False):
    print("INFO: Verifying all tools are installed. Please wait...")
    cmds = [["trim_galore", "--version"], ["samtools", "--version"], ["htseq-count", "--version"]]
    
    if not no_fastqc:
        cmds.extend([["fastqc", "--version"], ["multiqc", "--version"]])

    if bowtie2:
        cmds.append(["bowtie2", "--version"])
    else:
        cmds.extend(["bwa", ["samclip", "--version"]])
    
    missing = False
    for cmd in cmds:
        try:
            subprocess.run(cmd, capture_output=True)
        except KeyboardInterrupt:
            sys.exit(1)
        except SystemExit:
            raise SystemExit
        except FileNotFoundError:
            missing = True
            print("ERROR: Cannot find %s. Check if program is installed and try again." % cmd)
    
    if missing:
        sys.exit(1)


def _fastq_reads(file):
    if "q.gz" in file:
        f = gzip.open(file, "rt")
    else:
        f = open(file)
    i = 0
    for line in f:
        i += 1
    f.close()
    return i / 4

def parse_fastq(files):
    """Parse fastq files to ensure everything is hunky dorey"""
    if len(files) % 2 != 0:
        print("ERROR: Not an even number of fastq files. Expects paired reads.", file=sys.stderr)
        sys.exit(1)
    read1 = []
    read2 = []
    for file in files:
        if "_R1_" in file:
            read1.append(file)
        elif "_R2_" in file:
            read2.append(file)
        elif ".1.f" in file:
            read1.append(file)
        elif ".2.f" in file:
            read2.append(file)
        else:
            print("ERROR: Cannot determine if %s is read1 or read2" % file, file=sys.stderr)
            sys.exit(1)
    
    if len(read1) != len(read2):
        print("ERROR: Uneven read1 and read2 file count!", file=sys.stderr)
        sys.exit(1)
    
    read1 = sorted(read1)
    read2 = sorted(read2)

    bad_demux = []
    samples = {}
    for i, file in enumerate(read1):
        paired = file.replace("_R1", "_R2")
        if paired == read2[i]:
            n_1 = _fastq_reads(file)
            n_2 = _fastq_reads(paired)
            if n_1 != n_2:
                print("ERROR: %s and %s do not have the same number of reads!" % (file, paired), file=sys.stderr)
                sys.exit(1)
            elif n_1 == 0:
                bad_demux.append(file.split("_R1")[0])
            else:
                sample = file.split("_R1")[0]
                if sample in samples:
                    print("ERROR: %s is not unique sample name!" % sample, file=sys.stderr)
                    sys.exit(1)
                samples[sample] = (file, read2[i])
            continue
        elif paired in read2:
            print("ERROR: Read1 and read2 files are out of order! Are there missing files?", file=sys.stderr)
            sys.exit(1)
        
        paired = file.replace(".1.fastq", ".2.fastq").replace(".1.fq", ".2.fq")
        if paired == read2[i]:
            n_1 = _fastq_reads(file)
            n_2 = _fastq_reads(paired)
            if n_1 != n_2:
                print("ERROR: %s and %s do not have the same number of reads!" % (file, paired), file=sys.stderr)
                sys.exit(1)
            elif n_1 == 0:
                bad_demux.append(file.split(".1.fastq")[0].split(".1.fq")[0])
            else:
                sample = file.split(".1.fastq")[0].split(".1.fq")[0]
                if sample in samples:
                    print("ERROR: %s is not unique sample name!" % sample, file=sys.stderr)
                    sys.exit(1)
                samples[sample] = (file, read2[i])
            continue
        elif paired in read2:
            print("ERROR: Read1 and read2 files are out of order! Are there missing files?", file=sys.stderr)
            sys.exit(1)
        
        print("ERROR: Cannot determined paired file for %s" % file, file=sys.stderr)
        sys.exit(1)
    
    return samples, bad_demux


def run_fastqc(read1, read2, out):
    """Run FastQC on a set of paired reads"""
    cmd = "fastqc -o %s --noextract -f fastq %s %s" % (out, read1, read2)
    try:
        subprocess.run(cmd, check=True)
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except SystemExit:
        raise SystemExit
    except Exception as e:
        print(e, file=sys.stderr)
        return False
    
    return True


def qc_sample(sample, fwd, rvs, ref="reference.fasta", two_color=False, min_length=70, min_bq=20, max_N=1, no_fastqc=False):
    """Run a sample through QC and alignment"""

    if not no_fastqc and not run_fastqc(fwd, rvs, "fastqc_preqc"):
        print("WARNING: Could not run FastQC on pre-QC data!", file=sys.stderr)

    if two_color:
        cmd = "trim_galore --paired --length %d --2colour %d --max_n %d --dont_gzip -o qc %s %s" % (min_length, min_bq, max_N, fwd, rvs)
    else:
        cmd = "trim_galore --paired --length %d -q %d --max_n %d --dont_gzip -o qc %s %s" % (min_length, min_bq, max_N, fwd, rvs)
    
    output = subprocess.run(cmd, check=True, stdout="logs/%s.trim_galore.out" % sample, stderr="logs/%s.trim_galore.err" % sample)
    if output.returncode:
        print("ERROR: Trim-Galore failed on %s" % sample)
        return

    read1 = "qc/$%s_R1.fastq.gz" % sample
    read2 = "qc/$%s_R2.fastq.gz" % sample

    try:
        # TODO: make this more pythonic, less shell
        subprocess.run("gzip -c qc/%s*_R1_*.fq > %s" % (sample, read1), shell=True, check=True)
        subprocess.run("gzip -c qc/%s*_R2_*.fq > %s" % (sample, read2), shell=True, check=True)
        subprocess.run("rm -f qc/%s*.fq" % (sample), shell=True, check=True)
    except:
        print("ERROR: Failed to process QC-ed files for %s" % sample, file=sys.stderr)
        return

    if os.path.isfile(read1) and os.path.isfile(read2):
        n_1 = _fastq_reads(read1)
        n_2 = _fastq_reads(read2)
        if n_1 and n_1 == n_2:
            print("INFO: QC complete for %s" % sample, file=sys.stderr)
            if not no_fastqc and not run_fastqc(read1, read2, "fastqc_postqc"):
                print("WARNING: Could not run FastQC on post-QC data!", file=sys.stderr)
            return True
    
    print("WARNING: QC of %s failed" % sample, file=sys.stderr)
    return False


    

def align_sample(sample, ref="reference.fasta", max_insert_size=500, soft_clip=5, bowtie2=False, no_fastqc=False):
    read1 = "qc/$%s_R1.fastq.gz" % sample
    read2 = "qc/$%s_R2.fastq.gz" % sample

    w = open("logs/%s.alignment.err", "w")
    if bowtie2:
        pass
    else:
        subprocess.run("bwa mem -I '200,100,%d,50' %s %s %s" % (max_insert_size, ref, read1, read2), shell=True, check=True, stderr=w, stdout="alignments/%s.sam" % sample)
        if os.path.isfile("alignments/%s.sam"):
            output = subprocess.run("samtools flagstat alignments/%s.sam | head -n 1 | awk '{ print $1 }'" % sample, shell=True, check=True, stdout=subprocess.PIPE, stderr=w)
            if int(output.stdout) > 0:
                subprocess.run("samclip --ref {ref} --max {soft_clip} alignments/{sample}.sam \
                    | samtools fixmate - - | samtools view -bf 3 | samtools sort -T alignments/{sample} > alignments/{sample}.proper_pairs.bam && \
                        samtools index alignments/{sample}.proper_pairs.bam && samtools sort -T alignments/{sample} alignments/{sample}.sam > alignments/{sample}.bam".format(ref=ref, soft_clip=soft_clip, sample=sample), check=True, shell=True, stderr=w)
                os.remove("alignments/%s.sam" % sample)
                if os.path.isfile("alignments/%s.proper_pairs.bam" % sample):
                    output = subprocess.run("samtools flagstat alignments/%s.proper_pairs.bam | head -n 1 | awk '{ print $1 }'" % sample, shell=True, check=True, stdout=subprocess.PIPE, stderr=w)
                    if int(output.stdout) > 0:
                        print("INFO: Alignment of %s finished with %s reads." % (sample, output.stdout), file=sys.stderr)
                    else:
                        print("WARNING: Alignment of %s finished without any good aligned reads." % sample, file=sys.stderr)
                else:
                    print("ERROR: Did not generate good alignment for %s" % sample, file=sys.stderr)
            else:
                print("ERROR: Alignment has zero reads for %s!" % sample, file=sys.stderr)
        else:
            print("ERROR: No alignment file for %s!" % sample, file=sys.stderr)
            return False

    if not no_fastqc:
        try:
            subprocess.run("fastqc -o fastqc_aligned --noextract -f bam alignments/%s.proper_pairs.bam" % sample, shell=True, check=True, stderr=w)
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except SystemExit:
            raise SystemExit
        except Exception as e:
            print(e, file=sys.stderr)
            print("WARNING: Could not run FastQC on aligned data!", file=sys.stderr)
    w.close()

    return True


def _process_sample(args):
    return process_sample(*args)

def process_sample(sample, fwd, rvs, ref="reference.fasta", two_color=False, min_length=70, min_bq=20, max_N=1, max_insert_size=500, soft_clip=5, bowtie2=False, no_fastqc=False):
    """QC and align sample"""
    if not qc_sample(sample, fwd, rvs, ref=ref, two_color=two_color, min_length=min_length, min_bq=min_bq, max_N=max_N, no_fastqc=no_fastqc):
        return "bad_qc"
    if not align_sample(sample, ref=ref, max_insert_size=max_insert_size, soft_clip=soft_clip, bowtie2=bowtie2, no_fastqc=no_fastqc):
        return "bad_alignment"
    return "good_alignment"


def filter_count_data(count_data, min_sample_count=0, min_amplicon_count = 0):
    if not (min_amplicon_count or min_amplicon_count):
        return count_data
    
    amplicons = {}
    for sample in count_data:
        if min_sample_count > 0:
            total_count = sum([count_data[sample][amplicon] for amplicon in count_data[sample]])
            if total_count < min_sample_count:
                continue
        for amplicon in count_data[sample]:
            if amplicon not in amplicons:
                amplicons[amplicon] = {}
            amplicons[amplicon][sample] = count_data[sample][amplicon]
    
    if min_amplicon_count > 0:
        filtered_amplicons = {}
        for amplicon in amplicons:
            total_count = sum([amplicons[amplicon][sample] for sample in amplicons[amplicon]])
            if total_count >= min_amplicon_count:
                filtered_amplicons[amplicon] = amplicons[amplicon]
        amplicons = filtered_amplicons
    
    filtered_data = {}
    for amplicon in amplicons:
        for sample in amplicons[amplicon]:
            if sample not in filtered_data:
                filtered_data[sample] = {}
            filtered_data[sample][amplicon] = amplicons[amplicon][sample]
    return filtered_data



def run_multiqc():
    """Run MultiQC"""
    cmds = ["multiqc --force --dirs --outdir multiqc_report_full --export --interactive . &> logs/multiqc.log", \
            "multiqc --force --dirs --outdir multiqc_report_preqc --export --interactive fastqc_preqc/ &>> logs/multiqc.log", \
            "multiqc --force --dirs --outdir multiqc_report_postqc --export --interactive fastqc_postqc/ &>> logs/multiqc.log", \
            "multiqc --force --dirs --outdir multiqc_report_aligned --export --interactive fastqc_aligned/ &>> logs/multiqc.log"]
    
    for cmd in cmds:
        try:
            subprocess.run(cmd, check=True, shell=True)
        except KeyboardInterrupt:
            break
        except SystemExit:
            break
        except Exception as e:
            print("WARNING: Exception occurred running MultiQC...\n%s" % e, file=sys.stderr)
            break
    else:
        return True
    

def main():
    parser = ArgumentParser(description="Amplicon sequencing quality control pipeline", usage="%(prog)s [options] -c [read_counts] -r [ref.fasta] -a [annot.gff] fastq [fastq ...]")
    parser.add_argument("fastq", nargs="+", help="Fastq file(s) to analyze (expects paired-end reads in separate files)")
    parser.add_argument("-c", "--counts", default="read_counts.tsv", help="Read count tsv file (default: read_counts.tsv)")
    parser.add_argument("-r", "--ref", default="reference.fasta", help="Indexed reference fasta file to align to (default: reference.fasta)")
    parser.add_argument("-a", "--annot", default="amplicons.gff", help="Amplicon/gene gff3 file to generate read counts of (default: amplicons.gff)")
    parser.add_argument("--2color", dest="two_color", action="store_true", default=False, help="Run with 2 color chemistry (e.g. iSeq) QC parameters (default: False")
    parser.add_argument("-l", "--min_length", type=int, default=70, help="Minimum read length (in bp) to retain after trimming (default: 70)")
    parser.add_argument("-q", "--min_bq", type=int, default=20, help="Minimum base quality (PHRED score) to retain after trimming (default: 20)")
    parser.add_argument("-N", "--max_N", type=int, default=1, help="Maximum number of 'N' bases allowed in read (default: 1)")
    parser.add_argument("-I", "--max_insert_size", type=int, default=500, help="Maximum insert size (in bp) expected for aligner (default: 500)")
    parser.add_argument("-S", "--soft_clip", type=int, default=5, help="Maximum soft clipping (in bp) allowed for BWA-Mem (default: 5)")
    parser.add_argument("--bowtie2", action="store_true", default=False, help="Align with Bowtie2 instead of BWA-Mem (default: False)")
    parser.add_argument("--min_amplicon_count", type=int, default=0, help="Minimum total read count to retain an amplicon (default: 0)")
    parser.add_argument("--min_sample_count", type=int, default=0, help="Minimum total read count to retain a sample (default: 0)")
    parser.add_argument("--no_fastqc", action="store_true", default=False, help="Do not run FastQC or MultiQC")
    parser.add_argument("-p", "--procs", type=int, default=1, help="Number of processors to use (default: 1)")
    args = parser.parse_args()

    check_commands(no_fastqc=args.no_fastqc, bowtie2=args.bowtie2)

    samples, bad_demux = parse_fastq(args.fastq)

    for folder in ["qc", "alignments", "bad_qc", "bad_alignment", "good_bams", "logs", "fastqc_preqc", "fastqc_postqc", "fastqc_aligned"]:
        try:
            shutil.rmtree(folder)
        except FileNotFoundError:
            pass
        os.mkdir(folder)
    
    print("QC-ing %d samples. Please wait..." % len(samples), file=sys.stderr)
    
    if args.procs > 1:
        print("INFO: Using %d processors" % args.procs, file=sys.stderr)
        with Pool(processes=args.procs) as pool:
            pool_args = [(sample, samples[sample][0], samples[sample][1], args.ref, args.two_color, args.min_length, args.min_bq, args.max_N, args.max_insert_size, args.soft_clip, args.bowtie2, args.no_fastqc) for sample in samples]
            result = pool.apply_async(_process_sample, pool_args)
            qc_results = result.get()
    else:
        qc_results = {sample: process_sample(sample, samples[sample][0], samples[sample][1], ref=args.ref, two_color=args.two_color, min_length=args.min_length, min_bq=args.min_bq, max_N=args.max_N, max_insert_size=args.max_insert_size, soft_clip=args.soft_clip, bowtie2=args.bowtie2, no_fastqc=args.no_fastqc) for sample in samples}
    
    bad_qc = sorted([sample for sample in qc_results if qc_results[sample] == "bad_qc"])
    bad_alignment = sorted([sample for sample in qc_results if qc_results[sample] == "bad_alignment"])
    good_alignment = sorted([sample for sample in qc_results if qc_results[sample] == "good_alignment"])
    bams = " ".join(["alignments/%s.proper_pairs.bam" % sample for sample in good_alignment])

    print("""INFO: Good data: %d
    INFO: Bad demux (no sequencing data): %d
    INFO: Bad QC (failed Trim-Galore): %d
    INFO: Bad alignment (no properly aligned reads): %d""" % (len(good_alignment), len(bad_demux), len(bad_qc), len(bad_alignment)), file=sys.stderr)
    
    print("INFO: Running htseq-count to generate read counts per amplicon per sample. Please wait...", file=sys.stderr)
    count_data = {sample: {} for sample in good_alignment}
    amplicons = []
    output = subprocess.run("htseq-count -r pos -s no -t 'amplicon' -i 'ID' -n %d %s %s" % (args.procs, bams, args.annot), capture_output=True, check=True)
    for line in output.stdout.split("\n"):
        line = line.split()
        amplicon = line[0]
        amplicons.append(amplicon)
        counts = [int(x) for x in line[1:]]
        for i, sample in enumerate(good_alignment):
            count_data[sample][amplicon] = counts[i]
    
    if args.min_sample_count > 0 or args.min_amplicon_count > 0:
        count_data = filter_count_data(count_data, min_sample_count=args.min_sample_count, min_amplicon_count=args.min_amplicon_count)
    
    if args.min_amplicon_count > 0:
        filtered_amplicons = set()
        for sample in count_data:
            filtered_amplicons.update(list(count_data[sample].keys()))
        print("INFO: %d / %d amplicons remain after read count filtering" %(len(filtered_amplicons), len(amplicons)), file=sys.stderr)
        amplicons = sorted(list(filtered_amplicons))
    
    if args.min_sample_count:
        low_alignment = [sample for sample in good_alignment if sample not in count_data]
        print("INFO: Filtered %d / %d aligned samples due to low read counts (< %d total)" %(len(low_alignment), len(good_alignment), args.min_sample_count), file=sys.stderr)
    

    print("INFO: Writing read counts to %s" % args.counts, file=sys.stderr)
    zero_line = "\t0" * len(amplicons)
    with open(args.counts, "w") as w:
        w.write("sample\t" + "\t".join(amplicons)+"\n")
        for sample in sorted(count_data):
            w.write(sample+"\t"+"\t".join(["%d" % count_data[sample][amplicon] for amplicon in amplicons])+"\n")
        if args.min_sample_count <= 0:
            for sample in bad_demux:
                w.write("%s.bad_demux%s\n" %(sample, zero_line))
            for sample in bad_qc:
                w.write("%s.bad_qc%s\n" %(sample, zero_line))
            for sample in bad_alignment:
                w.write("%s.bad_alignment%s\n" %(sample, zero_line))

    print("INFO: Running MultiQC. Please wait...", file=sys.stderr)
    if not args.no_fastqc and not run_multiqc():
        print("ERROR: Error running MultiQC, exiting...", file=sys.stderr)
        sys.exit(1)
    
    print("INFO: All steps completed successfully. Goodbye!", file=sys.stderr)


if __name__ == "__main__":
    main()