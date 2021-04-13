#!/usr/bin/env python
import argparse
import sys
import os
import subprocess
import tempfile

from Bio import SeqIO, AlignIO

# full gene targets for the amplicons in the gt-seq panel (except vivax)
AMPLICON_DATABASE="/gsap/garage-protistvector/ampseq_data/AmpSeQC/amplicon_genes.fasta"
AMP_TO_GENE="/gsap/garage-protistvector/ampseq_data/AmpSeQC/amplicon_to_gene_id.txt"

# parse amplicon name to gene id
def parse_amp_to_gene(file=AMP_TO_GENE):
    gene2amp = {}
    with open(file) as f:
        for line in f:
            line=line.strip().split("\t")
            amp = line[0]
            gene = line[1]
            if gene not in gene2amp:
                gene2amp[gene] = []
            gene2amp[gene].append(amp)
    return gene2amp

# parse amplicon database
def parse_amp_db(fasta_file=AMPLICON_DATABASE, amp_to_gene_file=AMP_TO_GENE):
    gene2amp = parse_amp_to_gene(amp_to_gene_file)
    amplicons = {}
    for seq in SeqIO.parse(fasta_file, "fasta"):
        gene = seq.id.split("::")[0]
        for amp in gene2amp[gene]:
            amplicons[amp] = seq
    return amplicons


# parse asv to amplicon table
def parse_asv_table(file, min_reads=50, min_samples=2, max_dist=30):
    bins = {}
    with open(file) as f:
        f.readline()
        for line in f:
            line = line.strip().split("\t")
            nreads = int(line[1])
            if nreads < min_reads:
                continue
            nsamples = int(line[2])
            if nsamples < min_samples:
                continue
            dist = min(int(line[5]), int(line[7]))
            if dist > max_dist:
                continue
            ASV = line[0]
            amplicon = line[4]
            if amplicon not in bins:
                bins[amplicon] = []
            bins[amplicon].append(ASV)
    return bins


# parse ASV fasta file
def get_asv_seqs(file):
    return {seq.id: seq for seq in SeqIO.parse(file, "fasta")}


def write_bins(asvs, bins, amplicons, outdir="ASVs"):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    for amplicon in bins:
        if amplicon not in amplicons:
            print(f"WARNING: {amplicon} target not found in amplicon sequence database", file=sys.stderr)
            continue
        with open(os.path.join(outdir, f"{amplicon}.fasta"), "w") as w:
            SeqIO.write(amplicons[amplicon], w, "fasta")
            SeqIO.write([asvs[asv] for asv in bins[amplicon]], w, "fasta")


# run muscle for each amplicon
def run_muscle(bins, outdir="ASVs"):
    for amplicon in bins:
        fasta = os.path.join(outdir, f"{amplicon}.fasta")
        if not os.path.isfile(fasta):
            print(f"ERROR: Could not find {fasta}", file=sys.stderr)
            continue
        msa = os.path.join(outdir, f"{amplicon}.msa")
        subprocess.run(["muscle", "-in", fasta, "-out", msa], capture_output=True)


# get coords of amplicons (basically where there aren't gaps across all sequences)
def _find_asv_coords(alignment):
    gaps = {i:0 for i in range(alignment.get_alignment_length())}
    for seq in alignment[1:]:
        for i in gaps:
            if seq[i] == '-':
                gaps[i] += 1
    n_asvs = len(alignment) - 1
    non_gaps = [i for i in gaps if gaps[i] < n_asvs]
    return min(non_gaps), max(non_gaps)


# parse muscle alignment
def parse_alignment(alignment):
    aln = AlignIO.read(alignment, "fasta")
    aln.sort(key = lambda record: (record.id[:5] != "PF3D7", record.id))
    anchor = aln[0]
    if anchor.id[:5] != "PF3D7":
        print(f"ERROR: No anchor gene for {alignment}", file=sys.stderr)

    start, end = _find_asv_coords(aln)
    if len(aln[0].seq.lstrip("-")) != aln.get_alignment_length():
        print(f"WARNING: {alignment} extends beyond 5' end of reference gene. ASVs may include non-genic sequence.", file=sys.stderr)
    elif len(aln[0].seq.rstrip("-")) != aln.get_alignment_length():
        print(f"WARNING: {alignment} extends beyond 3' end of reference gene. ASVs may include non-genic sequence.", file=sys.stderr)

    asv_to_cigar = {}
    for seq in aln[1:]:
        pos = start + 1
        cigar = ""
        for i in range(start, end):
            if seq[i] != anchor[i]:
                if anchor[i] == "-":
                    if i == start or anchor[i-1] != "-":
                        cigar += f"{pos}I="
                    cigar += seq[i]
                    continue
                elif seq[i] == "-":
                    if i == start or seq[i-1] != "-":
                        cigar += f"{pos}D="
                    cigar += f"{anchor[i]}"
                else:
                    cigar += f"{pos}{seq[i]}"
            pos += 1

        if not cigar:
            cigar = f"M"
        asv_to_cigar[seq.id] = cigar
    return asv_to_cigar


# get variants per amplicon per position
def parse_alignments(bins, outdir="ASVs"):
    cigars = {}
    for amplicon in bins:
        msa = os.path.join(outdir, f"{amplicon}.msa")
        if not os.path.isfile(msa):
            print(f"ERROR: Could not find {msa}", file=sys.stderr)
            continue
        cigars[amplicon] = parse_alignment(msa)
    
    return cigars

# write table of asv -> amplicon/cigar
def write_cigar_strings(cigars, out="CIGARs.tsv"):
    with open(out, 'w') as w:
        w.write("ASV\tAmplicon\tCIGAR\n")
        for amplicon in sorted(cigars):
            for ASV in sorted(cigars[amplicon], key = lambda x: int(x[1:])):
                w.write(f"{ASV}\t{amplicon}\t{cigars[amplicon][ASV]}\n")


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", required=True, help="Fasta file of ASV sequences from DADA2 pipeline")
parser.add_argument("-t", "--table", required=True, help="ASV table from DADA2 pipeline")
parser.add_argument("-a", "--alignments", required=True, help="Directory to store ASV alignment files")
parser.add_argument("-o", "--out", required=True, help="Output file for ASV -> CIGAR string table")
parser.add_argument("--min_reads", type=int, default=50, help="Minimum total reads to include ASV (default: 50)")
parser.add_argument("--min_samples", type=int, default=50, help="Minimum samples to include ASV (default: 2)")
parser.add_argument("--max_dist", type=int, default=30, help="Maximum edit distance to include ASV (default: 30)")
parser.add_argument("--amp_db", default=AMPLICON_DATABASE, help=f"Amplicon sequence database (default: {AMPLICON_DATABASE})")
parser.add_argument("--amp_to_gene", default=AMP_TO_GENE, help=f"Amplicon -> gene table (default: {AMP_TO_GENE})")
args = parser.parse_args()

print(f"INFO: Loading {args.amp_db} and {args.amp_to_gene}")
amplicons = parse_amp_db(args.amp_db, args.amp_to_gene)
if not amplicons:
    print(f"ERROR: No amplicons in {args.amp_db}", file=sys.stderr)
    sys.exit(1)

print(f"INFO: Loading {args.fasta}")
asvs = get_asv_seqs(args.fasta)
if not asvs:
    print(f"ERROR: No ASV sequences in {args.fasta}", file=sys.stderr)
    sys.exit(1)

print(f"INFO: Parsing {args.table} with total reads >= {args.min_reads}, samples >= {args.min_samples}, dist <= {args.max_dist}")
bins = parse_asv_table(args.table, min_reads=args.min_reads, min_samples=args.min_samples, max_dist=args.max_dist)
if not bins:
    print(f"ERROR: No useable data in {args.table}", file=sys.stderr)
    sys.exit(1)

outdir = args.alignments
print(f"INFO: Writing amplicon fasta files to {outdir}")
if not os.path.isdir(outdir):
    os.mkdir(outdir)
write_bins(asvs, bins, amplicons, outdir=outdir)

print("INFO: Running MUSCLE aligner on amplicon fasta files. Please wait...", file=sys.stderr)
run_muscle(bins, outdir=outdir)

print("INFO: Parsing alignments to CIGAR strings")
cigars = parse_alignments(bins, outdir=outdir)
if not cigars:
    print("ERROR: could not determine CIGAR strings", file=sys.stderr)
    sys.exit(1)

write_cigar_strings(cigars, args.out)
print(f"INFO: Wrote ASV->CIGAR to {args.out}")