#!/usr/bin/env python
import argparse
import sys
import os
import subprocess

from Bio import SeqIO, AlignIO, Seq

# full gene targets for the amplicons in the gt-seq panel (except vivax)
AMPLICON_DATABASE="/gsap/garage-protistvector/ampseq_data/AmpSeQC/amplicons.fasta"
AMPLICON_MASK_INFO="/gsap/garage-protistvector/ampseq_data/AmpSeQC/amplicons.mask"

verbose = False

# parse amplicon dust mask info
def parse_dustmasker(mask_info=AMPLICON_MASK_INFO):
    if not mask_info:
        return
    mask = {}
    with open(mask_info) as f:
        for line in f:
            line = line.strip().split("\t")
            gene = line[0].split(":")[0][1:]
            if gene not in mask:
                mask[gene] = set()
            start = int(line[1])
            end = int(line[2])+1
            mask[gene].update(list(range(start, end)))
    return mask


# parse amplicon database
def parse_amp_db(fasta_file=AMPLICON_DATABASE):
    amplicons = {}
    for seq in SeqIO.parse(fasta_file, "fasta"):
        amplicons[seq.id] = seq
    return amplicons


# parse asv to amplicon table
def parse_asv_table(file, min_reads=0, min_samples=0, max_dist=-1):
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
            dist = min(int(line[5])+int(line[6]), int(line[8])+int(line[9]))
            if max_dist >= 0 and dist > max_dist:
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

# write amplicon fasta files
def wrte_amplicon_fastas(asvs, bins, amplicons, outdir="ASVs"):
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


# get coords of amplicons (where there aren't gaps across all ASVs)
def _find_asv_coords(alignment):
    gaps = {i:0 for i in range(alignment.get_alignment_length())}
    for seq in alignment[1:]:
        for i in gaps:
            if seq[i] == '-':
                gaps[i] += 1
    n_asvs = len(alignment) - 1
    non_gaps = [i for i in gaps if gaps[i] < n_asvs]
    return min(non_gaps), max(non_gaps)


# get coords of homopolymer runs
def _get_homopolymer_runs(seq, min_length=5):
    runs = set()
    prev = ""
    run = 0
    start = None
    last_non_gap = None
    for i in range(len(seq)):
        if seq[i] == "-":
            continue
        if seq[i] == prev:
            if not start:
                start = last_non_gap
            run += 1
        else:
            if run >= min_length:
                runs.update(list(range(start, i)))
            run = 0
            start = None
        prev = seq[i]
        last_non_gap = i
    
    return runs


# parse muscle alignment
def parse_alignment(alignment, mask={}, min_homopolymer_length=5, amplicon=None):
    aln = AlignIO.read(alignment, "fasta")
    aln.sort(key = lambda record: (record.id != amplicon, record.id))
    anchor = aln[0]
    if anchor.id != amplicon:
        print(f"ERROR: No anchor gene for {alignment}", file=sys.stderr)

    start, end = _find_asv_coords(aln)

    if min_homopolymer_length > 1:
        homopolymer_runs = _get_homopolymer_runs(aln[0], min_length=min_homopolymer_length)

    if len(aln[0].seq.lstrip("-")) != aln.get_alignment_length():
        print(f"WARNING: {os.path.basename(alignment)} extends beyond 5' end of reference gene.", file=sys.stderr)
    elif len(aln[0].seq.rstrip("-")) != aln.get_alignment_length():
        print(f"WARNING: {os.path.basename(alignment)} extends beyond 3' end of reference gene.", file=sys.stderr)

    masked = mask.get(aln[0].id, None)

    asv_to_cigar = {}
    for seq in aln[1:]:
        pos = start + 1
        cigar = ""
        indel = False
        for i in range(start, end):
            if masked and (pos-1) in masked:
                pass # masked position in gene. mask info is 0-based :(
            elif min_homopolymer_length > 1 and i in homopolymer_runs:
                if i and i-1 not in homopolymer_runs and seq.id == aln[1].id:
                    if verbose:
                        print(f"INFO: Skipping homopolymer run (poly-{anchor[i]}) beginning at position {pos} in {os.path.basename(alignment)}", file=sys.stderr)
            elif seq[i] != anchor[i]:
                if anchor[i] == "-":
                    if not indel:
                        indel = True
                        cigar += f"{pos}I="
                        if i:
                            for j in range(1,len(anchor)-i):
                                if anchor[i-j] != "-":
                                    cigar += anchor[i-j]
                                    break
                    cigar += seq[i]
                    continue
                if seq[i] == "-":
                    if not indel:
                        indel = True
                        cigar += f"{pos}D="
                    cigar += f"{anchor[i]}"
                else:
                    cigar += f"{pos}{seq[i]}"
                    indel = False
            else:
                indel = False
            pos += 1

        if not cigar:
            cigar = f"M={start+1}-{pos}"
        asv_to_cigar[seq.id] = cigar
    return asv_to_cigar


# get variants per amplicon per position
def parse_alignments(bins, mask={}, min_homopolymer_length=5, outdir="ASVs"):
    cigars = {}
    for amplicon in bins:
        msa = os.path.join(outdir, f"{amplicon}.msa")
        if not os.path.isfile(msa):
            print(f"ERROR: Could not find {msa}", file=sys.stderr)
            continue
        cigars[amplicon] = parse_alignment(msa, mask=mask, min_homopolymer_length=min_homopolymer_length, amplicon=amplicon)
    
    return cigars

# write table of asv -> amplicon/cigar
def write_cigar_strings(cigars, out="CIGARs.tsv"):
    with open(out, 'w') as w:
        w.write("ASV\tAmplicon\tCIGAR\n")
        for amplicon in sorted(cigars):
            for ASV in sorted(cigars[amplicon], key = lambda x: int(x[1:])):
                w.write(f"{ASV}\t{amplicon}\t{cigars[amplicon][ASV]}\n")


parser = argparse.ArgumentParser(usage="%(prog)s [options] fasta table alignments out",
                                 description="Convert ASVs from DADA2 pipeline to pseudo-CIGAR strings.",
                                 epilog="Contact tim.straub@broadinstitute.org for details.")
parser.add_argument("fasta", help="Fasta file of ASV sequences from DADA2 pipeline")
parser.add_argument("table", help="ASV table from DADA2 pipeline")
parser.add_argument("alignments", help="Directory to store ASV alignment files")
parser.add_argument( "out", help="Output file for ASV -> CIGAR string table")
parser.add_argument("-p", "--polyN", type=int, default=5, help="Mask homopolymer runs length >= polyN (default: 5; disabled < 2)")
parser.add_argument("--min_reads", type=int, default=0, help="Minimum total reads to include ASV (default: 0)")
parser.add_argument("--min_samples", type=int, default=0, help="Minimum samples to include ASV (default: 0)")
parser.add_argument("--max_dist", type=int, default=-1, help="Maximum edit distance to include ASV (default: -1, disabled)")
parser.add_argument("--amp_db", default=AMPLICON_DATABASE, help=f"Amplicon sequence fasta file (default: {AMPLICON_DATABASE})")
parser.add_argument("--amp_mask", default=AMPLICON_MASK_INFO, help=f"Amplicon low complexity mask info (default: {AMPLICON_MASK_INFO}, enter 'None' to disable)")
parser.add_argument("-v", "--verbose", default=False, action='store_true', help="Increase verobsity")
args = parser.parse_args()

if args.verbose:
    verbose = True

print(f"INFO: Loading {args.amp_db}", file=sys.stderr)
amplicons = parse_amp_db(args.amp_db)
if not amplicons:
    print(f"ERROR: No amplicons in {args.amp_db}", file=sys.stderr)
    sys.exit(1)


if args.amp_mask in ["None", 'none', 'NONE']:
    mask = {}
else:
    print(f"INFO: Loading {args.amp_mask}", file=sys.stderr)
    mask = parse_dustmasker(args.amp_mask)

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
wrte_amplicon_fastas(asvs, bins, amplicons, outdir=outdir)

print("INFO: Running MUSCLE aligner on amplicon fasta files. Please wait...", file=sys.stderr)
run_muscle(bins, outdir=outdir)

print("INFO: Parsing alignments to CIGAR strings")
cigars = parse_alignments(bins, mask=mask, min_homopolymer_length=args.polyN, outdir=outdir)
if not cigars:
    print("ERROR: could not determine CIGAR strings", file=sys.stderr)
    sys.exit(1)

write_cigar_strings(cigars, args.out)
print(f"INFO: Wrote ASV->CIGAR to {args.out}")