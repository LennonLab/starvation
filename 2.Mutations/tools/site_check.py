"""Check three indel calls against the deposited reads.

Clone 4 of the endpoint total fraction carries three of the reported indels --
a 49 bp deletion plus a 1 bp insertion in recX, and a frameshift in ywcC --
so its calls were checked against the raw reads, together with the other two
ywcC carriers (clones 13 and 38), which double as negative controls.

No aligner is needed. A read counts for a site when it contains both 20-bp
reference anchors that flank the variable stretch, on either strand; the
stretch between them is compared with the reference and mutant haplotypes.
Each anchor is unique in the genome (printed first), so repeats cannot
confuse the counts.

Inputs live in data/sra/, which git ignores (~470 MB). If the folder is
missing, fetch them from the 2.Mutations folder:

Reads (PRJNA525711; interleaved pairs, about 120-180 MB each):

    LTMP_4  SRR41306380    LTMP_13  SRR41306373    LTMP_38  SRR41307095

    curl -L -o data/sra/LTMP_4.fq.gz \
      "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR41306380"

Reference: NC_000964.3 (B. subtilis 168) as FASTA:

    curl -o data/sra/NC_000964.3.fa "https://eutils.ncbi.nlm.nih.gov/entrez/\
eutils/efetch.fcgi?db=nuccore&id=NC_000964.3&rettype=fasta&retmode=text"

data/sra/ also holds LTMP_90 (SRR41306341, 833 read pairs), the failed library
deposited in place of LTMP_1; it is kept for the record and not checked.

usage, from 2.Mutations:

    python3 tools/site_check.py data/sra/NC_000964.3.fa data/sra/LTMP_4.fq.gz \
        data/sra/LTMP_13.fq.gz data/sra/LTMP_38.fq.gz > output/read_check.txt

Standard library only. Takes about a minute per library.
"""
import gzip
import os
import sys
from collections import Counter

K = 20
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def load_genome(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def sites(g):
    """(name, start, end, {label: middle}) in 0-based half-open coordinates."""
    out = []
    # ywcC 3,922,952: G then an 8-T run (0-based 3922952..3922959)
    s, e = 3922952, 3922960
    out.append(("ywcC 3,922,952", s, e, {"ref T8": g[s:e], "+T (clone 4 call)": "T" * 9,
                                          "-T (clone 38 call)": "T" * 7}))
    # ywcC 3,922,748: C then an 8-T run (0-based 3922748..3922755)
    s, e = 3922748, 3922756
    out.append(("ywcC 3,922,748", s, e, {"ref T8": g[s:e], "-T (clone 13 call)": "T" * 7,
                                          "+T": "T" * 9}))
    # recX: 49-bp deletion after the C at 925,845 (0-based 925845..925893)
    # and a C inserted after the A at 925,897 (0-based 925896)
    s, e = 925845, 925900
    ref = g[s:e]
    dele = g[925894:e]
    ins = g[s:925897] + "C" + g[925897:e]
    both = g[925894:925897] + "C" + g[925897:e]
    out.append(("recX 925,845-925,897", s, e, {"ref": ref, "del49 only": dele, "insC only": ins,
                                                "del49 + insC (clone 4 call)": both}))
    return out


def anchor_hits(g, a):
    n, i = 0, g.find(a)
    while i != -1:
        n, i = n + 1, g.find(a, i + 1)
    j = g.find(rc(a))
    while j != -1:
        n, j = n + 1, g.find(rc(a), j + 1)
    return n


def main():
    g = load_genome(sys.argv[1])
    st = sites(g)
    anchors = [(g[s - K:s], g[e:e + K]) for _, s, e, _ in st]
    for (name, _, _, _), (a, b) in zip(st, anchors):
        print(f"# {name}: left anchor in genome x{anchor_hits(g, a)}, right anchor x{anchor_hits(g, b)}")
    for fq in sys.argv[2:]:
        counts = [Counter() for _ in st]
        n = 0
        with gzip.open(fq, "rt") as fh:
            for k, line in enumerate(fh):
                if k % 4 != 1:
                    continue
                n += 1
                r = line.strip()
                for idx, (a, b) in enumerate(anchors):
                    for read, A, B in ((r, a, b), (r, rc(b), rc(a))):
                        i = read.find(A)
                        if i == -1:
                            continue
                        j = read.find(B, i + K)
                        if j == -1:
                            continue
                        mid = read[i + K:j]
                        if A != a:
                            mid = rc(mid)
                        lab = next((l for l, m in st[idx][3].items() if m == mid), None)
                        if lab is None:
                            lab = f"other (T{len(mid)})" if set(mid) == {"T"} else "other"
                        counts[idx][lab] += 1
        print(f"\n== {os.path.basename(fq)}: {n:,} reads")
        for (name, _, _, labs), c in zip(st, counts):
            tot = sum(c.values())
            parts = ", ".join(f"{l}: {c[l]}" for l in list(labs) + sorted(set(c) - set(labs)) if c[l])
            print(f"  {name:22s} spanning reads {tot:4d} | {parts}")


if __name__ == "__main__":
    main()
