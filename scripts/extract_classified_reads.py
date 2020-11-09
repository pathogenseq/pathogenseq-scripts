import sys
import argparse
from collections import defaultdict
from tqdm import tqdm
import fastq2matrix as fm
from uuid import uuid4


def main(args):

    nodes = defaultdict(set)
    sys.stderr.write("Loading taxonomy\n")
    for l in tqdm(open(fm.filecheck(args.tax_dump))):
        row = l.strip().split()
        nodes[row[2]].add(row[0])

    def flatten(d):
        v = [[i] if not isinstance(i, list) else flatten(i) for i in d]
        return [i for b in v for i in b]

    def get_tax(t):
        if len(nodes[t])==0:
            return [t]
        return [t] + flatten([get_tax(sub_t) for sub_t in nodes[t]])

    sys.stderr.write("Extracting read names\n")
    args.tmp_file = str(uuid4())
    with open(args.tmp_file,"w") as O:
        if args.exclude:
            tax_tree = set(flatten([get_tax(x) for x in args.exclude.split(",")]))
            for l in tqdm(open(fm.filecheck(args.kraken2_output))):
                row = l.strip().split()
                if row[2] not in tax_tree:
                    O.write("%s\n" % row[1])
        else:
            tax_tree = set(flatten([get_tax(x) for x in args.extract.split(",")]))
            for l in tqdm(open(fm.filecheck(args.kraken2_output))):
                row = l.strip().split()
                if row[2] in tax_tree:
                    O.write("%s\n" % row[1])


    sys.stderr.write("Writing filtered fastq files\n")
    fm.filecheck(args.R1)
    args.R1_filt = args.R1.replace(".fastq.gz","").replace(".fq.gz","").replace(".fastq","") + "kraken_filt.fastq.gz"
    fm.run_cmd("seqtk subseq %(R1)s %(tmp_file)s | gzip -c > %(R1_filt)s" % vars(args))

    if args.R2:
        fm.filecheck(args.R2)
        args.R2_filt = args.R2.replace(".fastq.gz","").replace(".fq.gz","").replace(".fastq","") + "kraken_filt.fastq.gz"
        fm.run_cmd("seqtk subseq %(R2)s %(tmp_file)s | gzip -c > %(R2_filt)s" % vars(args))



parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-1','--R1',type=str,help='Forward reads',required=True)
parser.add_argument('-2','--R2',type=str,help='Reverse reads')
parser.add_argument('--kraken2-output',type=str,help='Output file from kraken2',required=True)
parser.add_argument('--tax-dump',type=str,help='Tax dump file',required=True)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--extract',type=str,help='Tax ID to extract')
group.add_argument('--exclude',type=str,help='Tax ID to exclude')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
