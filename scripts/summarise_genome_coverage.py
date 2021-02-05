import sys
import argparse
import fastq2matrix as fm
from collections import defaultdict
import csv
import os.path


try:
    import statsmodels.stats.weightstats
except:
    quit("\nCan't find python module statsmodels. Please install with using:\n\nconda install statsmodels\n")

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    quit("\nCan't find %s. Please install with:\n\nconda install bedtools.\n" % program)

def main(args):
    which("bedtools")
    if not os.path.isfile(args.bam+".genomecov.txt"):
        fm.run_cmd("bedtools genomecov -ibam %(bam)s > %(bam)s.genomecov.txt" % vars(args))
    dp = defaultdict(dict)
    for l in open(args.bam+".genomecov.txt"):
        row = l.strip().split()
        dp[row[0]][int(row[1])] = {"freq":int(row[2]),"fraction":float(row[4])}

    with open(args.out, "w") as O:
        writer = csv.DictWriter(O,fieldnames=["chrom","mean","std","dp_0","dp_5","dp_10"])
        writer.writeheader()
        print()
        for chrom in dp:
            d1 = statsmodels.stats.weightstats.DescrStatsW(list(dp[chrom].keys()),[x["freq"] for x in dp[chrom].values()])
            # import pdb; pdb.set_trace()
            res = {
                "chrom":chrom,
                "mean": d1.mean,
                "std": d1.std,
                "dp_0": (1-sum([dp[chrom][x]["fraction"] for x in [0]])) * 100,
                "dp_5": (1-sum([dp[chrom][x]["fraction"] for x in range(6)])) * 100,
                "dp_10": (1-sum([dp[chrom][x]["fraction"] for x in range(11)])) * 100
            }
            writer.writerow(res)


parser = argparse.ArgumentParser(description='Coverage summary',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam',help='BAM/CRAM file',required=True)
parser.add_argument('--out',help='Ouput file (csv format)',required=True)
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
