import sys
import argparse
import os
def main(args):
    with open(args.out,"w") as O:
        O.write("new\n")
        O.write("genome Mtb\n")
        O.write("goto %s\n" % args.region)

        for bam in args.bams:
            O.write("\n")
            O.write("new\n")
            O.write("load %s/%s\n" % (os.getcwd(),bam))
            O.write("snapshot %s.png\n" % bam)

parser = argparse.ArgumentParser(description='XXX pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--region',help='VCF file',required=True)
parser.add_argument('--out',help='VCF file',required=True)
parser.add_argument('bams',nargs='+',help='VCF file')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
