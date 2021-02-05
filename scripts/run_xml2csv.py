import csv
import argparse
import xmltodict


def main(args):
    data = xmltodict.parse(open(args.xml).read())
    for run in data["ROOT"]["RUN"]:
        print(run["@accession"])

parser = argparse.ArgumentParser(description='XXX pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--xml',help='VCF file',required=True)
parser.add_argument('--out',help='VCF file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
