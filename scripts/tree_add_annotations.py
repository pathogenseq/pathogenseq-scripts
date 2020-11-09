import argparse
import csv
from ete3 import Tree

def main(args):
	t = Tree(args.tree)
	annotation_dict = {}
	for row in csv.DictReader(open(args.csv)):
		annotation_dict[row["id"]] = row
	if args.annotations:
		annotations = args.annotations.split(",")
	else:
		annotations = list(list(annotation_dict.values())[0].keys())
		if "id" in annotations:	annotations.remove("id")
	for node in t.get_leaves():
		if args.strict:
			if node.name not in annotations:
				quit("%s not in annotations" % node.name)
		node.name = node.name+args.sep+args.sep.join([annotation_dict[node.name][d] for d in annotations])
		print(node.name)
		
	t.write(outfile=args.out)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--tree',type=str,help='File to search in',required=True)
parser.add_argument('--csv',type=str,help='File to search in',required=True)
parser.add_argument('--out',type=str,help='File to search in',required=True)
parser.add_argument('--annotations',type=str,help='Comma separated list of annotations to include')
parser.add_argument('--strict',action="store_true",help='Quit if ID not found')
parser.add_argument('--sep',type=str,default="_",help='Seperator for file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
