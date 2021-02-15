import argparse
import csv
from collections import OrderedDict
class fasta:
	"""
	Class to represent fasta seuqnces in a python dict.

	Args:
		filename(str): Location of the fasta file

	Returns:
		fasta: A fasta class object
	"""
	def __init__(self,filename):
		fa_dict = OrderedDict()
		seq_name = ""
		self.fa_file = filename
		for l in open(filename):
			line = l.rstrip()
			if line.startswith(">"):
				seq_name = line[1:].split()[0]
				fa_dict[seq_name] = []
			else:
				fa_dict[seq_name].append(line)
		result = {}
		counter = 0
		sum_length = {}
		for seq in fa_dict:
			result[seq] = "".join(fa_dict[seq])
			result[seq] = result[seq].upper()
			sum_length[(counter+1,counter+len(result[seq]))] = seq
			counter = counter+len(result[seq])
		self.sum_length = sum_length
		self.fa_dict = result


def main(args):
	fa = fasta(args.fasta).fa_dict
	annotation_dict = {}
	for row in csv.DictReader(open(args.csv)):
		annotation_dict[row[args.id_key]] = row
	if args.annotations:
		annotations = args.annotations.split(",")
	else:
		annotations = list(list(annotation_dict.values())[0].keys())
		if "id" in annotations:	annotations.remove("id")
	new_fa = {}

	for seq in fa:
		new_name = seq+args.sep+args.sep.join([annotation_dict[seq][d] for d in annotations])
		new_fa[new_name] = fa[seq]

	with open(args.out,"w") as O:
		for seq in new_fa:
			O.write(">%s\n%s\n" % (seq,new_fa[seq]))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fasta',type=str,help='File to search in',required=True)
parser.add_argument('--csv',type=str,help='File to search in',required=True)
parser.add_argument('--out',type=str,help='File to search in',required=True)
parser.add_argument('--annotations',type=str,help='Comma separated list of annotations to include')
parser.add_argument('--sep',type=str,default="_",help='Seperator for file')
parser.add_argument('--id-key',type=str,default="id",help='ID key in metadata')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
