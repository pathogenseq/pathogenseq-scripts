#! /usr/bin/env python
import subprocess
import argparse
from collections import OrderedDict
import re
import sys
def run_cmd(cmd,verbose=1,target=None):
	"""
	Wrapper to run a command using subprocess with 3 levels of verbosity and automatic exiting if command failed
	"""
	if target and filecheck(target): return True
	cmd = "set -u pipefail; " + cmd
	if verbose==2:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stdout = open("/dev/stdout","w")
		stderr = open("/dev/stderr","w")
	elif verbose==1:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stdout = open("/dev/null","w")
		stderr = open("/dev/null","w")
	else:
		stdout = open("/dev/null","w")
		stderr = open("/dev/null","w")

	res = subprocess.call(cmd,shell=True,stderr = stderr,stdout = stdout)
	stderr.close()
	if res!=0:
		sys.stderr.write("Command Failed! Please Check!")
		exit(1)


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
	def get_ref_variants(self,refseq,prefix,gff=None):
		self.refseq = refseq
		self.prefix = prefix
		self.gff = gff
		self.csq_cmd = "| bcftools csq -f %(refseq)s -g %(gff)s" % vars(self) if gff else ""
		run_cmd("minimap2 %(refseq)s %(fa_file)s --cs | sort -k6,6 -k8,8n | paftools.js call -l 100 -L 100 -f %(refseq)s -s %(prefix)s - %(csq_cmd)s |  bcftools view -Oz -o %(prefix)s.vcf.gz" % vars(self))
		return "%s.vcf.gz" % self.prefix


def main(args):
	assembly = fasta(args.assembly)
	if args.assembly[-6:]==".fasta":
		prefix = args.assembly.replace(".fasta","")
	elif args.assembly[-3:]==".fa":
		prefix = args.assembly.replace(".fa","")
	else:
		prefix = args.assembly
	assembly.get_ref_variants(args.ref,prefix)

parser = argparse.ArgumentParser(description='Assembly to VCF',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('ref',help='BAM file')
parser.add_argument('assembly',help='GFF file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
