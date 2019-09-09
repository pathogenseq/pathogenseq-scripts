#! /usr/bin/python
import sys
import subprocess
import os.path
import argparse
from collections import defaultdict
import random
rand_generator = random.SystemRandom()

def get_random_file(prefix = None,extension=None):
	randint = rand_generator.randint(1,999999)
	if prefix:
		if extension:
			return "%s.%s%s" % (prefix,randint,extension)
		else:
			return "%s.%s.txt" % (prefix,randint)
	else:
		if extension:
			return "%s.tmp%s" % (randint,extension)
		else:
			return "%s.tmp.txt" % (randint)


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


def fa2dict(filename):
        fa_dict = {}
        seq_name = ""
        for l in open(filename):
                line = l.rstrip()
                if line[0] == ">":
                        seq_name = line[1:].split()[0]
                        fa_dict[seq_name] = []
                else:
                        fa_dict[seq_name].append(line)
        result = {}
        for seq in fa_dict:
                result[seq] = "".join(fa_dict[seq])
        return result

def revcom(x):
	return "".join([{"A":"T","C":"G","G":"C","T":"A"}[x[i]] for i in range(len(x)-1,-1,-1)])



def main(args):
	params = vars(args)
	motifs = fa2dict(args.primers)
	results = {}
	refgenome = fa2dict(args.genome)
	for mname in motifs:
		results[mname] = {}
		params["tmp_results"] = get_random_file()
		params["motif"] = motifs[mname]
		run_cmd("fuzznuc -sequence %(genome)s -pattern %(motif)s -outfile /dev/stdout -complement | grep -A1 '# Sequence:' > %(tmp_results)s" % params,0)
		tmp_chrom = None
		for l in open(params["tmp_results"]):
			row = l.rstrip().split()
			if l[0]!="#": continue
			if row[1]=="Sequence:":
				tmp_chrom = row[2]
			elif row[1]=="HitCount:":
				results[mname][tmp_chrom] = row[2]

	sys.stdout.write("Primer\t%s\n" % ("\t".join(refgenome.keys())))
	for m in motifs:
		sys.stdout.write("%s\t%s\n" % (m,"\t".join([results[m][s] if s in results[m] else "0" for s in refgenome])))


parser = argparse.ArgumentParser(description='Assembly to VCF',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--primers','-p',help='Fasta file with kmers',required=True)
parser.add_argument('--genome','-g',help='Reference genome',required=True)
parser.add_argument('--out','-o',help='Output file',required=True)
parser.add_argument('--threads','-t',default=2,help='Number of threads required file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
