import sys
import argparse
import random
import subprocess
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


def main(args):

	params = {"bam":args.bam,"ref":args.ref,"prefix":args.bam.split("/")[-1].replace(".bam","")}
	params["tmpfile"] = get_random_file()
	run_cmd("samtools depth -aa --reference %(ref)s %(bam)s | datamash mean 3 > %(tmpfile)s" % params)
	params["mean_depth"] = float(open(params["tmpfile"]).readline().strip())
	print("%(prefix)s\t%(mean_depth)s" % params)

parser = argparse.ArgumentParser(description='Assembly to VCF',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam',help='BAM file',required = True)
parser.add_argument('--ref',help='BAM file',required = True)

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
