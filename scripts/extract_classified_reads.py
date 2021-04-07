import sys
import argparse
from collections import defaultdict
from tqdm import tqdm
import fastq2matrix as fm
from uuid import uuid4
import os

def download_files(directory=None):
    sys.stderr.write("Downloading required files\n")
    import urllib.request
    if not directory:
        directory = "%s/.taxonkit/" % os.path.expanduser("~")
    if not os.path.isdir(directory):
        os.mkdir(directory)

    urllib.request.urlretrieve('ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz ','%s/taxdump.tar.gz' % directory )
    fm.run_cmd("tar -C %s -xvf %s/taxdump.tar.gz" % (directory,directory))


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

    return None

def check_programs(lst):
    for p in lst:
        sys.stderr.write("Checking if %s is in path..." % p)
        if which(p)==None:
            sys.stderr.write("No\n\nPlease install %s and try again\n" % p)
            quit()
        else:
            sys.stderr.write("Yes\n")

def main(args):
    check_programs(["taxonkit","seqtk"])
    if not os.path.isdir("%s/.taxonkit/" % os.path.expanduser("~")) or not os.path.isfile("%s/.taxonkit/nodes.dmp" % os.path.expanduser("~")):
        download_files()
    nodes = set()
    sys.stderr.write("Loading taxonomy\n")
    cmd = "taxonkit list --ids %s" % (args.extract if args.extract else args.exclude)
    for l in fm.cmd_out(cmd):
        if l=="": continue
        row = l.strip().split()
        nodes.add(row[0])
        

    sys.stderr.write("Extracting read names\n")
    args.tmp_file = str(uuid4())
    total_reads = 0
    kept_reads = 0
    
    with open(args.tmp_file,"w") as O:
        if args.exclude:
            for l in tqdm(open(fm.filecheck(args.kraken2_output))):
                total_reads+=1
                row = l.strip().split()
                if row[2] not in nodes:
                    O.write("%s\n" % row[1])
                    kept_reads+=1
        else:
            for l in tqdm(open(fm.filecheck(args.kraken2_output))):
                total_reads+=1
                row = l.strip().split()
                if row[2] in nodes:
                    O.write("%s\n" % row[1])
                    kept_reads+=1


    sys.stderr.write("Writing filtered fastq files\n")
    fm.filecheck(args.R1)
    args.R1_filt = args.R1.replace(".fastq.gz","").replace(".fq.gz","").replace(".fastq","") + ".kraken2_filt.fastq.gz"
    fm.run_cmd("seqtk subseq %(R1)s %(tmp_file)s | gzip -c > %(R1_filt)s" % vars(args))

    if args.R2:
        fm.filecheck(args.R2)
        args.R2_filt = args.R2.replace(".fastq.gz","").replace(".fq.gz","").replace(".fastq","") + ".kraken2_filt.fastq.gz"
        fm.run_cmd("seqtk subseq %(R2)s %(tmp_file)s | gzip -c > %(R2_filt)s" % vars(args))

    fm.rm_files([args.tmp_file])
    sys.stderr.write("\nKept %s/%s reads\n" % (kept_reads,total_reads))

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-1','--R1',type=str,help='Forward reads',required=True)
parser.add_argument('-2','--R2',type=str,help='Reverse reads')
parser.add_argument('--kraken2-output',type=str,help='Output file from kraken2',required=True)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--extract',type=str,help='Tax ID to extract')
group.add_argument('--exclude',type=str,help='Tax ID to exclude')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
