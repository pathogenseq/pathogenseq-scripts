import sys
import os
import subprocess
from collections import defaultdict
import random
import argparse
rand_generator = random.SystemRandom()

def cmd_out(cmd,verbose=1):
    cmd = "set -u pipefail; " + cmd
    if verbose==2:
        sys.stderr.write("\nRunning command:\n%s\n" % cmd)
        stderr = open("/dev/stderr","w")
    elif verbose==1:
        sys.stderr.write("\nRunning command:\n%s\n" % cmd)
        stderr = open("/dev/null","w")
    else:
        stderr = open("/dev/null","w")
    try:
        res = subprocess.Popen(cmd,shell=True,stderr = stderr,stdout=subprocess.PIPE)
        for l in res.stdout:
            yield l.decode().rstrip()
    except:
        sys.stderr.write("Command Failed! Please Check!")
        exit(1)
    stderr.close()

def get_random_file(prefix = None,extension=None):
    randint = rand_generator.randint(1,999999)
    if prefix:
        if extension:
            return "%s.%s.%s" % (prefix,randint,extension)
        else:
            return "%s.%s.txt" % (prefix,randint)
    else:
        if extension:
            return "%s.tmp.%s" % (randint,extension)
        else:
            return "%s.tmp.txt" % (randint)

def log(msg,ext=False):
    sys.stderr.write("\n"+str(msg)+"\n")
    if ext:
        exit(1)

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
        print("Command Failed! Please Check!")
        exit(1)


def write_spacers():
    tmp = get_random_file()
    open(tmp,"w").write("""> Spacer1
ATAGAGGGTCGCCGGCTCTGGATCA
> Spacer2
CCTCATGCTTGGGCGACAGCTTTTG
> Spacer3
CCGTGCTTCCAGTGATCGCCTTCTA
> Spacer4
ACGTCATACGCCGACCAATCATCAG
> Spacer5
TTTTCTGACCACTTGTGCGGGATTA
> Spacer6
CGTCGTCATTTCCGGCTTCAATTTC
> Spacer7
GAGGAGAGCGAGTACTCGGGGCTGC
> Spacer8
CGTGAAACCGCCCCCAGCCTCGCCG
> Spacer9
ACTCGGAATCCCATGTGCTGACAGC
> Spacer10
TCGACACCCGCTCTAGTTGACTTCC
> Spacer11
GTGAGCAACGGCGGCGGCAACCTGG
> Spacer12
ATATCTGCTGCCCGCCCGGGGAGAT
> Spacer13
GACCATCATTGCCATTCCCTCTCCC
> Spacer14
GGTGTGATGCGGATGGTCGGCTCGG
> Spacer15
CTTGAATAACGCGCAGTGAATTTCG
> Spacer16
CGAGTTCCCGTCAGCGTCGTAAATC
> Spacer17
GCGCCGGCCCGCGCGGATGACTCCG
> Spacer18
CATGGACCCGGGCGAGCTGCAGATG
> Spacer19
TAACTGGCTTGGCGCTGATCCTGGT
> Spacer20
TTGACCTCGCCAGGAGAGAAGATCA
> Spacer21
TCGATGTCGATGTCCCAATCGTCGA
> Spacer22
ACCGCAGACGGCACGATTGAGACAA
> Spacer23
AGCATCGCTGATGCGGTCCAGCTCG
> Spacer24
CCGCCTGCTGGGTGAGACGTGCTCG
> Spacer25
GATCAGCGACCACCGCACCCTGTCA
> Spacer26
CTTCAGCACCACCATCATCCGGCGC
> Spacer27
GGATTCGTGATCTCTTCCCGCGGAT
> Spacer28
TGCCCCGGCGTTTAGCGATCACAAC
> Spacer29
AAATACAGGCTCCACGACACGACCA
> Spacer30
GGTTGCCCCGCGCCCTTTTCCAGCC
> Spacer31
TCAGACAGGTTCGCGTCGATCAAGT
> Spacer32
GACCAAATAGGTATCGGCGTGTTCA
> Spacer33
GACATGACGGCGGTGCCGCACTTGA
> Spacer34
AAGTCACCTCGCCCACACCGTCGAA
> Spacer35
TCCGTACGCTCGAAACGCTTCCAAC
> Spacer36
CGAAATCCAGCACCACATCCGCAGC
> Spacer37
CGCGAACTCGTCCACAGTCCCCCTT
> Spacer38
CGTGGATGGCGGATGCGTTGTGCGC
> Spacer39
GACGATGGCCAGTAAATCGGCGTGG
> Spacer40
CGCCATCTGTGCCTCATACAGGTCC
> Spacer41
GGAGCTTTCCGGCTTCTATCAGGTA
> Spacer42
ATGGTGGGACATGGACGAGCGCGAC
> Spacer43
CGCAGAATCGCACCGGGTGCGGGAG
""")
    return tmp

def main(args):
    spacers = ["0" for i in range(43)]
    spacersfile = write_spacers()
    for l in cmd_out("blastn -task blastn -query %s -subject %s -word_size 25 -outfmt 6" % (spacersfile,args.fasta)):
        row = l.rstrip().split()
        spacers[int(row[0].replace("Spacer",""))-1] = "1"


    octal = []
    for i in range(0,40,3):
        tmp = "".join([spacers[i],spacers[i+1],spacers[i+2]])
        if tmp=="000":octal.append("0")
        elif tmp=="001":octal.append("1")
        elif tmp=="010":octal.append("2")
        elif tmp=="011":octal.append("3")
        elif tmp=="100":octal.append("4")
        elif tmp=="101":octal.append("5")
        elif tmp=="110":octal.append("6")
        elif tmp=="111":octal.append("7")
        else:ps.log("Don't know what to do with %s" % tmp,ext=T)

    octal.append("0" if spacers[42]=="0" else "1")
    sitvit_str = "".join(["n" if x=="1" else "o" for x in spacers])
    binary_str = "".join(spacers)
    octal_str = "".join(octal)

    print("%s\t%s\t%s\t%s" % (args.fasta,sitvit_str,binary_str,octal_str))
    os.remove(spacersfile)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fasta',help='VCF file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
