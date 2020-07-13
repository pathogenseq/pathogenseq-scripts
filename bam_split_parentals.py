#! /usr/bin/env python
import sys
from bisect import bisect_left, bisect_right
import random
import subprocess
import argparse

_version = "3.0.1"


def run_cmd(cmd,verbose=1):
    sys.stderr.write("\nRunning command:\n%s\n" % cmd)
    stdout = open("/dev/null","w")
    stderr = open("/dev/null","w")
    res = subprocess.call(cmd,shell=True,stderr = stderr,stdout = stdout)
    stderr.close()
    if res!=0:
        sys.stderr.write("Command Failed! Please Check!")
        exit(1)

def main(args):
    prefix = args.bam.replace(".bam","")
    run_cmd("samtools view -F 2048 %s.bam -h | awk '$1 ~/^@/ || ($9<2000 && $9>-2000)' | paftools.js sam2paf -L - > %s.paf" % (prefix,prefix))
    run_cmd("paftools.js view -f maf %s.paf > %s.maf" % (prefix,prefix))
    mixed_calls_dict = {}
    run_cmd("freebayes -f %s %s.bam | bcftools view -Oz -o %s.freebayes.vcf.gz" % (args.ref,prefix,prefix))
    run_cmd("bcftools view -a -M2 -m2 %s.freebayes.vcf.gz | bcftools query -f '%%POS\\t%%REF,%%ALT\\t[%%AD]\\n' > %s.af.txt" % (prefix,prefix))

    for l in open("%s.af.txt" % prefix):
        row = l.rstrip().split()
        pos = int(row[0])
        alts = row[1].split(",")
        ads = [int(x) for x in row[2].split(",")]
        if sum(ads)==0: continue
        afs = [x/sum(ads) for x in ads]#int(row[3])/(int(row[2])+int(row[3])) if len(alts)==1 else int(row[4])/(int(row[3])+int(row[4]))
        if afs[0]<0.9 and afs[0]>0.1:
            mixed_calls_dict[pos] = {alts[i]:afs[i] for i in range(len(alts))}
    mixed_call_list = list(mixed_calls_dict.keys())
    mixed_call_list_len = len(mixed_call_list)
    F = open("%s.maf" % prefix)
    l = F.readline()
    i = 0
    O = {0:open("%s.minor.txt" % prefix,"w"),1:open("%s.major.txt" % prefix,"w")}
    major_reads = set()
    minor_reads = set()
    other_reads = set()
    while 1:
        i+=1
        F.readline()
        F.readline()
        l1 = F.readline().rstrip()
        if l1 == "": break
        l2 = F.readline().rstrip()
        row1 = l1.split()
        row2 = l2.split()
        read_name = row2[1].replace("/1","").replace("/2","")
        start = int(row1[2])+1
        end = start+int(row1[3])
        start_index = bisect_left(mixed_call_list,start)
        end_index = bisect_left(mixed_call_list,end)
        mixed_calls = mixed_call_list[start_index:end_index]

        #if read_name=="M01637:69:000000000-J2465:1:2102:9090:13026":
        #    import pdb; pdb.set_trace()
        if len(mixed_calls)==0:
            other_reads.add(read_name)
        else:

            offset = mixed_calls[0] - start
            allele = row2[6][offset]

            if allele not in mixed_calls_dict[mixed_calls[0]]:
                O[random.randint(0,1)].write(read_name+"\n")
            elif mixed_calls_dict[mixed_calls[0]][allele]<0.5:
                minor_reads.add(read_name)
            else:
                major_reads.add(read_name)

    for read_name in major_reads:
        O[1].write(read_name+"\n")
    for read_name in minor_reads:
        O[0].write(read_name+"\n")
    for read_name in other_reads-major_reads-minor_reads:
        O[random.randint(0,1)].write(read_name+"\n")


    O[0].close()
    O[1].close()

    run_cmd("gatk FilterSamReads -I %s.bam -O %s.major.bam --FILTER includeReadList --READ_LIST_FILE %s.major.txt" % (prefix,prefix,prefix))
    run_cmd("gatk FilterSamReads -I %s.bam -O %s.minor.bam --FILTER includeReadList --READ_LIST_FILE %s.minor.txt" % (prefix,prefix,prefix))
    run_cmd("samtools index %s.major.bam" % prefix)
    run_cmd("samtools index %s.minor.bam" % prefix)
    run_cmd("freebayes -f %s %s.major.bam | bcftools view -Oz -o %s.major.freebayes.vcf.gz" % (args.ref,prefix,prefix))
    run_cmd("freebayes -f %s %s.minor.bam | bcftools view -Oz -o %s.minor.freebayes.vcf.gz" % (args.ref,prefix,prefix))

parser = argparse.ArgumentParser(description='Bam splitting pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam',type=str,help='The bam file you would like to split',required=True)
parser.add_argument('--ref',type=str,help='The reference file',required=True)
#parser.add_argument('--no-clean',action="store_true",help='Don\'t clean up temp files')
parser.add_argument('--version', action='version', version=_version)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
