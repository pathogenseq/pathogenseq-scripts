import fastq2matrix as fm
import argparse
import random
rand_generator = random.SystemRandom()

def main(args):
    randint = rand_generator.randint(1,999999)

    window_cmd = "bedtools makewindows -n %(chunks)s -g %(ref)s.fai | awk '{print $1\":\"$2+1\"-\"$3\" \"$1\"_\"$2+1\"_\"$3}'" % vars(args)

    cmd_to_run = "\"bcftools view --threads %s -r {1} %s -Ou | %s | bcftools view  -Oz -o %s_{2}.vcf.gz\"" % (args.compression_threads,args.vcf,args.cmd,randint)
    fm.run_cmd(f"{window_cmd} | parallel --bar -j {args.threads} --col-sep \" \" {cmd_to_run}", verbose=2)
    fm.run_cmd("bcftools concat -Oz -o %s `%s | awk '{print \"%s_\"$2\".vcf.gz\"}'`" % (args.out,window_cmd,randint))
    fm.run_cmd("rm `%s | awk '{print \"%s_\"$2\".vcf.gz*\"}'`" % (window_cmd,randint))


parser = argparse.ArgumentParser(description='XXX pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',help='VCF file',required=True)
parser.add_argument('--ref',help='reference',required=True)
parser.add_argument('--cmd',help='Command',required=True)
parser.add_argument('--out',help='Output vcf file',required=True)
parser.add_argument('--chunks',default=20,help='VCF file')
parser.add_argument('--threads',default=10,help='VCF file')
parser.add_argument('--compression-threads',default=4,help='VCF file')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
