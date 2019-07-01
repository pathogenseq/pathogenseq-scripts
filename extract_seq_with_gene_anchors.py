#! /usr/bin/env python
import sys
import subprocess
import argparse
import re
import random
from Bio.Blast import NCBIXML
rand_generator = random.SystemRandom()

def run_cmd(cmd,verbose=1):
	sys.stderr.write("\nRunning command:\n%s\n" % cmd)
	stdout = open("/dev/null","w")
	stderr = open("/dev/null","w")
	res = subprocess.call(cmd,shell=True,stderr = stderr,stdout = stdout)
	stderr.close()
	if res!=0:
		sys.stderr.write("Command Failed! Please Check!")
		exit(1)

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

class gene_class:
	def __init__(self,name,locus_tag,strand,start,end,length,chrom):
		self.name = name
		self.locus_tag = locus_tag
		self.strand = strand
		self.start = start
		self.end = end
		self.length = length
		self.chrom = chrom

class blast_aln:
	def __init__(self,records):
		best_rec = None
		best_aln = None
		best_hsp = None
		best_eval = 100
		for rec in records:
			for aln in rec.alignments:
				for hsp in aln.hsps:
					if hsp.expect<best_eval:
						best_rec = rec
						best_aln = aln
						best_hsp = hsp
						best_eval = hsp.expect
		self.seq = best_aln.hit_def
		self.start = best_hsp.sbjct_start
		self.end = best_hsp.sbjct_end
		self.identities = best_hsp.identities
		self.proportion = best_hsp.align_length/rec.query_length

def load_gff(gff,gene_name_prefix="Name",gene_id_prefix="gene"):
	genes = []
	for l in open(gff):
		#ID=gene:Rv0150c;biotype=protein_coding;description=Conserved hypothetical protein;gene_id=Rv0150c;logic_name=ena
		if l[0]=="#": continue
		fields = l.rstrip().split()
		if fields[2]!="gene": continue
		strand = fields[6]
		p1 = int(fields[3])
		p2 = int(fields[4])
		gene_length = p2-p1+1
		re_obj = re.search("%s=([a-zA-Z0-9\.\-\_]+)" % gene_id_prefix,l)
		gene_name = re_obj.group(1) if re_obj else "NA"
		re_obj = re.search("%s=([a-zA-Z0-9\.\-\_]+)" % gene_name_prefix,l)
		locus_tag = re_obj.group(1) if re_obj else "NA"
		start = p1 if strand=="+" else p2
		end =  p2 if strand=="+" else p1
		chrom = fields[0]
		tmp = gene_class(gene_name,locus_tag,strand,start,end,gene_length,chrom)
		genes.append(tmp)
	return genes

def extract_gene(ref,gene,outfile=None):
	outfile = get_random_file() if not outfile else outfile
	run_cmd("samtools faidx %s" % ref)
	if gene.strand=="+":
		run_cmd("samtools faidx %s %s:%s-%s > %s" % (ref,gene.chrom,gene.start,gene.end,outfile))
	else:
		run_cmd("samtools faidx %s %s:%s-%s > %s" % (ref,gene.chrom,gene.end,gene.start,outfile))
	return outfile

def extract_seq(ref,chrom,start,end,outfile=None,revcom=False):
	outfile = get_random_file() if not outfile else outfile
	revcom_flag = " --reverse-complement " if revcom else " "
	run_cmd("samtools faidx %s" % ref)
	run_cmd("samtools faidx %s \"%s:%s-%s\" %s > %s" % (ref,chrom,start,end,revcom_flag,outfile))
	return outfile

def blast_seq(subject,query,result):
	run_cmd("blastn -task blastn -subject %s -query %s -outfmt 5 > %s" % (subject,query,result))
	blast_records = NCBIXML.parse(open(result))
	return blast_aln(blast_records)

def main(args):
	genes = load_gff(args.gff,gene_name_prefix=args.gene_name_prefix,gene_id_prefix=args.gene_id_prefix)
	tmp = [g for g in genes if g.locus_tag==args.gene1]
	if len(tmp)==0: exit("Can't find %s in gff" % args.gene1)
	gene1 = tmp[0]
	tmp = [g for g in genes if g.locus_tag==args.gene2]
	if len(tmp)==0: exit("Can't find %s in gff" % args.gene2)
	gene2 = tmp[0]
	gene1_fasta = extract_gene(args.subject,gene1)
	gene2_fasta = extract_gene(args.subject,gene2)

	gene1_blast_result = get_random_file()
	gene2_blast_result = get_random_file()

	gene1_hsp = blast_seq(args.query,gene1_fasta,gene1_blast_result)
	gene2_hsp = blast_seq(args.query,gene2_fasta,gene2_blast_result)

	if gene1_hsp.seq!=gene2_hsp.seq:
		quit("Genes are on different sequences: %s and %s" % (gene1_hsp.seq,gene2_hsp.seq))
	elif (gene1_hsp.proportion<0.9 or gene1_hsp.proportion>1.1):
		quit("Alignment length to query length proportion is %s for %s" % (gene1_hsp.proportion,args.gene1))
	elif (gene2_hsp.proportion<0.9 or gene2_hsp.proportion>1.1):
		quit("Alignment length to query length proportion is %s for %s" % (gene2_hsp.proportion,args.gene2))

	region_start = min([gene1_hsp.start,gene1_hsp.end,gene2_hsp.start,gene2_hsp.end])
	region_end = max([gene1_hsp.start,gene1_hsp.end,gene2_hsp.start,gene2_hsp.end])

	prefix = args.query.replace(".fasta","").replace(".fa","")
	outfile = "%s_%s_%s.fasta" % (prefix,args.gene1,args.gene2)
	extract_seq(args.query,gene1_hsp.seq,region_start,region_end,outfile=outfile,revcom=args.revcom)

	if args.mafft:
		region_start = min([gene1.start,gene1.end,gene2.start,gene2.end])
		region_end = max([gene1.start,gene1.end,gene2.start,gene2.end])
		ref_fasta = extract_seq(args.subject,gene1.chrom,region_start,region_end,revcom=args.revcom)
		unaln = get_random_file()
		aln = "%s_%s_%s.aln" % (prefix,args.gene1,args.gene2)
		run_cmd("cat %s %s > %s" % (ref_fasta,outfile,unaln))
		run_cmd("mafft --clustalout %s  > %s" % (unaln,aln))


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('subject',type=str,help='NGS Platform')
parser.add_argument('gff',default="pp-results",type=str,help='NGS Platform')
parser.add_argument('gene1',type=str,help='NGS Platform')
parser.add_argument('gene2',type=str,help='NGS Platform')
parser.add_argument('query',type=str,help='NGS Platform')
parser.add_argument('--revcom',action="store_true",help='NGS Platform')
parser.add_argument('--mafft',action="store_true",help='NGS Platform')
parser.add_argument('--gene-name-prefix',default="gene",help='Gene name prefix')
parser.add_argument('--gene-id-prefix',default="Name",help='Gene ID prefix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
