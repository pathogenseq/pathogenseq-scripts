import argparse
import csv
import pathogenprofiler as pp
import os
from tqdm import tqdm
from collections import defaultdict

###### Functions #######

def nofile(filename):
    """
    Return True if file does not exist
    """
    if not os.path.isfile(filename):
        return True
    else:
        return False

def filecheck(filename):
    """
    Check if file is there and quit if it isn't
    """
    if filename=="/dev/null":
        return filename
    elif not os.path.isfile(filename):
        sys.stderr.write("Can't find %s\n" % filename)
        exit(1)
    else:
        return filename

def index_bcf(bcffile,threads=1,overwrite=False):
    """
    Indexing a bam file
    """
    cmd = "bcftools index --threads %s -f %s" % (threads,bcffile)
    if filecheck(bcffile):
        if nofile(bcffile+".csi"):
            pp.run_cmd(cmd)
        elif os.path.getmtime(bcffile+".csi")<os.path.getmtime(bcffile) or overwrite:
            pp.run_cmd(cmd)


###### Classes #######

##### VCF class #####
class vcf_class:
    def __init__(self,filename,prefix=None,threads=1):
        self.samples = []
        self.filename = filename
        self.prefix = prefix
        self.threads = threads
        if prefix==None:
            if filename[-4:] == ".bcf":
                self.prefix = filename[:-4]
            elif filename[-5:] == ".gbcf":
                self.prefix = filename[:-5]
            elif filename[-7:] == ".vcf.gz":
                self.prefix = filename[:-7]
            elif filename[-8:] == ".gvcf.gz":
                self.prefix = filename[:-8]
            elif filename[-4:] == ".vcf":
                self.prefix = filename[:-4]
            else:
                self.prefix = filename
        else:
            self.prefix = prefix
        index_bcf(filename,self.threads)
        for l in pp.cmd_out("bcftools query -l %(filename)s" % vars(self)):
            self.samples.append(l.rstrip())

    def get_mean_genotype(self,outfile=None):
        self.outfile = outfile
        if self.outfile==None:
            self.outfile = self.prefix+".geno"
        O = open(self.outfile,"w")
        for l in tqdm(pp.cmd_out("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%TGT]\\n' %(filename)s" % vars(self))):
            row = l.rstrip().split()
            alts = row[3].split(",")
            for alt in alts:
                ref = "%s/%s" % (row[2],row[2])
                tmp = "%s/%s" % (alt,alt)
                genos = []
                for x in row[4:]:
                    if x==ref:
                        genos.append("0")
                    elif x==tmp:
                        genos.append("1")
                    else:
                        genos.append("NA")
                O.write("%s, %s, %s, %s\n" % (row[0]+"_"+row[1]+"_"+alt,row[2],alt,", ".join(genos)))
        O.close()

    def get_genesum(self,outfile=None):
        self.outfile = outfile
        if self.outfile==None:
            self.outfile = self.prefix+".gensum"
        genesum = defaultdict(lambda:defaultdict(int))
        O = open(self.outfile,"w")
        for l in tqdm(pp.cmd_out("bcftools query -f '[%%SAMPLE\\t%%GT\\t%%TBCSQ\\n]' %(filename)s" % vars(self))):
            row = l.split()
            #por4A    1/1    synonymous|Rv0002|gene1|protein_coding|+|109L|2378G>A    synonymous|Rv0002|gene1|protein_coding|+|109L|2378G>A
            info = row[2].split("|")
            if info[0]=="synonymous": continue
            if info[0][0]=="@": continue
            genesum[info[1]][row[0]]+=1
        for gene in genesum:
            O.write("%s\tNA\tNA\t%s\n" % (gene,"\t".join(str(genesum[gene][s]) for s in self.samples)))
        O.close()




###### Main #######

def main(args):
    vcf = vcf_class(args.vcf)
    # vcf.get_mean_genotype()
    if args.genes:
        vcf.get_genesum()
    geno_file = vcf.prefix+".geno"
    genesum_file = vcf.prefix+".genesum"
    meta = {}
    for s in vcf.samples:
        meta[s] = {}
    for row in csv.DictReader(open(args.csv)):
        for pheno in row.keys():
            if pheno=="id": continue
            if row['id'] not in meta: continue
            meta[row["id"]][pheno] = row[pheno]
    phenos = [x.rstrip() for x in open(args.phenos).readlines()]
    cmd_file = pp.get_random_file()
    X = open(cmd_file,"w")
    for pheno in phenos:
        pheno_file = "%s.pheno" % pheno
        if pheno not in row:
            pp.log("%s not in CSV file"%pheno,True)
        P = open(pheno_file,"w")
        P.write("\n".join([meta[s][pheno] if pheno in meta[s] else "NA" for s in vcf.samples]))
        P.close()
        X.write("gemma -p %s -g %s -gk 1 -o %s -maf 0.00005 -miss 0.99 && gemma  -lmm 1 -p %s -g %s  -k output/%s.cXX.txt  -o %s -maf 0.00005 -miss 0.99 && gemma  -lmm 1 -p %s -g %s  -k output/%s.cXX.txt  -o %s.genesum -notsnp\n" % (pheno_file,geno_file,pheno,pheno_file,geno_file,pheno,pheno,pheno_file,genesum_file,pheno,pheno))
    X.close()

    if args.preprocess:
        pp.log("Preprocessing finished\n", True)
    else:
        pp.run_cmd("cat %s | parallel -j %s" % (cmd_file,args.threads))


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('vcf', help='vcf file')
parser.add_argument('csv', help='CSV file')
parser.add_argument('phenos', help='Columns file')
parser.add_argument('--preprocess',default=False,action='store_true', help='Columns file')
parser.add_argument('--genes',default=False,action='store_true', help='Columns file')
parser.add_argument('--threads','-t',default=4,type=int, help='Columns file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
