import os
import sys
import argparse
from tqdm import tqdm
import fastq2matrix as fm
import uuid

def c(x):
    return {"A":"T","C":"G","G":"C","T":"A","I":"D","D":"I"}[x]

def main(args):
    chromosomes = list(range(1,23))
    # chromosomes = [1]
    snp_data = {}
    for f in tqdm(["%s/chr%s.1kg.phase3.v5a.pvar" % (args.ref_dir,x) for x in chromosomes]):
        for l in open(f):
            if l[0]=="#": continue
            row = l.strip().split()
            snp_data[row[2]] = row


    tmp_prefix = str(uuid.uuid4())
    exclude_file = "%s.exclude.txt" % args.out
    new_bim_file = "%s.bim" % tmp_prefix
    log_file = "%s.fill_bim.log" % args.out

    EXCLUDE = open(exclude_file,"w")
    BIM = open(new_bim_file,"w")
    LOG = open(log_file,"w")

    for l in tqdm(open(args.bfile+".bim")):
        row = l.strip().split()
        rid = row[1]
        ref_snp_data = snp_data[row[1]] if rid in snp_data else None

        if args.remove_exm and "ex" in rid:
            LOG.write("%s\tExcluded: Variant starts with exm\n" % row[1])
            BIM.write("\t".join(row)+"\n")
            EXCLUDE.write(row[1]+"\n")
            continue


        if row[4]!="0" and row[5]!="0":
            if row[4]!="I" and row[4]!="D":
                LOG.write("%s\tOK: No change\n" % row[1])
                BIM.write("\t".join(row)+"\n")
            elif (row[5]=="I" or row[5]=="D") and ref_snp_data==None:
                LOG.write("%s\tExcluded: Indel not in ref\n" % row[1])
                BIM.write("\t".join(row)+"\n")
                EXCLUDE.write(row[1]+"\n")
            elif (row[4]=="I" or row[5]=="I") and ("," in ref_snp_data[4] or "," in ref_snp_data[3]):
                LOG.write("%s\tExcluded: More than one alt allele in ref\n" % row[1])
                BIM.write("\t".join(row)+"\n")
                EXCLUDE.write(row[1]+"\n")
            elif row[5]=="I" and len(ref_snp_data[3])>1:
                row[5]=ref_snp_data[3]
                row[4]=ref_snp_data[4]
                BIM.write("\t".join(row)+"\n")
                LOG.write("%s\tFilled indel to ref\t\n" % row[1])
            elif row[5]=="I" and len(ref_snp_data[4])>1:
                row[5]=ref_snp_data[4]
                row[4]=ref_snp_data[3]
                BIM.write("\t".join(row)+"\n")
                LOG.write("%s\tFilled indel to ref\t\n" % row[1])
            elif row[4]=="I" and len(ref_snp_data[4])>1:
                row[5]=ref_snp_data[3]
                row[4]=ref_snp_data[4]
                BIM.write("\t".join(row)+"\n")
                LOG.write("%s\tFilled indel to ref\t\n" % row[1])
            elif row[4]=="I" and len(ref_snp_data[3])>1:
                row[5]=ref_snp_data[4]
                row[4]=ref_snp_data[3]
                BIM.write("\t".join(row)+"\n")
                LOG.write("%s\tFilled indel to ref\t\n" % row[1])

            else:
                import pdb; pdb.set_trace()
        elif row[4]=="0" and row[5]=="0":
            LOG.write("%s\tExcluded: No ref or alt present\n" % row[1])
            BIM.write("\t".join(row)+"\n")
            EXCLUDE.write(row[1]+"\n")
        elif row[4]=="0":
            if not ref_snp_data:
                LOG.write("%s\tExcluded: SNP not present in ref\n" % row[1])
                BIM.write("\t".join(row)+"\n")
                EXCLUDE.write(row[1]+"\n")
            elif set([ref_snp_data[3],ref_snp_data[4]])==set(["A","T"]) or set([ref_snp_data[3],ref_snp_data[4]])==set(["C","G"]):
                LOG.write("%s\tExcluded: Ambiguous ref/alt strand\n" % row[1])
                BIM.write("\t".join(row)+"\n")
                EXCLUDE.write(row[1]+"\n")
            elif "," in ref_snp_data[4] or "," in ref_snp_data[3]:
                LOG.write("%s\tExcluded: More than one alt allele in ref\n" % row[1])
                BIM.write("\t".join(row)+"\n")
                EXCLUDE.write(row[1]+"\n")
            elif row[5]==ref_snp_data[3]:
                row[4]=ref_snp_data[4]
                BIM.write("\t".join(row)+"\n")
                LOG.write("%s\tOK: All ref\n" % row[1])
            elif row[5]==ref_snp_data[4]:
                row[4]=ref_snp_data[3]
                BIM.write("\t".join(row)+"\n")
                LOG.write("%s\tOK: All alt\n" % row[1])
            elif c(row[5])==ref_snp_data[3]:
                row[5]=ref_snp_data[3]
                row[4]=ref_snp_data[4]
                BIM.write("\t".join(row)+"\n")
                LOG.write("%s\tFlipped to be ref\t\n" % row[1])
            elif c(row[5])==ref_snp_data[4]:
                row[5]=ref_snp_data[4]
                row[4]=ref_snp_data[3]
                BIM.write("\t".join(row)+"\n")
                LOG.write("%s\tFlipped to be alt\t\n" % row[1])
            elif row[5]=="I" and len(ref_snp_data[3])>1:
                row[5]=ref_snp_data[3]
                row[4]=ref_snp_data[4]
                BIM.write("\t".join(row)+"\n")
                LOG.write("%s\tFilled indel to ref\t\n" % row[1])
            elif row[5]=="D" and len(ref_snp_data[4])>1:
                row[5]=ref_snp_data[3]
                row[4]=ref_snp_data[4]
                BIM.write("\t".join(row)+"\n")
                LOG.write("%s\tFilled indel to ref\t\n" % row[1])
            elif (row[5]!="I" and row[5]!="D") and len(ref_snp_data[3])>1:
                LOG.write("%s\tExcluded: Ref says indel but gt is SNP\n" % row[1])
                BIM.write("\t".join(row)+"\n")
                EXCLUDE.write(row[1]+"\n")
            else:
                quit(row)
        else:
            quit(row)

    EXCLUDE.close()
    BIM.close()
    LOG.close()

    fm.run_cmd("cp %s.bed %s.bed" % (args.bfile,tmp_prefix))
    fm.run_cmd("cp %s.fam %s.fam" % (args.bfile,tmp_prefix))
    fm.run_cmd("plink --bfile %s --exclude %s --make-bed --out %s" % (tmp_prefix,exclude_file,args.out))

parser = argparse.ArgumentParser(description='XXX pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bfile',help='VCF file',required=True)
parser.add_argument('--ref-dir',help='VCF file',required=True)
parser.add_argument('--out',help='VCF file',required=True)
parser.add_argument('--remove-exm',action="store_true",help='VCF file')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
