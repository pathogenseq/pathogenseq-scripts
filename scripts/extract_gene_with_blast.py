import argparse
import fastq2matrix as fm
import json


def revcom(s):
    """Return reverse complement of a sequence"""
    def complement(s):
            basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N','-':'-'}
            letters = list(s)
            letters = [basecomplement[base] if base in basecomplement else "N" for base in letters]
            return ''.join(letters)
    return complement(s[::-1])



def parse_blast(filename,min_identity):
    result_json = json.load(open(filename))
    hsps = []
    for d in result_json["BlastOutput2"][0]["report"]["results"]["bl2seq"]:
        for hit in d["hits"]:
            for hsp in hit["hsps"]:
                h = {
                    "identity":hsp["identity"],
                    "subject_start":hsp["hit_from"],
                    "subject_end":hsp["hit_to"],
                    "subject_seq":hit["description"][0]["id"],
                    "subject_strand":hsp["hit_strand"]
                }
                if h["identity"]<min_identity: continue
                hsps.append(h)
    return hsps

def main(args):
    fm.filecheck(args.query)
    fm.filecheck(args.subject)

    ref_gene_seq = list(fm.fasta(args.query).fa_dict.values())[0]

    start_anchor = ref_gene_seq[:args.anchor_size]
    end_anchor = ref_gene_seq[-args.anchor_size:]

    tmp_in = fm.get_random_file()
    tmp_out = fm.get_random_file()
    with open(tmp_in,"w") as O:
        O.write(">tmp\n%s" %start_anchor)
    fm.run_cmd("blastn -task blastn -query %s -subject %s -outfmt 15 > %s" % (tmp_in,args.subject,tmp_out),verbose=0)
    start_hits = parse_blast(tmp_out,args.anchor_size*0.9)

    with open(tmp_in,"w") as O:
        O.write(">tmp\n%s" % end_anchor)
    fm.run_cmd("blastn -task blastn -query %s -subject %s -outfmt 15 > %s" % (tmp_in,args.subject,tmp_out),verbose=0)
    end_hits = parse_blast(tmp_out,args.anchor_size*0.9)

    fm.rm_files([tmp_in,tmp_out])

    result_type = ""
    if args.strict_one_hit and (len(start_hits)>1 or len(end_hits)>1):
        result_type = "NA"
    else:
        if start_hits[0]["subject_seq"]==end_hits[0]["subject_seq"]:
            result_type = "OK"
            start_hit = start_hits[0]
            end_hit = end_hits[0]
        else:
            result_type = "Fragmented"


    with open("%s.result.txt" % args.prefix,"w") as O:
        O.write("%s\t%s\n" % (args.prefix,result_type))

    if result_type!="OK":
        quit()

    print(start_hit,end_hit)
    subject_seqs = fm.fasta(args.subject).fa_dict
    if start_hit["subject_strand"]=="Plus" and end_hit["subject_strand"]=="Plus":
        hit_seq = subject_seqs[start_hit["subject_seq"]][start_hit["subject_start"]-1:end_hit["subject_end"]]
    elif start_hit["subject_strand"]=="Minus" and end_hit["subject_strand"]=="Minus":
        hit_seq = revcom(subject_seqs[start_hit["subject_seq"]][end_hit["subject_end"]-1:start_hit["subject_start"]])

    # import pdb; pdb.set_trace()
    with open("%s.extracted_seq.fa" % args.prefix,"w") as O:
        O.write(">%s\n%s\n" % (args.prefix,hit_seq))





parser = argparse.ArgumentParser(description='XXX pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--prefix',help='Prefix for the result files',required=True)
parser.add_argument('--query',help='Fasta with gene of interest',required=True)
parser.add_argument('--subject',help='Fasta assembly to extract from',required=True)
parser.add_argument('--anchor-size',default=40,help='VCF file')
parser.add_argument('--strict-one-hit',action="store_true",help='VCF file')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
