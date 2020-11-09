import sys
from ete3 import Tree
import subprocess
from collections import defaultdict
import argparse
import functools


def main(args):

    def is_cluster(n,num_isolate_cutoff=10,max_dist_cutoff=0.01):
        num_isolates = len(n.get_leaf_names())
        max_dist = n.get_farthest_leaf()[1]
        if num_isolates>num_isolate_cutoff and max_dist<max_dist_cutoff:
            return True
        else:
            return False

    is_cluster = functools.partial(is_cluster,num_isolate_cutoff=args.num_samples,max_dist_cutoff=args.max_dist)

    clusters = []
    t = Tree(args.tree)
    for n in t.iter_leaves(is_leaf_fn=is_cluster):
        print(n)
        clusters.append(n.get_leaf_names())
    with open(args.output,"w") as O:
        for i in range(len(clusters)):
            for s in clusters[i]:
                O.write("%s\t%s\n" % (s,i))

    if args.itol:
        for i in range(len(clusters)):
            with open(f"{args.output}.{i}.itol.txt","w") as O:
                O.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tlabel1\nCOLOR\t#ff0000\nDATA\n")
                for s in clusters[i]:
                    O.write("%s\t%s\n" % (s,"black"))


parser = argparse.ArgumentParser(description='VCF mergin pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--tree',help='sample file',required=True)
parser.add_argument('--max-dist',type=float,help='sample file',required=True)
parser.add_argument('--num-samples',type=int,help='sample file',required=True)
parser.add_argument('--output',type=str,help='sample file',required=True)
parser.add_argument('--itol',action="store_true")
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
