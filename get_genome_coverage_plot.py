import sys
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import subprocess
import argparse
import math
from collections import defaultdict
from tqdm import tqdm

def main(args):
    x = defaultdict(list)
    y = defaultdict(list)
    num_windows = None
    for l in tqdm(subprocess.Popen("bedtools makewindows -w %(window_size)s -g %(ref)s.fai | wc -l" % vars(args),shell=True, stdout=subprocess.PIPE).stdout):
        num_windows = int(l.decode().strip())
    for l in tqdm(subprocess.Popen("bedtools makewindows -w %(window_size)s -g %(ref)s.fai | bedtools coverage -mean -a - -b %(bam)s -sorted -nobuf" % vars(args),shell=True, stdout=subprocess.PIPE).stdout,total=num_windows):
        row = l.decode().strip().split()
        x[row[0]].append((int(row[1]) + int(row[2]))/2)
        y[row[0]].append(float(row[3]))

    fig = make_subplots(rows=len(x), cols=1, subplot_titles=list(x))
    for i,chrom in enumerate(x):
        fig.append_trace(go.Scatter(go.Scatter(x=x[chrom], y=y[chrom])), row=i+1, col=1)

    fig.update_layout(height=600*len(x), title_text="Coverage across chromosomes in %s" % args.bam)
    fig.write_html(args.out)

parser = argparse.ArgumentParser(description='Get coverage plots',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam',help='Bam file',required=True)
parser.add_argument('--ref',help='Reference file',required=True)
parser.add_argument('--out',default="coverage_plot.html",help='Output file name')
parser.add_argument('--window-size',type=int,default=1000,help='Window size')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
