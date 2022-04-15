import argparse
import re
parser = argparse.ArgumentParser()
parser.add_argument('--partition',action='store',default=1,required=False)
args = parser.parse_args()
partition = args.partition
fname = 'CCAPart%02d.csv.gz'% partition
outfname = 'locusZoomPart%02d.csv'% partition
import pandas as pd
df = pd.read_csv(fname)
df.rename(columns={
    'chromosome': 'Chromosome',
    'P_value':'p-value',
    'A1':'Ref allele',
    'A2':'Alt allele',
    'position':'Position'
},inplace=True)
df = df[['Chromosome','p-value', 'Ref allele', 'Alt allele', 'Position']]
df.sort_values(['Chromosome','Position'],inplace=True)
df.to_csv(outfname, index=False,sep='\t')
import gzip
with open(outfname, 'r') as inp:
    with gzip.open(outfname + '.gz', 'wb') as output:
        output.write(inp.read())
import os
os.remove(outfname)