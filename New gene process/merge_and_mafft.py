#!/usr/bin/env python3
import os
import sys
import subprocess

if len(sys.argv)<2:
    print("Usage: merge_and_mafft.py <species>")
    sys.exit(1)

sp=sys.argv[1]

faa1=f"./spaln_results_{sp}/all.pep.fa"
faa2="Dpse_denovo.faa"
mapping_file="denovolist.txt"

outdir=f"merged_fa_{sp}"
os.makedirs(outdir,exist_ok=True)

def read_faa(f):
    d={}
    name=None
    seq=[]
    with open(f) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    d[name]="".join(seq)
                name=line.strip()[1:]
                seq=[]
            else:
                seq.append(line.strip())
        if name:
            d[name]="".join(seq)
    return d

dm=read_faa(faa1)
dp=read_faa(faa2)

with open(mapping_file) as fh:
    for line in fh:
        if not line.strip():
            continue
        c1,c2=line.strip().split()

        if c1 not in dp or c2 not in dm:
            continue

        outfile=os.path.join(outdir,f"{c2}.fa")

        with open(outfile,"w") as out:
            out.write(f">{c1}\n{dp[c1]}\n")
            out.write(f">{c2}\n{dm[c2]}\n")

        aln_file=outfile+".aln"
        subprocess.run(["mafft","--auto",outfile],stdout=open(aln_file,"w"))

print("Done.")
