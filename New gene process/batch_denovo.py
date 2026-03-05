{\rtf1\ansi\ansicpg936\cocoartf2707
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fnil\fcharset0 LucidaGrande;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/usr/bin/env python3\
import os\
import sys\
import numpy as np\
import scipy.stats\
\
# =======================\
# BLOSUM\
# =======================\
blosum62 = \{('A','A'):4, ('A','R'):-1, ('A','N'):-2\}\
\
# =======================\
# FASTA parser\
# =======================\
def fasta2seq(fafile):\
    seqs = \{\}\
    with open(fafile) as f:\
        key = None\
        for line in f:\
            line=line.strip()\
            if line.startswith(">"):\
                key=line[1:]\
                seqs[key]=""\
            elif key:\
                seqs[key]+=line\
    return seqs\
\
# =======================\
# similarity\
# =======================\
def similarity_seq(seq1, seq2, prolen):\
    if len(seq1)!=len(seq2):\
        return -1\
    nsim=0\
    for s1,s2 in zip(seq1,seq2):\
        if s1!="-" and s2!="-":\
            if s1==s2:\
                nsim+=1\
            else:\
                nsim+=0\
    return min(1, nsim/prolen)\
\
# =======================\
# probability\
# =======================\
def probability_denovo(score):\
    return 1-score   \
\
# =======================\
# single aln\
# =======================\
def process_aln(file_path):\
    seqs=fasta2seq(file_path)\
    keys=list(seqs.keys())\
    if len(keys)<2:\
        return None\
\
    ref_seq=seqs[keys[0]]\
    target_seq=seqs[keys[1]]\
    prolen=len(ref_seq.replace("-",""))\
\
    score=sum(1 for a,b in zip(ref_seq,target_seq) if a==b)/prolen\
    sim=similarity_seq(ref_seq,target_seq,prolen)\
    prob=probability_denovo(score)\
\
    description="DeNovo" if prob>0.5 else "Conserved"\
\
    return [keys[1],score,sim,prob,prolen,description]\
\
# =======================\
# batch\
# =======================\
def main(input_dir, output_dir):\
    os.makedirs(output_dir,exist_ok=True)\
    out_file=os.path.join(output_dir,"denovo_summary.txt")\
\
    with open(out_file,"w") as fo:\
        fo.write("Gene\\tScore\\tSimilarity\\tP_denovo\\tLength\\tDescription\\n")\
        for f in os.listdir(input_dir):\
            if f.endswith(".aln"):\
                res=process_aln(os.path.join(input_dir,f))\
                if res:\
                    fo.write("\\t".join(map(str,res))+"\\n")\
\
    print("Finished 
\f1 \uc0\u8594 
\f0 ",out_file)\
\
if __name__=="__main__":\
    if len(sys.argv)<2:\
        print("Usage: batch_denovo.py <aln_dir> [output]")\
        sys.exit(1)\
\
    input_dir=sys.argv[1]\
    output=sys.argv[2] if len(sys.argv)>2 else "./denovo_results"\
    main(input_dir,output)\
}