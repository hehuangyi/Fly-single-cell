{\rtf1\ansi\ansicpg936\cocoartf2707
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/usr/bin/env python3\
\
"""\
Pipeline to identify de novo genes.\
\
Final output = intersection of:\
    v1: BLAST-based filters\
    v2: Spaln-based filters\
\
No ID conversion.\
"""\
\
import os\
import glob\
import pandas as pd\
\
\
#############################\
# Parameters (edit here)\
#############################\
\
BASE = "/Volumes/Elements/hhy/project/03.testis single cell new/gene list/Orthofinder2/denovo_blast"\
\
gene_length_file = f"\{BASE\}/denovo_Length.txt"\
\
spaln_folder = f"\{BASE\}/spaln result"\
merged_gene_file = f"\{spaln_folder\}/merged_uniq_gene.txt"\
\
output_file = f"\{BASE\}/Dpse_denovo_final.txt"\
\
# species requiring BOTH blastn & blastp absence\
species = [\
    "dmel","dwil","dyak","dana","dere",\
    "dgri","dmoj","dsim","dsuz"\
]\
\
\
#############################\
# Functions\
#############################\
\
def read_gene_length():\
    gl = pd.read_table(gene_length_file, header=None)\
    gl.columns = ["gene", "length"]\
    return gl\
\
\
def blast_filter(file, gl):\
    df = pd.read_table(file, header=None)\
    df = pd.merge(df, gl, left_on=0, right_on="gene", how="right")\
    df = df.fillna(0)\
\
    df[3] = pd.to_numeric(df[3])\
    df[2] = pd.to_numeric(df[2])\
\
    df["ratio"] = df[3] / df["length"]\
\
    # keep BAD hits (short + low identity)\
    df = df[(df["ratio"] < 0.30) & (df[2] < 30)]\
\
    return set(df["gene"])\
\
\
#############################\
# Step 1. v1 from BLAST\
#############################\
\
print("Running BLAST filters...")\
\
gl = read_gene_length()\
\
# start from all genes\
v1 = set(gl["gene"])\
\
for sp in species:\
    print("  ", sp)\
\
    cds = f"\{BASE\}/Dpse.denovo.blastn_to\{sp\}.uniq"\
    prot = f"\{BASE\}/Dpse.denovo.blastp_to\{sp\}.uniq"\
\
    cds_set = blast_filter(cds, gl)\
    prot_set = blast_filter(prot, gl)\
\
    v1 = v1.intersection(cds_set)\
    v1 = v1.intersection(prot_set)\
\
\
# uniprot\
uniprot = f"\{BASE\}/Dpse.denovo.blastp_tounip.uniq"\
uniprot_set = blast_filter(uniprot, gl)\
v1 = v1.intersection(uniprot_set)\
\
\
# remove good self hits\
print("  removing Dpse strong matches")\
dpse = pd.read_table(f"\{BASE\}/Dpse.denovo.blastp_todpse.uniq", header=None)\
dpse = pd.merge(dpse, gl, left_on=0, right_on="gene", how="left").fillna(0)\
dpse[3] = pd.to_numeric(dpse[3])\
dpse[2] = pd.to_numeric(dpse[2])\
dpse["ratio"] = dpse[3] / dpse["length"]\
\
good_self = set(dpse[(dpse[2] > 30) & (dpse["ratio"] > 0.30)][0])\
\
v1 = v1 - good_self\
\
print("v1 size =", len(v1))\
\
\
#############################\
# Step 2. v2 from Spaln\
#############################\
\
print("Running Spaln filters...")\
\
files = glob.glob(os.path.join(spaln_folder, "denovo_summary.*"))\
\
conserved = set()\
\
for f in files:\
    df = pd.read_table(f)\
    genes = df.loc[df["Description"] == "Conserved", "Gene"]\
    conserved.update(genes)\
\
raw = pd.read_table(merged_gene_file, header=None)[0]\
\
v2 = set(raw) - conserved\
\
print("v2 size =", len(v2))\
\
\
#############################\
# Step 3. Final\
#############################\
\
final = sorted(v1.intersection(v2))\
\
print("Final size =", len(final))\
\
with open(output_file, "w") as f:\
    for g in final:\
        f.write(g + "\\n")\
\
print("Done.")\
print("Output:", output_file)\
}