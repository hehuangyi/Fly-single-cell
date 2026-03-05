{\rtf1\ansi\ansicpg936\cocoartf2707
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 cellranger mkref --genome=cellranger_Dpse --fasta=Dpse.UCI_MV25_Y.fa --genes=Dpse.Convert.gtf\
cellranger count --id=Dpse_larvae_1 --fastqs=Dpse_larvae_1/ --chemistry=SC3Pv2 --transcriptome=cellranger_Dpse --localcores=20 --localmem=500 > log.Dpse_larvae1 &\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 cellranger count --id=Dpse_larvae_2 --fastqs=Dpse_larvae_2/ --chemistry=SC3Pv2 --transcriptome=cellranger_Dpse --localcores=20 --localmem=500 > log.Dpse_larvae2 &\
cellranger count --id=Dpse_adult_1 --fastqs=Dpse_adult_1/ --chemistry=SC3Pv2 --transcriptome=cellranger_Dpse --localcores=20 --localmem=500 > log.Dpse_adult1 &\
cellranger count --id=Dpse_adult_2 --fastqs=Dpse_adult_2/ --chemistry=SC3Pv2 --transcriptome=cellranger_Dpse --localcores=20 --localmem=500 > log.Dpse_adult2 &}