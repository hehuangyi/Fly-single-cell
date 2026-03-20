cellranger mkref --genome=cellranger_Dpse --fasta=Dpse.UCI_MV25_Y.fa --genes=Dpse.Convert.gtf
cellranger count --id=Dpse_larvae_1 --fastqs=Dpse_larvae_1/ --chemistry=SC3Pv2 --transcriptome=cellranger_Dpse --localcores=20 --localmem=500 > log.Dpse_larvae1 &
cellranger count --id=Dpse_larvae_2 --fastqs=Dpse_larvae_2/ --chemistry=SC3Pv2 --transcriptome=cellranger_Dpse --localcores=20 --localmem=500 > log.Dpse_larvae2 &
cellranger count --id=Dpse_adult_1 --fastqs=Dpse_adult_1/ --chemistry=SC3Pv2 --transcriptome=cellranger_Dpse --localcores=20 --localmem=500 > log.Dpse_adult1 &
cellranger count --id=Dpse_adult_2 --fastqs=Dpse_adult_2/ --chemistry=SC3Pv2 --transcriptome=cellranger_Dpse --localcores=20 --localmem=500 > log.Dpse_adult2 &
