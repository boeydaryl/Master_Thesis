#Salmon mass quantification bash script

cd /Volumes/scRNAseq_1/SS2_15_0150-0151/raw_fastq_test

for element in *.fastq.gz

do

salmon quant -i /Volumes/scRNAseq_1/Homo_sapiens/homo_index -l A -r $element -o /Volumes/scRNAseq_1/SS2_15_0150-0151/raw_fastq_test/$element"_quant"

done