#Bash script to quantify all FASTA files in local folder with salmon 

for element in *.fastq

do

salmon quant -i /Volumes/scRNAseq_1/Homo_sapiens/homo_index -l A -r $element -o $element"_quant"

done