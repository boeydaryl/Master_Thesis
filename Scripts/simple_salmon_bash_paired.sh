#Bash script to quantify all FASTA files in local folder with salmon with paired end libraries

for element in *_1.fastq

do

var1=$element
replace=""
var2=${var1/_1.fastq/$replace}

echo $var2

salmon quant -i /Volumes/scRNAseq_1/Homo_sapiens/homo_index -l IU -1 $var2"_1.fastq" -2 $var2"_2.fastq"  -o $var2"_quant"

done