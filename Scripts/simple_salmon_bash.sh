#Bash script to quantify all FASTA files in local folder with salmon 

mkkdir Quant

for element in *.fastq

do

    var1=$element
    replace=""
    var2=${var1/.fastq/$replace}
    
    echo $var2
    

salmon quant -i /Volumes/scRNAseq_1/Homo_sapiens/homo_index -l A -r $element -o Quant/$var2

done