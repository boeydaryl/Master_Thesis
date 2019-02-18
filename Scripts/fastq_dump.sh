# Batch download of SRA FASTQ files

cd /Volumes/scRNAseq_1/GSE66053

file_name=$(<GSE66053_to_retrieve.txt)

for element in $file_name; do
    echo $element
    /Volumes/JetStream/sratoolkit.2.9.2-mac64/bin/./fastq-dump -I --split-files -O /Volumes/scRNAseq_1/GSE66053/ $element

done

    #/Volumes/JetStream/sratoolkit.2.9.2-mac64/bin/./fastq-dump -I --split-files -O /Volumes/scRNAseq_1/GSE66053/ $element
