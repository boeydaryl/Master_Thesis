#!/usr/bin/env bash
#Bash Script for Salmon pipeline quantification
echo Important! Remember to activate salmon with: conda activate salmon
cd /Volumes/Jetstream/Thesis_Raw_Data
echo Input path to FASTQ files
read vardir1
cd $vardir1
echo Input read1 FASTQ
read varfq1
echo Input read2 FASTQ file
read varfq2
#echo $varfq1
echo Input path to save quantification folder
read vardir
cd vardir
echo Beginning Salmon Run
quant="_quant"
output="$varfq1$quant"
#echo $output
salmon quant -i /Volumes/Jetstream/Thesis_Raw_Data/Homo_sapiens/homo_index/ -l IU -1 $varfq1 -2 $varfq2 -o output
