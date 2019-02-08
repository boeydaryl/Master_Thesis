#Python script to export Smartseq2 count.tab
import pandas as pd
import os

path = "/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/Separated_counts/"

df = pd.read_csv("/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/counts.tab", sep='\t', lineterminator='\n')
for column in df:
    #print (df[column])
    output_name = path + column
    os.mkdir(output_name)
    #print(output_name)
    output_file_name = output_name + "/" + "quant.csv"
    #print(output_file_name)
    header = ["gene", column]
    df.to_csv(output_file_name, columns = header)