#Python script to export Smartseq2 count.tab
import pandas as pd

df1 = pd.read_csv('/Volumes/scRNAseq_1/SS2_15_0150-0151/patient_1_layout.csv')
Well_number_list = list(df1['Well number'])
Well_number_list.insert(0, "gene")

df = pd.read_csv("/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/counts.tab", sep='\t', lineterminator='\n')
df2 = pd.DataFrame(df[Well_number_list])
df2.to_csv("/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/patient1.csv", index=False, sep = "\t")