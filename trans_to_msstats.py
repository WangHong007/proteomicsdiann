import pandas as pd 
import re
import os
import sys
import time

start = time.time()
report = sys.argv[1]
UNIMOD = sys.argv[2]
experiment = sys.argv[3]
rep = pd.read_csv(report, sep = "\t", header = 0, dtype = 'str')
unimod = pd.read_csv(UNIMOD, sep = ",", header = 0, dtype = 'str')

rep1 = rep[['Protein.Names', 'Modified.Sequence', 'Precursor.Charge', 'PG.Quantity', 'File.Name']]
rep1.columns = ['ProteinName', 'PeptideSequence', 'PrecursorCharge', 'Intensity', 'Reference']
col_name = rep1.columns.tolist() 
add_columns = ['FragmentIon', 'ProductCharge', 'IsotopeLabelType', 'Condition', 'BioReplicate', 'Run']
for i in range(3, 9):
    col_name.insert(i, add_columns[i-3])
rep1 = rep1.reindex(columns = col_name)
for i in range(0, len(rep1)):
    rep1.loc[i, 'Reference'] = rep1.loc[i, 'Reference'].split('/')[-1] 
    rep1.loc[i, 'FragmentIon'] = 'NA'
    rep1.loc[i, 'ProductCharge'] = '0'
    rep1.loc[i, 'IsotopeLabelType'] = 'L'
    PTM = re.findall(re.compile(r'[(](.*?)[)]', re.S), rep1.loc[i, 'PeptideSequence'])               
    if len(PTM) > 0:
        for j in range(0, len(unimod)):
            for m in range(0, len(PTM)):
                if(PTM[m] == unimod.loc[j, 'id']):
                    rep1.loc[i, 'PeptideSequence'] = re.sub(r'(?<=\().+?(?=\))', unimod.loc[j, 'name'], rep1.loc[i, 'PeptideSequence'])
                if(rep1.loc[i, 'PeptideSequence'][0] == '('):
                    rep1.loc[i, 'PeptideSequence'] = '.' + rep1.loc[i, 'PeptideSequence']

    
for i in range(0, len(rep1)):
    file_extension = 'file_extension'
    while(file_extension != ''):
        file_name, file_extension = os.path.splitext(rep1.loc[i, 'Reference'])
        rep1.loc[i, 'Reference'] = file_name
    rep1.loc[i, 'Reference'] = rep1.loc[i, 'Reference'] + '.mzML'

##
## 
##
# experiment = "C:\\Users\\Mr.HG\\Desktop\\expr.tsv"
exper = pd.read_csv(experiment, sep = "\t", header = 0, dtype = 'str')
for i in range(0, len(exper)):
    if(exper.loc[i, 'Fraction_Group'] == 'Sample'):
        len_df1 = i
        len_df2 = len(exper) - len_df1
df1 = exper.iloc[:len_df1,:]
df2 = exper.iloc[len_df1:,:]


df1_col_name = df1.columns.tolist() 
df1_add_columns = ['MSstats_Condition', 'MSstats_BioReplicate']
for i in range(5, 7):
    df1_col_name.insert(i, df1_add_columns[i-5])
df1 = df1.reindex(columns = df1_col_name)
df1 = df1.dropna(axis=0, subset=['Label'])
for i in range(0, len(df1)):
    for j in range(len_df1, len(exper)):
        if(df2.loc[j, 'Fraction_Group'] == df1.loc[i, 'Sample']):
            df1.loc[i, 'MSstats_Condition'] = df2.loc[j, 'Fraction']
            df1.loc[i, 'MSstats_BioReplicate'] = df2.loc[j, 'Spectra_Filepath']
    file_extension = 'file_extension'
    while(file_extension != ''):
        file_name, file_extension = os.path.splitext(df1.loc[i, 'Spectra_Filepath'])
        df1.loc[i, 'Spectra_Filepath'] = file_name
    df1.loc[i, 'Spectra_Filepath'] = df1.loc[i, 'Spectra_Filepath'] + '.mzML'


##
## 
##
for i in range(0, len(rep1)):
    for j in range(0, len(df1)):
        if df1.loc[j, 'Spectra_Filepath'] == rep1.loc[i, 'Reference']:
            rep1.loc[i, 'Condition'] = df1.loc[j, 'MSstats_Condition']
            rep1.loc[i, 'BioReplicate'] = df1.loc[j, 'MSstats_BioReplicate']
            rep1.loc[i, 'Run'] = df1.loc[j, 'Fraction_Group']


rep1 = rep1.dropna(axis=0, subset=['ProteinName'])

rep1.to_csv(path_or_buf='./MSstats.csv', sep=',', index=False)
print(rep1)

end = time.time()
print("Run time:{t}".format(t=end - start))
