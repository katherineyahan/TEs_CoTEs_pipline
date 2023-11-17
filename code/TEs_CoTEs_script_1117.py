#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import numpy as np


def visualize(TEs,CoTEs,chromosomes,plot_name):
    bar_TEs=[]
    bar_CoTEs=[]
    total_TEs=len(TEs)
    total_CoTEs=len(CoTEs)
    for chromosome in chromosomes:
        num_TEs= len(TEs[TEs['chr'] == chromosome])
        percent_TE=num_TEs/total_TEs
        bar_TEs.append(percent_TE)
        num_CoTEs=len(CoTEs[CoTEs['chr'] == chromosome])
        percent_CoTEs=num_CoTEs/total_CoTEs
        bar_CoTEs.append(percent_CoTEs)
    bar_width = 0.4
    index = np.arange(len(chromosomes))

    fig, ax = plt.subplots(figsize=(12, 8))

    cotes_bars = ax.bar(index , bar_CoTEs, bar_width, label='CoTEs')
    tes_bars = ax.bar(index+bar_width , bar_TEs, bar_width, label='TEs')

    ax.set_xlabel('Chromosomes')
    ax.set_ylabel('Percentage')
    ax.set_title(plot_name)
    ax.set_xticks(index + bar_width / 2)
    ax.set_xticklabels(chromosomes, rotation=60, ha='right')
    ax.legend()

    #cannot include the .txt
    plot_name=plot_name+'.png'
    output_path = os.path.join('../output_1117/', plot_name)

    plt.savefig(output_path)
    plt.show()

    return bar_CoTEs,bar_TEs

#this function find the outliers of chromosomes percentage differences between TEs and CoTEs
def find_outliers_IQR(df):

   q1=df.quantile(0.25)

   q3=df.quantile(0.75)

   IQR=q3-q1

   outliers = df[((df<(q1-1.5*IQR)) | (df>(q3+1.5*IQR)))]

   return outliers

#read TEs bed files
TEs_path = sys.argv[1]
#read CoTEs bed files
CoTE_path = sys.argv[2]

TEs_df = pd.read_csv(TEs_path,delimiter='\t',header=None)
CoTE_df=pd.read_csv(CoTE_path,delimiter='\t',header=None)


# Assign columns to the two dataframes
TEs_df.columns=['chr','start','end','N']
del TEs_df['N']
CoTE_df.columns=['chr','start','end']

unique_chr=['chr1','ch2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chY']

# Extract the file name and extension from the file path
file_name, file_extension = os.path.splitext(os.path.basename(TEs_path))

#run the function to generate plots
per_CoTEs,per_TEs=visualize(TEs_df,CoTE_df,unique_chr,file_name)

#modify the dataframe to include a difference column
df = pd.DataFrame({'chr': unique_chr, 'TEs': per_TEs, 'CoTEs': per_CoTEs,'difference':[x - y for x, y in zip(per_TEs, per_CoTEs)]})

outliers=find_outliers_IQR(df['difference'])
selected_rows = df.iloc[outliers.index]

# Change the extension to '.txt'
new_file_name = file_name + '.txt'

# Concatenate with '../output/'
output_path = os.path.join('../output_1117/', new_file_name)

selected_rows.to_csv(output_path, sep='\t', index=False)




