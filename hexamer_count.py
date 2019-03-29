

import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib as mpl
# if os.environ.get('DISPLAY','') == '':
#     print('Using non-interactive Agg backend')
#     mpl.use('Agg')
# mpl.use('TkAgg')
import matplotlib.pyplot as plt
import operator


def count_kmers(reads, k):
    counts = {}
    for read in reads:
        num_kmers = len(read) - k + 1
        for i in range(num_kmers):
            # Slice the string to get the kmer
                kmer = read[i:i+k]
                # Add the kmer to the dictionary if it's not there
                if kmer not in counts:
                    counts[kmer] = 0
                # Increment the count for this kmer
                counts[kmer] += 1
    return counts

def plothexamers(filename):
    
    tail_data = pd.read_csv(filename, sep=" ", error_bad_lines=False,
                                names=['ReadID', 'Length', 'Sequence', 'Chrom', 'Start', 'End', 'Strand'])
    tail_data['Sequence'] = [z.strip().replace('T', 'U') for z in tail_data['Sequence']]
    seqlist = list(tail_data['Sequence'])
    func = count_kmers(seqlist, 6)
    sorted_x = sorted(func.items(), key=operator.itemgetter(1), reverse=True)
    hexamers_df = pd.DataFrame.from_dict(sorted_x)
    hexamers_df.columns = ['Hexamers', 'Frequency']
    hexamers_df['Percentage'] = ''
    total = hexamers_df['Frequency'].sum()

    for i, freq in enumerate(hexamers_df['Frequency']):
        hexamers_df.at[i, 'Percentage'] = freq / total * 100
    hexamers_df.to_csv(filename[:-4] + "_hexamer_analysis.csv",mode='a')
    tophexamers = list(hexamers_df['Hexamers'][0:15])
    freq = list(hexamers_df['Percentage'][0:15])
    sns.set(style="darkgrid")
    ax = sns.barplot(x=tophexamers, y=freq)
    plt.xticks(rotation=36)
    plt.xlabel('Hexamers', fontweight='bold')
    plt.ylabel('Frequency(%)', fontweight='bold')
    plt.tight_layout()
    plt.savefig(filename[:-4] + "_top_Hexamers.pdf")




# dict=pd.DataFrame({'ReadID':[1231412,1231123,1243124123],'Sequence':['AAUAAAAUGCAUGCAUGC','AAUGAAGUAGUCGAGUCG','AAGAAAAUCGCUAACUAAAAAUGAA']})
# # sampledf=pd.DataFrame(dict)
# hexamerplot(dict)

