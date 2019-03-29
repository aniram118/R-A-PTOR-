# Find the composition of each base in the sequence per line

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def base_comp(filename):
    
    tail_data = pd.read_csv(filename, sep=" ", error_bad_lines=False,
                                names=['ReadID', 'Length', 'UMRSequence', 'Chrom', 'Start', 'End', 'Strand'])
    umr_df=tail_data[['ReadID','UMRSequence']].copy()
    umr_df['UMRSequence'] = [z.strip().replace('T', 'U') for z in umr_df['UMRSequence']]
    umr_df['A'] = ''
    umr_df['U'] = ''
    umr_df['G'] = ''
    umr_df['C'] = ''
    A=[]
    U=[]
    G=[]
    C=[]

    for i,seq in enumerate(umr_df['UMRSequence']):
        n = len(seq)
        a = 0
        u = 0
        g = 0
        c = 0
        for x in seq:
            if x == 'A':
                a+=1
            elif x == 'U':
                u+=1
            elif x == 'G':
                g+=1
            else:
                c+=1
        umr_df.at[i,'A'] = (np.around(a/n, decimals=2))
        umr_df.at[i,'U'] = (np.around(u/n, decimals=2))
        umr_df.at[i,'G'] = (np.around(g/n, decimals=2))
        umr_df.at[i,'C'] = (np.around(c/n, decimals=2))
        A.append(seq.count('A'))
        U.append(seq.count('U'))
        G.append(seq.count('G'))
        C.append(seq.count('C'))

    A=sum(A)
    U=sum(U)
    G=sum(G)
    C=sum(C)
    s= A + U + G + C
    A=A/s
    U=U/s
    G=G/s
    C=C/s

    y=[A,U,G,C]
    bars=['A','U','G','C']
    y_pos = np.arange(len(bars))

    plt.grid(zorder=5)
    plt.bar(y_pos, y, color=['red', 'green', 'blue', 'black'], width=0.7, align='center',zorder=3)
    plt.xticks(y_pos, bars)
    plt.ylabel('Percentage', fontweight='bold')
    plt.xlabel('Nucleotide', fontweight='bold')
    plt.savefig('Nucleotide_composition.pdf')
    plt.show()

    # print(umr_df)
    umr_df.to_csv(filename[:-4] + "_Nucleotidecomposition.csv",header=True, sep=',', mode='a', index=False)
    return umr_df


# new = pd.DataFrame({'ReadID': [12343,1231231],'Sequence':['AUGCUGCAUGCAUGCGAUGCAAAAAAAAA','AUGCAUGCAUGCAGCUGAGCUAG']})
# base_comp(new)