# Find the composition of each base in the sequence per line

import pandas as pd
import numpy as np


def base_comp(dataframe):
    data1=dataframe[['ReadID','Sequence']].copy()
    data1['Sequence'] = [z.strip().replace('T', 'U') for z in data1['Sequence']]
    data1['A'] = ''
    data1['U'] = ''
    data1['G'] = ''
    data1['C'] = ''
    
    for i,seq in enumerate(data1['Sequence']):
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
        data1.at[i,'A'] = (np.around(a/n, decimals=2))
        data1.at[i,'U'] = (np.around(u/n, decimals=2))
        data1.at[i,'G'] = (np.around(g/n, decimals=2))
        data1.at[i,'C'] = (np.around(c/n, decimals=2))
        # A.append(np.around(a/n,decimals=2))
        # U.append(np.around(u/n,decimals=2))
        # G.append(np.around(g/n,decimals=2))
        # C.append(np.around(c/n,decimals=2))


    # print(data1)
    data1.to_csv("base_composition.csv",header=True, sep=',', mode='a')
    return data1


# new = pd.DataFrame({'ReadID': [12343,1231231],'Sequence':['AUGCUGCAUGCAUGCGAUGCAAAAAAAAA','AUGCAUGCAUGCAGCUGAGCUAG']})
# base_comp(new)