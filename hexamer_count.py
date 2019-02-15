

import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt

def hexamerplot(dataframe):
    #  This function generates plots for Hexamer frequency

    data1 = dataframe[['ReadID', 'Sequence']].copy()
    data1['Sequence'] = [z.strip().replace('T', 'U') for z in data1['Sequence']]
    data1['AAUAAA']=''
    data1['AUUAAA']=''
    data1['AGUAAA']=''
    data1['UAUAAA']=''
    data1['CAUAAA']=''
    data1['GAUAAA']=''
    data1['AAUAUA']=''
    data1['AAUACA']=''
    data1['AAUAGA']=''
    data1['ACUAAA']=''
    data1['AAGAAA']=''
    data1['AAUGAA']=''
    dlen=len(data1['Sequence'])

    for i,val in enumerate(data1['Sequence']):
        data1.at[i,'AAUAAA']=int(val.count("AAUAAA"))
        data1.at[i,'AUUAAA']=int(val.count("AUUAAA"))
        data1.at[i,'AGUAAA']=int(val.count("AGUAAA"))
        data1.at[i,'UAUAAA']=int(val.count("UAUAAA"))
        data1.at[i,'CAUAAA']=int(val.count("CAUAAA"))
        data1.at[i,'GAUAAA']=int(val.count("GAUAAA"))
        data1.at[i,'AAUAUA']=int(val.count("AAUAUA"))
        data1.at[i,'AAUACA']=int(val.count("AAUACA"))
        data1.at[i,'AAUAGA']=int(val.count("AAUAGA"))
        data1.at[i,'ACUAAA']=int(val.count("ACUAAA"))
        data1.at[i,'AAGAAA']=int(val.count("AAGAAA"))
        data1.at[i,'AAUGAA']=int(val.count("AAUGAA"))
    # print(data1)
    data1.to_csv("hexamer_composition.csv",header=True,mode='a')

    s1=data1['AAUAAA'].sum()
    s2=data1['AUUAAA'].sum()
    s3=data1['AGUAAA'].sum()
    s4=data1['UAUAAA'].sum()
    s5=data1['CAUAAA'].sum()
    s6=data1['GAUAAA'].sum()
    s7=data1['AAUAUA'].sum()
    s8=data1['AAUACA'].sum()
    s9=data1['AAUAGA'].sum()
    s10=data1['ACUAAA'].sum()
    s11=data1['AAGAAA'].sum()
    s12=data1['AAUGAA'].sum()
    s=[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12]

    # print(s)
    ypos=np.arange(12)
    signals=['AAUAAA','AUUAAA','AGUAAA','UAUAAA','CAUAAA','GAUAAA','AAUAUA','AAUACA','AAUAGA','ACUAAA','AAGAAA','AAUGAA']
    # val=[]
    # for hexamer in s:
    #         val.append(sum(hexamer)/dlen)

    plt.bar(ypos,s,align='center',alpha=0.5)
    plt.xticks(ypos,signals, rotation=38.5)
    # plt.show()
    plt.savefig('hexamer_distribution.jpg')

    return data1


# dict=pd.DataFrame({'ReadID':[1231412,1231123,1243124123],'Sequence':['AAUAAAAUGCAUGCAUGC','AAUGAAGUAGUCGAGUCG','AAGAAAAUCGCUAACUAAAAAUGAA']})
# # sampledf=pd.DataFrame(dict)
# hexamerplot(dict)

