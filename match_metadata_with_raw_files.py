import pandas as pd
import os
from os import path
#import seaborn as sns

#from IPython.display import display


DATA_PATH=r'Z:\data\gqc954\Natera_Transfer'


samplesheet_pathlist=[r'metadata\spectrum_meta_data_2018.xlsx',
                      r'metadata\spectrum_meta_data_2019.xlsx',
                      r'metadata\spectrum_meta_data_2020.xlsx'
                     ]

names = [ '2018',
         '2019',
         '2020'
]
samplesheets=[pd.read_excel(os.path.join(DATA_PATH, i),engine="openpyxl") for i in samplesheet_pathlist]

#all_samplesheets=pd.concat(samplesheets)

#for g,df in all_samplesheets.groupby('casefile_id'):
#print(g,len(df))

for n,s in zip(names,samplesheets):
    for ar in  s['array']:
        if path.exists(os.path.join(DATA_PATH, n,ar+'.xy.gz')) and path.exists(os.path.join(DATA_PATH, n, ar + '.xy.gz')):
            continue
        else:
            raise Exception('problem with ' + ar)



#print(samplesheet.describe(include='all'))