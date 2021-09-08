'''
Script takes metadata and gzip locations and creates dataframe per family
'''
import pandas as pd
import os
from os import path
import re
from  main_data import process_xy_file,process_ballele_file
import numpy as np
from main_data import SNP_MAP
from sys import platform
import multiprocessing
from joblib import Parallel, delayed
from multiprocessing import Process, Manager


'''
ind is familyid
'''

def main(ind,batchid,samplesheet):
    #get arrray list
    arlist=samplesheet.loc[samplesheet.casefile_id == ind, 'array'].values
    namelist=samplesheet.loc[samplesheet.casefile_id == ind, ['casefile_id','family_position']].apply(
          lambda x: str(x[0]) + '_' + str(x[1]), axis=1).values

    ar=[]
    #for xydata,bdata in zip(sorted(xydatalist),sorted(bdatalist)):
    for aritem,nameitem in zip(arlist,namelist):
      xydata=os.path.join(DATA_PATH,batchid,aritem + '.xy.gz')
      bdata=os.path.join(DATA_PATH,batchid,aritem + '.b.gz')
      samplename=nameitem
      #print(ind,xysamplename,samplename)
      status_xy,xy=process_xy_file(xydata)
      status_baf,baf=process_ballele_file(bdata)
      if status_xy==False:
          q.put('{0},{1},{2},{3},{4}'.format(str(os.getpid()),str(batchid),str(ind),str(xydata),'empty'))
          #print('{0},{1},{2},{3}'.format(str(batchid),str(ind),str(xydata),'empty'))
      if status_baf==False:
          q.put('{0},{1},{2},{3},{4}'.format(str(os.getpid()),str(batchid), str(ind), str(bdata), 'empty'))
          #print('{0},{1},{2},{3}'.format(str(batchid), str(ind), str(bdata), 'empty'))
      df=pd.DataFrame(np.concatenate(xy), columns=[(samplename,'x_raw'), (samplename,'y_raw')])
      df[(samplename,'baf')] = np.concatenate(baf)
      ar.append(df)
    output=pd.concat(ar,axis=1)
    output.columns=pd.MultiIndex.from_tuples(output.columns,names=['individual','feature'])
    #print(output.columns)
    final=pd.concat([output,pd.DataFrame(np.concatenate(SNP_MAP))],axis=1).set_index(['position','rsid'])
    final.to_pickle(os.path.join(OUTPUT_PATH,str(batchid),str(ind)+ '.pkl'))
    return(ind)

def worker(familylist,batchid,samplesheet):
    for family in familylist:
        main(family,batchid,samplesheet)
    return(len(familylist))

def save_to_file(q):
    with open('empty_files.txt', 'w') as out:
        while True:
            val = q.get()
            if val is None: break
            out.write(val + '\n')


if platform == "linux" or platform == "linux2":
    DATA_PATH=r'/data/gqc954/Natera_Transfer'
    OUTPUT_PATH = r'/data/gqc954/Natera_Transfer/Serialized'
else:
    DATA_PATH=r'Z:\data\gqc954\Natera_Transfer'
    OUTPUT_PATH=r'Z:\data\gqc954\Natera_Transfer\Serialized'


samplesheet_pathlist=[os.path.join(r'metadata',r'spectrum_meta_data_2018.xlsx'),
                      os.path.join(r'metadata',r'spectrum_meta_data_2019.xlsx'),
                      os.path.join(r'metadata',r'spectrum_meta_data_2020.xlsx')
                     ]

names = [ '2018',
         '2019',
         '2020'
]
if platform!='linux' and platform!='linux2':
  samplesheets=[pd.read_excel(os.path.join(DATA_PATH, i),engine='openpyxl') for i in samplesheet_pathlist]
else:
    samplesheets = [pd.read_excel(os.path.join(DATA_PATH, i)) for i in samplesheet_pathlist]


#fest 1 family

#g=1591706
#num_cores = multiprocessing.cpu_count()
num_cores=22


#main(g,'2018',samplesheets[0])


#g=2785298
#main(g,'2020',samplesheets[2])
m = Manager()
q = m.Queue()
p = Process(target=save_to_file, args=(q,))
p.start()
####
for batchid,samplesheet in zip(names,samplesheets):
    #get nr of families
    nrfamilies=len(set(samplesheet.casefile_id))
    Parallel(n_jobs=num_cores)(delayed(main)(i, batchid,samplesheet) for i in set(samplesheet.casefile_id))
q.put(None)
p.join()

#s_g=samplesheets[0][samplesheets[0].casefile_id==g]
#get list of arrays


'''
for n,s in zip(names,samplesheets):
    nr_of_families=len(set(s.casefile_id))
    #divide families into 24 CP

    for g,s_g in samplesheets.groupby('casefile_id'):
        s_g['array']
        if path.exists(os.path.join(DATA_PATH, n,ar+'.xy.gz')) and path.exists(os.path.join(DATA_PATH, n, ar + '.xy.gz')):
            continue
        else:
            raise Exception('problem with ' + ar)
'''