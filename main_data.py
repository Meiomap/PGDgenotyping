import os
import itertools
import numpy as np
import glob
import pandas as pd
import re
import gzip
from sys import platform

# Change path to point to example code directory
#EXAMPLE_CODE_PATH = '/mnt/ntapstor1/vol/spectrum_phasing/example_code/example_code'
if platform =='win32':
  SNP_MAP_FILE = r'C:\Users\gqc954\Documents\__WORK__\natera_trial\pilot_data\snp_map_cyto12b_f004.txt'
else:
    SNP_MAP_FILE = r'/data/gqc954/Natera_Transfer/pilot_data/example_code/snp_map_cyto12b_f004.txt'



def load_snp_map(filename):
    with open(filename, 'r') as fin:
        data = [line.rstrip('\n').split('\t') for line in fin][1:]

    snp_map = []
    for _, data_iter in itertools.groupby(data, key=lambda x: int(x[0])):

        snp_map.append(np.array([(pos, rsid) for _, pos, rsid in data_iter],
                                dtype=np.dtype([('position', 'int'), ('rsid', 'U32')])))

    return snp_map



def process_xy_file(array_file):
    status=True
    with gzip.open(array_file, 'rb') as fin:
        d=fin.read()
        if len(d)==0:
            status=False
        line = d.rstrip(b' ')

    array = np.fromstring(line, dtype=float, sep=' ')
    array = array.reshape([2, int(array.size / 2 / 25), 25], order='F')
    array = np.moveaxis(array, 1, 0)

    return status,[array[:snp_map.size, :, index] for index, snp_map in enumerate(SNP_MAP)]


def process_ballele_file(array_file):
    status=True
    with gzip.open(array_file, 'rb') as fin:
        if (len(fin.readlines())==0):
            array=np.full((23740, 25), np.nan)
            status=False

        #d=fin.read()
        #d_dec=d.decode()
        #array = np.vstack([np.fromstring(line.rstrip(' '), dtype=float, sep=' ') for line in d_dec]).T
        else:
            fin.seek(0)
            array = np.vstack([np.fromstring(a.decode(),dtype=float,sep=' ') for a in fin.readlines()]).T

    return (status,[array[:snp_map.size, index] for index, snp_map in enumerate(SNP_MAP)])

SNP_MAP = load_snp_map(SNP_MAP_FILE)


if __name__== "__main__":
    DATA_PATH = r'Z:\data\gqc954\Natera_Transfer'

    ARRAY_PATH = os.path.join(DATA_PATH, '2018')
    samplesheet=pd.read_excel(os.path.join(DATA_PATH,r'metadata\spectrum_meta_data_2018.xlsx'),engine="openpyxl")

    ar=[]

    for xydata,bdata in zip(glob.glob(ARRAY_PATH+ '/*.xy.gz'),glob.glob(ARRAY_PATH+ '/*.b.gz')):
        xysamplename=re.search('([0-9]+_R[0-9]{2}C[0-9]{2})', xydata).group(1)
        bdatasamplename=re.search('([0-9]+_R[0-9]{2}C[0-9]{2})', bdata).group(1)
        assert(xysamplename==bdatasamplename)
        sample_sel=samplesheet.loc[samplesheet['array']==xysamplename,['casefile_id','family_position']].apply(lambda x: str(x[0])+'_' +  str(x[1]),axis=1)
        assert(sample_sel.empty!=True)
        samplename=sample_sel.iloc[0]
        print(xysamplename,samplename)
        xy=process_xy_file(xydata)
        baf=process_ballele_file(bdata)
        df=pd.DataFrame(np.concatenate(xy), columns=[(samplename,'x_raw'), (samplename,'y_raw')])
        df[(samplename,'baf')] = np.concatenate(baf)
        ar.append(df)

    output=pd.concat(ar,axis=1)

    output.columns=pd.MultiIndex.from_tuples(output.columns,names=['individual','feature'])
    print(output.columns)
    pd.concat([output,pd.DataFrame(np.concatenate(SNP_MAP))],axis=1)\
        .set_index(['position','rsid'])\
        .to_csv('Natera_xybaf_with_index.txt')
    #output.to_csv('Natera_xybaf_dump.txt')

        #with gzip.open(i, 'rb') as f:
        #    file_content = f.read()
        #    print(1)


    #with gzip.open('/home/joe/file.txt.gz', 'rb') as f:
    #    file_content = f.read()
    '''
    arrayid = '204030020070_R05C01'
    xydata = process_xy_file(os.path.join(ARRAY_PATH, arrayid + '.xy'))
    bdata = process_ballele_file(os.path.join(ARRAY_PATH, arrayid + '.b'))
    '''