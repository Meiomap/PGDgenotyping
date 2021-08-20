import sys
import numpy as np
import os
from pathlib import Path
import pandas as pd
import seaborn as sns; sns.set_theme(color_codes=True)
import SureTypeSC as sc
import logging
import SureTypeSC.MachineLearning as MachineLearning
import SureTypeSC.Config as Config
import pickle5 as pickle
#import pickle
from SureTypeSC.DataLoader import Patterns
from itertools import product
import random
from collections import Counter
import re


####
def load(f='../Family_774.pkl'):
  with open(f, "rb") as fh:
    df = pickle.load(fh)

    df=df[df.apply(lambda x: 'cnv' not in x['Name'],axis=1)]
    dfs=sc.Data.create_from_frame(df)

    gdna=dfs.slice('gdna')
    scc=dfs.slice('sc')

    gdna.apply_NC_threshold_3(0.15,inplace=True)
    mother=gdna.slice("gm07224")
    father=gdna.slice("gm07225")
    ref_proband=gdna.slice("gm07228")

    for d in [mother,father]:
        consistent_calls=(d.df.loc[:,(slice(None),'gtype')].eq(d.consensus_genotype(),axis=0).all(axis=1)) & (d.consensus_genotype()!='NC')
        d.df=d.df.loc[consistent_calls]


    for el in [mother,father,scc]:
        el.calculate_transformations_2()

    scc.compare_against_reference_lightweight(ref_proband)

  return (mother,father,scc,ref_proband)

def generate_fold(data,sample_combinations):
    mother,father,scc,ref_proband=data
    features_parents = ['b_allele_freq', 'x_raw', 'y_raw', 'm_raw', 'a_raw', 'gtype']
    features_embryos = ['b_allele_freq', 'x_raw', 'y_raw', 'm_raw', 'a_raw', 'gtype', 'output']
    features_proband = ['gtype']
    list_to_concat = []
    for (mother_s, father_s, embryo_s, proband_s) in sample_combinations:
        _df = pd.concat([mother.df.loc[:, (mother_s, features_parents)],
                         father.df.loc[:, (father_s, features_parents)],
                         scc.df.loc[:, (embryo_s, features_embryos)],
                         ref_proband.df.loc[:, (proband_s, features_proband)]], axis=1)

        # print(mother.df.loc[:,(mother_s,features_parents)].columns)

        # _df.columns=['_'.join(list(c)) for c in ren_columns]
        ren_columns = ['mother'] * len(features_parents) + ['father'] * len(features_parents) + ['scc'] * len(
            features_embryos) + ['ref_proband'] * len(features_proband)
        assert (len(_df.columns) == len(ren_columns))
        _df.columns = [sample + '_' + orig_c[1] for sample, orig_c in zip(ren_columns, _df.columns)]
        list_to_concat.append(_df)

    dataset = pd.concat(list_to_concat, axis=0)

    dataset = dataset[
        (dataset.mother_gtype != 'NC') & (dataset.father_gtype != 'NC') & (dataset.scc_gtype != 'NC') & (
                    dataset.ref_proband_gtype != 'NC')].dropna()
    dataset.scc_output = dataset.scc_output.apply(int)
    g = dataset.groupby(['mother_gtype', 'father_gtype', 'scc_output'])
    a = g.apply(lambda x: x.sample(g.size().min()).reset_index(drop=True))
    a = a.reset_index(drop=True)
    return a

def fold_iterator(data,cross_val=True):
    mother, father, scc, ref_proband = data
    sc_individuals = set(scc.df.columns.get_level_values(0))
    mother_individuals = set(mother.df.columns.get_level_values(0))
    father_individuals = set(father.df.columns.get_level_values(0))
    ref_proband_individuals = set(ref_proband.df.columns.get_level_values(0))
    random.seed(2)
    a = [(m, f, s, r) for m, f, s, r in zip(random.sample(mother_individuals, 9),
                                            random.sample(father_individuals, 9),
                                            random.sample(sc_individuals, 9),
                                            random.sample(ref_proband_individuals, 9))]


    #split into training and testing
    random.seed(12)
    indices=[i for i in range(9)]
    random.shuffle(indices)
    if cross_val == True:
      train=indices[0:7]
      test=indices[7:]
    else:
        train=indices
        test=indices
    yield (generate_fold(data,np.array(a)[train]),generate_fold(data,np.array(a)[test]))


if __name__=='__main__':
  for run in range(0,9):
    family_data=load()
    for train,test in fold_iterator(family_data):
        print(run)
        print(train)
        print(test)
