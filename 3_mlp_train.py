import pandas as pd
import numpy as np
from sklearn.preprocessing import QuantileTransformer
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.model_selection import KFold
from sklearn.linear_model import RidgeClassifier
from sklearn.metrics import roc_auc_score
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score
import xgboost as xgb
from sklearn.model_selection import StratifiedKFold
import cross_validation
import pickle5 as pickle
# {'activation': 'tanh', 'alpha': 0.0001, 'hidden_layer_sizes': (50, 50, 50), 'learning_rate': 'adaptive', 'solver': 'adam'}


clf_instances = [RandomForestClassifier(n_jobs=20),
                 MLPClassifier(),
                 xgb.XGBClassifier(objective='multi:softprob', random_state=42,n_jobs=20)
                 ]

clf_names = ['RF',
             'MLP',
             'XGB']


features_to_exclude=['father_m_raw',
                    'mother_m_raw',
                    'scc_m_raw',
                    'father_gtype',
                    'mother_gtype',
                    'scc_gtype',
                    'scc_output']

res_ar=[]

family_data = cross_validation.load()

clf=MLPClassifier()

for train, _ in cross_validation.fold_iterator(family_data,cross_val=False):
        train=train[[c for c in train.columns if c not in features_to_exclude]]
        ## normalize data
        qt = QuantileTransformer(n_quantiles=10, random_state=0)
        for c in train.columns:
            if 'a_raw' in c or 'x_raw'  in c or 'y_raw' in c:
                 train[c]=qt.fit_transform(train[[c]])

        X_train=train[[c for c in train.columns if c!='ref_proband_gtype']]
        y_train=train['ref_proband_gtype']

        print(X_train.columns)
        clf.fit(X_train, y_train)

        with open('MLP.pickle', 'wb') as f:
            pickle.dump(clf, f)



