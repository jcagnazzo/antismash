# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" SANDPUMA. """

import logging
import os, sys
import re
from typing import Any, Dict, List, Optional, Tuple

sys.path.append(os.path.abspath('../../..')) ## just for testing
from antismash.common import fasta, subprocessing, pplacer, module_results
#from antismash.common.secmet import Record

from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier

import numpy as np

from collections import OrderedDict

def get_feature_matrix(spec: str, i2s: List[str]) -> List:
    f = []
    for i in i2s:
        if i == spec:
            f.append(1)
        else:
            f.append(0)
    return f

def run_classifier_test(jackknife_data: str):
    ## Load jackknife data
    jk = {}
    allspec = {}
    with open(jackknife_data, "r") as j:
        next(j) ## skip header
        for line in j:
            line = line.strip()
            l = line.split("\t")
            if l[9] not in jk:
                jk[l[9]] = {'true': l[4],
                            'pid': l[3],
                            'shuf': l[0],
                            'jk': l[1],
                            'query': l[2]}
            called_spec = l[5]
            if l[8] == 'N':
                called_spec = 'no_call'
            if 'method' not in jk[l[9]]:
                jk[l[9]]['method'] = {}
            jk[l[9]]['method'][l[6]] = called_spec
            allspec[l[4]] = -1
            allspec[l[5]] = -1
    ## Map specificities to integers
    i2s = []
    i = 0
    for spec in sorted(allspec, key=allspec.get):
        allspec[spec] = i
        i2s.append(spec)
        i += 1
    ## Prepare features and labels
    allmethods = ['prediCAT', 'forced_prediCAT', 'svm', 'asm', 'pHMM']
    for j in range(1,11):
        test_shuf = 'jk'+str(j)
        for k in range(1,11):
            test_jk = 'k'+str(k)+'_as_query'
            features = []
            labels = []
            test_queries = {}
            for uname in jk:
                if jk[uname]['shuf'] != test_shuf: ## only test within the same shuffle
                    continue
                if jk[uname]['shuf'] == test_shuf:
                    if jk[uname]['jk'] == test_jk: ## test
                        continue
                for m in allmethods:
                    if m in jk[uname]['method']:
                        continue
                    else:
                        jk[uname]['method'][m] = 'no_call'
                labels.append( allspec[jk[uname]['true']] )
                feature_matrix = [ jk[uname]['pid'] ]
                for m in allmethods:
                    feature_matrix.extend(get_feature_matrix(jk[uname]['method'][m], i2s))
                features.append(feature_matrix)
            features = np.array(features).astype(np.float)
            ## Train the ML algorithms
            clf_dt2 = DecisionTreeClassifier(min_samples_leaf=2, max_depth=40)
            clf_dt2 = clf_dt2.fit(features, labels)
            clf_dt5 = DecisionTreeClassifier(min_samples_leaf=5, max_depth=40)
            clf_dt5 = clf_dt5.fit(features, labels)
            clf_dt10 = DecisionTreeClassifier(min_samples_leaf=10, max_depth=40)
            clf_dt10 = clf_dt10.fit(features, labels)
            clf_dt15 = DecisionTreeClassifier(min_samples_leaf=15, max_depth=40)
            clf_dt15 = clf_dt15.fit(features, labels)
            clf_dt20 = DecisionTreeClassifier(min_samples_leaf=20, max_depth=40)
            clf_dt20 = clf_dt20.fit(features, labels)
            
            clf_rf2 = RandomForestClassifier(min_samples_leaf=2, max_depth=40)
            clf_rf2 = clf_rf2.fit(features, labels)
            clf_rf5 = RandomForestClassifier(min_samples_leaf=5, max_depth=40)
            clf_rf5 = clf_rf5.fit(features, labels)
            clf_rf10 = RandomForestClassifier(min_samples_leaf=10, max_depth=40)
            clf_rf10 = clf_rf10.fit(features, labels)
            clf_rf15 = RandomForestClassifier(min_samples_leaf=15, max_depth=40)
            clf_rf15 = clf_rf15.fit(features, labels)
            clf_rf20 = RandomForestClassifier(min_samples_leaf=20, max_depth=40)
            clf_rf20 = clf_rf20.fit(features, labels)
            
            clf_ab = AdaBoostClassifier()
            clf_ab = clf_ab.fit(features, labels)
            
            clf_nb = GaussianNB()
            clf_nb = clf_nb.fit(features, labels)
            
            clf_nn0 = MLPClassifier(alpha=1)
            clf_nn0 = clf_nn0.fit(features, labels)
            clf_nn1 = MLPClassifier(alpha=.1)
            clf_nn1 = clf_nn1.fit(features, labels)
            clf_nn2 = MLPClassifier(alpha=.01)
            clf_nn2 = clf_nn2.fit(features, labels)
            clf_nn3 = MLPClassifier(alpha=.001)
            clf_nn3 = clf_nn3.fit(features, labels)
            
            for uname in jk:
                if jk[uname]['shuf'] == test_shuf:
                    if jk[uname]['jk'] == test_jk: ## test
                        query_features = [jk[uname]['pid']]
                        #print(jk[uname])
                        for m in allmethods:
                            if m not in jk[uname]['method']:
                                jk[uname]['method'][m] = 'no_call'
                            query_features.extend(get_feature_matrix(jk[uname]['method'][m], i2s))
                        query_features = np.array(query_features).astype(np.float).reshape(1, -1)
                        dt2 =  i2s[clf_dt2.predict(query_features)[0]]
                        dt5 =  i2s[clf_dt5.predict(query_features)[0]]
                        dt10 =  i2s[clf_dt10.predict(query_features)[0]]
                        dt15 =  i2s[clf_dt15.predict(query_features)[0]]
                        dt20 =  i2s[clf_dt20.predict(query_features)[0]]
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'DecisionTree_MS02', dt2]))
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'DecisionTree_MS05', dt5]))
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'DecisionTree_MS10', dt10]))
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'DecisionTree_MS15', dt15]))
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'DecisionTree_MS20', dt20]))
                        
                        rf2 =  i2s[clf_rf2.predict(query_features)[0]]
                        rf5 =  i2s[clf_rf5.predict(query_features)[0]]
                        rf10 =  i2s[clf_rf10.predict(query_features)[0]]
                        rf15 =  i2s[clf_rf15.predict(query_features)[0]]
                        rf20 =  i2s[clf_rf20.predict(query_features)[0]]
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'RandomForest_MS02', rf2]))
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'RandomForest_MS05', rf5]))
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'RandomForest_MS10', rf10]))
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'RandomForest_MS15', rf15]))
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'RandomForest_MS20', rf20]))
                        
                        ab =  i2s[clf_ab.predict(query_features)[0]]
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'AdaBoost', ab]))
                        
                        nb =  i2s[clf_nb.predict(query_features)[0]]
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'NaiveBayes', nb]))
                        
                        nn0 =  i2s[clf_nn0.predict(query_features)[0]]
                        nn1 =  i2s[clf_nn1.predict(query_features)[0]]
                        nn2 =  i2s[clf_nn2.predict(query_features)[0]]
                        nn3 =  i2s[clf_nn3.predict(query_features)[0]]
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'NeuralNetwork_A1e0', nn0]))
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'NeuralNetwork_A1e-1', nn1]))
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'NeuralNetwork_A1e-2', nn2]))
                        print("\t".join([test_shuf, test_jk, jk[uname]['query'], jk[uname]['true'], 'NeuralNetwork_A1e-3', nn3]))


run_classifier_test(sys.argv[1])
