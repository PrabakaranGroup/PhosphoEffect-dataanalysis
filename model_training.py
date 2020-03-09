# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 14:52:34 2019

@author: steph
"""

import sklearn
import pandas as pd
import random
import numpy as np
import re
import urllib
import urllib.request

benign_vars = pd.read_csv("humvar/humvar-2011_12.neutral.humvar.output",sep="\t").iloc[:,:4].values
scores = pd.read_csv("humvar/humvar-2011_12.neutral.humvar.output",sep="\t").iloc[:,15:16].values
benign_vars = np.hstack((benign_vars,scores))

rows = random.sample(list(range(len(benign_vars))),10000)
variants = benign_vars[rows,:]
variants_to_keep = []


for i in range(len(variants)):
    url = "http://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{0}".format(variants[i,0])
    try:
        pdbfile = urllib.request.urlopen(urllib.request.Request(url)).read()
        variants_to_keep.append(i)
    except urllib.error.HTTPError:
        continue
variants = variants[variants_to_keep,:]

a=[]
for v in range(len(variants)):
    i = variants[v]
    for j in [0,2,3]:
        i[j] = re.sub("\s","",i[j])
    try:
        i[4] = float(i[4])
    except ValueError:
        a.append(v)             
        
variants = np.delete(variants,a,0)

class_1=[]
class_2=[]
class_3=[]
for i in range(len(variants)):
    v = variants[i]
    if v[0] in uniprot_IDs_named:
        index = list(uniprot_IDs_named[:,0]).index(v[0])
        gene = uniprot_IDs_named[index,2]
        position = v[1]
        if int(position) in gene_phospho.get(gene,[]):
            class_1.append(i)
        elif any(i in gene_phospho.get(gene,[]) for i in range(int(position)-5,int(position)+6)):
            class_2.append(i)
        elif gene in gene_phospho:
            class_3.append(i)
        

        
rows_1 = class_1
rows_2 = random.sample(class_2,424)
rows_3 = random.sample(class_3,250)
rows = rows_1+rows_2+rows_3
variants_final = variants[rows,:]

train_rows = random.sample(list(range(750)),600)
test_rows = [i for i in range(750) if not i in train_rows]

variants_train = variants_final[train_rows,:]
variants_test = variants_final[test_rows,:]

np.savetxt("benign_variants_training.txt",variants_train,delimiter="\t",fmt="%s")

np.savetxt("benign_variants_test.txt",variants_test,delimiter="\t",fmt="%s")



        

rows = random.sample(list(range(len(benign_vars))),5000)
training_rows = rows[:4000]
test_rows = rows[4000:]

benign_train = benign_vars[training_rows,:]
benign_test = benign_vars[test_rows,:]

pathogenic_vars = pd.read_csv("pathogenic_variants_formatted.txt",sep="\t").iloc[:,:].values
n = len(pathogenic_vars)
training_rows = random.sample(list(range(n)),int(n*0.8))
test_rows = [i for i in range(n) if not i in training_rows]


pathogenic_train = pathogenic_vars[training_rows,:]
pathogenic_test = pathogenic_vars[test_rows,:]

all_train = np.vstack((benign_train,pathogenic_train))

all_train = all_train[:,:4]
for i in all_train:
    for j in [0,2,3]:
        i[j] = re.sub("\s","",i[j])
np.savetxt("training_set.txt", all_train, fmt="%s",delimiter="\t")
