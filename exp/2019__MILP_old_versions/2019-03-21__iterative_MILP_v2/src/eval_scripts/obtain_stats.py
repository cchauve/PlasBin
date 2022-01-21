import os
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.switch_backend('agg')
output_dir = '/home/amane/projects/rrg-chauvec/wg-anoph/Plasmids-Assembly/exp/2019-02-19__plasmid_MILP/output'
ratio_tests = os.path.join(output_dir,'ratio_test')

files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(ratio_tests) for f in filenames if "eval.csv" in f]

mean_scores = []
best_scores = []
precs = defaultdict(list)
recs = defaultdict(list)
f1s = defaultdict(list)
mean, best = {}, {}  
for file in files:
    with open(file, 'r') as f:
        for line in f:
            if "precision" in line:
                precision = line.split(" ")[-1]
                precs[file.split('/')[11]].append(float(precision))
            if "recall" in line:
                recall = line.split(" ")[-1]
                recs[file.split('/')[11]].append(float(recall))
            if "f1" in line:                                                          
                f1 = line.split(" ")[-1]
                f1s[file.split('/')[11]].append(float(f1)) 

for ratio in precs:
    mean[ratio] = {}
    mean[ratio]['precision'] = sum(precs[ratio])/len(precs[ratio])
    mean[ratio]['recall'] = sum(recs[ratio])/len(recs[ratio])
    mean[ratio]['f1_score'] = sum(f1s[ratio])/len(f1s[ratio])
    mean_scores.append([ratio, mean[ratio]['precision'], mean[ratio]['recall'], mean[ratio]['f1_score']])

    best[ratio] = {}
    best[ratio]['precision'] = max(precs[ratio])
    best[ratio]['recall'] = max(recs[ratio])
    best[ratio]['f1_score'] = max(f1s[ratio])
    best_scores.append([ratio, best[ratio]['precision'], best[ratio]['recall'], best[ratio]['f1_score']])

mean_scores = pd.DataFrame(mean_scores)
mean_scores.rename(columns = {0: 'Ratio', 1: 'Precision', 2: 'Recall', 3: 'F1 score'}, inplace = True)

best_scores = pd.DataFrame(best_scores)
best_scores.rename(columns = {0: 'Ratio', 1: 'Precision', 2: 'Recall', 3: 'F1 score'}, inplace = True)

with open('exp1_scores.csv', 'w') as f:
    mean_scores.to_csv(f, sep = '\t', encoding='utf-8', index=False)
with open('exp1_scores.csv', 'a') as f:
    best_scores.to_csv(f, sep = '\t', encoding='utf-8', index=False)

#print(mean)
#print(best)

files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(output_dir) for f in filenames if "eval.csv" in f]
precs = defaultdict(list)
recs = defaultdict(list)
f1s = defaultdict(list)
for file in files:
    if 'sample' in file.split('/')[10]:
        sample_id = file.split('/')[10].split('_')[1]
        with open(file, 'r') as f:
            for line in f:
                if "precision" in line:
                    precision = line.split(" ")[-1]
                    precs[sample_id].append(float(precision))
                if "recall" in line:
                    recall = line.split(" ")[-1]
                    recs[sample_id].append(float(recall))
                if "f1" in line:
                    f1 = line.split(" ")[-1]
                    f1s[sample_id].append(float(f1))

mean, best = {}, {}
mean_scores, best_scores = [], [] 
for sample_id in precs:
        mean[sample_id] = {}
        mean[sample_id]['precision'] = sum(precs[sample_id])/len(precs[sample_id])
        mean[sample_id]['recall'] = sum(recs[sample_id])/len(recs[sample_id])
        mean[sample_id]['f1_score'] = sum(f1s[sample_id])/len(f1s[sample_id])
        mean_scores.append([sample_id, mean[sample_id]['precision'], mean[sample_id]['recall'], mean[sample_id]['f1_score']])

        best[sample_id] = {}
        best[sample_id]['precision'] = max(precs[sample_id])
        best[sample_id]['recall'] = max(recs[sample_id])
        best[sample_id]['f1_score'] = max(f1s[sample_id])
        best_scores.append([sample_id, best[sample_id]['precision'], best[sample_id]['recall'], best[sample_id]['f1_score']])

mean_scores = pd.DataFrame(mean_scores)
mean_scores.rename(columns = {0: 'Sample', 1: 'Precision', 2: 'Recall', 3: 'F1 score'}, inplace = True)

best_scores = pd.DataFrame(best_scores)
best_scores.rename(columns = {0: 'Sample', 1: 'Precision', 2: 'Recall', 3: 'F1 score'}, inplace = True)

#with open('exp1_scores.csv', 'w') as f:
#    mean_scores.to_csv(f, sep = '\t', encoding='utf-8', index=False)
#with open('exp1_scores.csv', 'a') as f:
#    best_scores.to_csv(f, sep = '\t', encoding='utf-8', index=False)

print(mean_scores['Recall'].values.tolist())

N = 10
ind = np.arange(N)
width = 0.27

fig = plt.figure()
ax = fig.add_subplot(111)
#fig, ax = plt.subplot(nrows = 1, figsize = (45, 30))

pvals = mean_scores['Precision'].values.tolist()
rects1 = ax.bar(ind, pvals, width, color='r')
rvals = mean_scores['Recall'].values.tolist()
rects2 = ax.bar(ind+width, rvals, width, color='g')
fvals = mean_scores['F1 score'].values.tolist()
rects3 = ax.bar(ind+width*2, fvals, width, color='b')
ids = mean_scores['Sample'].values.tolist()
ax.set_ylabel('Scores')
ax.set_xticks(ind+width)
ax.set_xticklabels( (ids[0], ids[1], ids[2], ids[3], ids[4], ids[5], ids[6], ids[7], ids[8], ids[9]) )
ax.legend( (rects1[0], rects2[0], rects3[0]), ('Prec', 'Rec', 'F1') )
'''
def autolabel(rects):
    for rect in rects:
        h = rect.get_height()
        #ax.text(rect.get_x()+rect.get_width()/2., 1.05*h, '%d'%float("{0:.2f}".format(h)),ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)
autolabel(rects3)
'''
plt.savefig('mean_scores_MILP_exp1.pdf', format = 'pdf', dpi = 1200, bbox_inches = 'tight')

fig = plt.figure()
ax = fig.add_subplot(111)
#fig, ax = plt.subplot(nrows = 1, figsize = (45, 30))

pvals = best_scores['Precision'].values.tolist()
rects1 = ax.bar(ind, pvals, width, color='r')
rvals = best_scores['Recall'].values.tolist()
rects2 = ax.bar(ind+width, rvals, width, color='g')
fvals = best_scores['F1 score'].values.tolist()
rects3 = ax.bar(ind+width*2, fvals, width, color='b')
ids = best_scores['Sample'].values.tolist()
ax.set_ylabel('Scores')
ax.set_xticks(ind+width)
ax.set_xticklabels( (ids[0], ids[1], ids[2], ids[3], ids[4], ids[5], ids[6], ids[7], ids[8], ids[9]) )
ax.legend( (rects1[0], rects2[0], rects3[0]), ('Prec', 'Rec', 'F1') )

plt.savefig('best_scores_MILP_exp1.pdf', format = 'pdf', dpi = 1200, bbox_inches = 'tight')
