import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.spatial.distance import cosine
from sklearn.metrics import adjusted_rand_score, auc, roc_curve







def auc_curve(y,prob,name):
    fpr,tpr,threshold = roc_curve(y,prob) ###计算真正率和假正率
    roc_auc = auc(fpr,tpr) ###计算auc的值
    plt.figure()
    lw = 2
    plt.figure(figsize=(5,5))
    plt.plot(fpr, tpr, color='darkorange',lw=lw, label='ROC curve (area = %0.3f)' % roc_auc) ###假正率为横坐标，真正率为纵坐标做曲线
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic ' + name )
    plt.legend(loc="lower right")
    plt.show()

def softmax(x):
    return np.exp(x) / np.sum(np.exp(x))
    
def cal_acroc(embs,label):
    labels = pd.get_dummies(label)
    n_label = labels.columns
    for lb in n_label:
        y = labels[lb].values.astype(bool)
        t_emb = (embs[y]).mean(axis = 0)
        f_emb = (embs[~y]).mean(axis = 0)
        prob = [softmax([cosine(i,f_emb),cosine(i,t_emb)])[0] for i in embs]
        auc_curve(y,prob,lb)        


