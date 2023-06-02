import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.spatial.distance import cosine
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, roc_curve, auc

def evaluate_clustering(y_true, y_pred):
    ari = adjusted_rand_score(y_true, y_pred)
    nmi = normalized_mutual_info_score(y_true, y_pred)
    return {'ARI': ari, 'NMI': nmi}

def cal_clustering(embs,label_list,cluster_method,evaluate_method):
    label_dict = {label:i for i,label in enumerate(sorted(set(label_list)))}
    labels = [label_dict[label] for label in label_list]
    n_label = len(label_dict)
    pred = cluster_method(embs,n_label)
    evals = evaluate_method(labels,pred)
    return evals

def auc_curve(y,prob,name,plot=True):
    fpr,tpr,threshold = roc_curve(y,prob) ###计算真正率和假正率
    roc_auc = auc(fpr,tpr) ###计算auc的值
    if plot:
        plt.figure()
        lw = 2
        plt.figure(figsize=(5,5))
        plt.plot(fpr, tpr, color='darkorange',lw=lw, label='ROC curve (area = %0.3f)' % roc_auc) ###假正率为横坐标，真正率为纵坐标做曲线
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(name)
        plt.legend(loc="lower right")
        plt.show()
    return roc_auc,fpr,tpr

def softmax(x):
    return np.exp(x) / np.sum(np.exp(x))
    
def cal_acroc(embs,label,plot=True):
    """
    cal_acroc Calculate Average Circular ROC while ROC can not directly apply on Multi 

    Parameters
    ----------
    embs : _type_
        _description_
    label : _type_
        _description_
    plot : bool, optional
        _description_, by default True

    Returns
    -------
    _type_
        _description_
    """
    labels = pd.get_dummies(label)
    acroc,acfpr,actpr = {},{},{}
    n_label = labels.columns
    for lb in n_label:
        y = labels[lb].values.astype(bool)
        t_emb = (embs[y]).mean(axis = 0)
        f_emb = (embs[~y]).mean(axis = 0)
        prob = [softmax([cosine(i,f_emb),cosine(i,t_emb)])[0] for i in embs]
        roc_auc,fpr,tpr = auc_curve(y,prob,lb,plot=False)        
        acroc[lb] = roc_auc
        acfpr[lb] = fpr
        actpr[lb] = tpr
    if plot :
        plt.figure()
        lw = 2
        plt.figure(figsize=(10,10),dpi=120)
        for lb in n_label:
            plt.plot(acfpr[lb], actpr[lb],lw=lw, label= lb+' (acroc = %0.3f)' % acroc[lb]) 
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.legend(loc="lower right")
    return acroc

