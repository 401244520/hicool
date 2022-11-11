import numpy as np
import pandas as pd
import cooler
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA,LatentDirichletAllocation,FastICA,NMF

from sklearn.cluster import KMeans,AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score, pairwise_distances
from sklearn.preprocessing import QuantileTransformer

from sklearn.manifold import TSNE

def node_degree(pixels):
    p1 = pixels.value_counts('bin1_id') 
    p2 = pixels.value_counts('bin2_id')
    p = p1.add(p2,fill_value = 0)
    return p
def node_strength(pixels):
    # node weighted degree
    p1 = pixels.groupby("bin1_id").sum()['count']
    p2 = pixels.groupby("bin2_id").sum()['count']
    p = p1.add(p2,fill_value = 0)
    return p

def all_strata(cell_list,chrom = "chr1",n_strata = 20,field = "count",concat = False):
    all_strata = [[] for i in range(n_strata)]
    for cell in cell_list:
        # When resolution is 1KB, the matrix is too large to store in numpy.array. Suggest use MinHash.
        mat = cooler.Cooler(cell).matrix(balance=False,field=field).fetch(chrom)[:]
        #mat = matrix_operation(mat,['oe','randomwalk','convolution'])
        for i in range(n_strata):
            all_strata[i].append(np.diag(mat,i))
    all_strata = [np.array(strata) for strata in all_strata]
    if concat :
        concat_strata = all_strata[0]
        for i in all_strata[1:]:
            concat_strata = np.concatenate((concat_strata,i),axis = 1)
        return concat_strata
    return all_strata
