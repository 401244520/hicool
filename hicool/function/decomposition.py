import numpy as np
import pandas as pd
import cooler
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA,LatentDirichletAllocation,FastICA,NMF

from sklearn.cluster import KMeans,AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score, pairwise_distances
from sklearn.preprocessing import QuantileTransformer
from functools import partial
from multiprocessing import Pool
from tqdm import tqdm
from sklearn.manifold import TSNE

#class Scool_Handler():


def pixels_operation(cool_path,
                     operation="bin_strength",
                     chrom="chr1",
                     comp = None,
                     balance = False,
                     field = "count") -> pd.Series:
    """
    pixels_operation transform 2D matrix to 1D pixels.

    Parameters
    ----------
    cool_path : str
        path to cool file
    operation : list, optional
        A list including various predefined methods, by default ["bin_strength"]
    chrom : str, optional
        Which chrom you want to calculate, by default "chr1"
    comp : list, optional
        A sorted list including boundaries of compartments or tads which can separate bins, by default None
    balance : bool, optional
        Using "weight" column to balance matrix, by default False
    field : str, optional
        Extract which column values, by default "count"

    Returns
    -------
    pd.Series
        A 1D Series including values for each bin 
    """    
    operation_dict = {"bin_degree" : _node_degree,
                      "bin_strength" : _node_strength,
                      "comp_strength" : partial(_node_aggregation_strength,bin=comp),
                      "comp_conservation" : partial(_node_conservation,bin=comp)
    }
    ops = operation_dict[operation]
    clr = cooler.Cooler(cool_path)
    pixels = clr.matrix(balance=balance,as_pixels=True,field=field).fetch(chrom)
    result = ops(pixels)
    return result



def _node_degree(pixels):
    p1 = pixels.value_counts('bin1_id') 
    p2 = pixels.value_counts('bin2_id')
    p = p1.add(p2,fill_value = 0)
    return p

def _node_strength(pixels):
    # node weighted degree
    p1 = pixels.groupby("bin1_id").sum()['count']
    p2 = pixels.groupby("bin2_id").sum()['count']
    p = p1.add(p2,fill_value = 0)
    return p


def _node_aggregation_strength(pixels,bin):
    strength = pd.DataFrame(_node_strength(pixels))
    strength["interval"] = 0
    for b in bin:
        strength.loc[strength.index > b,"interval"] = b
    return strength.groupby("interval").mean()

def _node_conservation(pixels,bin):
    # TODO: A value is trying to be set on a copy of a slice from a DataFrame
    comp = 1
    pixels["bin1_comp"] = 0
    pixels["bin1"] = 0
    pixels["bin2"] = 0
    pixels["bin2_comp"] = 0
    for b in bin:
        pixels.loc[pixels.bin1_id > b,"bin1"] = b
        pixels.loc[pixels.bin2_id > b,"bin2"] = b
        pixels.loc[pixels.bin1_id > b,"bin1_comp"] = comp
        pixels.loc[pixels.bin2_id > b,"bin2_comp"] = comp
        comp *= -1
    # calculate counts of each bin type
    pixels["bin_type"]  = pixels["bin1_comp"] * pixels["bin2_comp"]
    bin1 = pixels.groupby(["bin1","bin_type"]).sum()["count"]
    bin2 = pixels.groupby(["bin2","bin_type"]).sum()["count"]
    bin2.index.names = bin1.index.names
    bins = pd.merge(bin1,bin2,on=bin1.index.names,how="outer").fillna(0).sum(axis=1)
    index_set = set(bins.reset_index("bin_type").index)
    # calculate conservation : AA(or BB) / AB 
    conserve = pd.DataFrame(index = index_set,columns=["conservation"])
    for i in index_set:
        conserve["conservation"][i] = bins[i].get(1,1) / bins[i].get(-1,1)
    return conserve.sort_index()


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
