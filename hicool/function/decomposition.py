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

def pixels_operation(cool_path,
                     operation="bin_strength",
                     chrom="chr1",
                     regions = None,
                     balance = False,
                     field = "count",
                     **kwargs) -> pd.Series:
    """
    pixels_operation transform 2D matrix to 1D pixels.

    Parameters
    ----------
    cool_path : str
        path to cool file
    operation : str, optional
        Choose one method from various predefined methods or customize, by default "bin_strength"
    chrom : str, optional
        Which chrom you want to calculate, by default "chr1"
    regions : list, optional
        A sorted list including genome regions such as  ["chr1:1001000-1002000","chr2:1234000-1244000"], by default None
    balance : bool, optional
        Using "weight" column to balance matrix, by default False
    field : str, optional
        Extract which column values, by default "count"

    Returns
    -------
    pd.Series
        A 1D Series including values for each bin 
    """    
    operation_dict = {"bin_degree" :  partial(node_degree,chrom=chrom,balance=balance,field=field),
                      "bin_strength" :  partial(node_strength,chrom=chrom,balance=balance,field=field),
                      "region_strength" : partial(strength_region,regions=regions,balance=balance,field=field,**kwargs),
                      "region_conservation" : partial(conservation_region,regions=regions,balance=balance,field=field,**kwargs),
                      "comp_strength" : partial(_node_aggregation_strength,**kwargs),
                      "comp_conservation" : partial(_node_conservation,**kwargs),
                      "strata_strength" : partial(strata_strength,chrom=chrom,balance=balance,field=field,**kwargs),
    }
    ops = operation_dict[operation]
    # if operation in ["bin_degree","bin_strength"]:
    #     clr = cooler.Cooler(cool_path)
    #     pixels = clr.matrix(balance=balance,as_pixels=True,field=field).fetch(chrom)
    #     result = ops(pixels)
    # else:
    result = ops(cool_path) 
    return result

def loop_APA(cool_path,loop_anchors,balance,field,width=20,pos="lower-right"):
    clr = cooler.Cooler(cool_path)
    binsize = clr.binsize
    APA = []
    for loop in loop_anchors:
        chromosome, interval = loop.split(":")
        start, end = map(int, interval.split("-"))
        len = end - start
        # min_len = binsize * width * 1.5
        min_len = 200000
        if len < min_len:
            APA.append(0)
            break
        else:
            r1 = chromosome + ":" + str(start - width*binsize) + "-" + str(start + (1+width)*binsize)
            r2 = chromosome + ":" + str(end - width*binsize) + "-" + str(end + (1+width)*binsize)
            mat = clr.matrix(balance=balance,field=field).fetch(r1,r2)
        center = mat[width,width]
        n = width // 2
        if pos == None :
            APA.append(center)
            break
        elif pos == "lower-right" :
            backgroud = mat[-n:,-n:].mean()
        elif pos == "lower-left" :
            backgroud = mat[:n,-n:].mean()
        elif pos == "upper-left" :
            backgroud = mat[:n,:n].mean()
        elif pos == "upper-right" :
            backgroud = mat[-n:,:n].mean()
        elif pos == "around" :
            backgroud = (mat[-n:,-n:]+mat[:n,-n:]+mat[:n,:n]+mat[-n:,:n]).mean()
        else :
            raise Exception(f"Not supported position {pos} !")
        if backgroud == 0 :
            APA.append(center)
            break
        else :
            APA.append(center/backgroud)
            break
    return pd.DataFrame(APA)

def strata_strength(cool_path,chrom,balance,field,n_strata=20):
    all_strata = []
    mat = cooler.Cooler(cool_path).matrix(balance=balance,field=field).fetch(chrom)[:]
    for i in range(n_strata):
        strata = pd.DataFrame(np.diag(mat,i))
        strata.index = chrom + "-" + str(i) + "-" + strata.index.astype(str)
        all_strata.append(strata)
    return pd.concat(all_strata)

def conservation_region(cool_path,regions,balance,field,intra_chrom=True):
    result = []
    clr = cooler.Cooler(cool_path)
    for region in regions:
        intra_region = clr.matrix(balance=balance,field=field,as_pixels=True).fetch(region)[field].sum()
        if intra_chrom :
            chrom = region.split(":")[0]
            background_region = clr.matrix(balance=balance,field=field,as_pixels=True).fetch(region,chrom)[field].sum()
        else:
            # backgroud set as genomewide without balance
            background_region = clr.pixels().fetch(region)[field].sum()
        if background_region != 0 :
            conservation = intra_region/background_region
        else:
            conservation = 0
        result.append(conservation)
    return pd.DataFrame(result)

def strength_region(cool_path,regions,balance,field,intra_chrom=True,percentage=False):
    result = []
    clr = cooler.Cooler(cool_path)
    for region in regions:
        if intra_chrom:
            result.append(clr.matrix(balance=balance,field=field).fetch(region).mean())
        else:
            chrom = region.split(":")[0]
            intra = clr.matrix(balance=balance,field=field,as_pixels=True).fetch(region,chrom)[field].sum()
            background = clr.pixels().fetch(region)[field].sum()
            if background == 0:
                result.append(0)
            else:
                result.append(intra/background)
    return pd.DataFrame(result)

def node_degree(cool_path,chrom,balance,field):
    clr = cooler.Cooler(cool_path)
    pixels = clr.matrix(balance=balance,as_pixels=True,field=field).fetch(chrom)
    p1 = pixels.value_counts('bin1_id') 
    p2 = pixels.value_counts('bin2_id')
    p = p1.add(p2,fill_value = 0)
    return p

def node_strength(cool_path,chrom,balance,field):
    # node weighted degree
    clr = cooler.Cooler(cool_path)
    pixels = clr.matrix(balance=balance,as_pixels=True,field=field).fetch(chrom)
    p1 = pixels.groupby("bin1_id").sum()['count']
    p2 = pixels.groupby("bin2_id").sum()['count']
    p = p1.add(p2,fill_value = 0)
    return p



############################################################################################

def _node_strength_region(cool_path,regions,balance,field):
    result = []
    clr = cooler.Cooler(cool_path)
    for region in regions:
        result.append(clr.matrix(balance=balance,field=field).fetch(region).mean())
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

def _node_conservation(cool,bin,chrom="chr1"):
    # TODO: A value is trying to be set on a copy of a slice from a DataFrame
    pixels = cooler.Cooler(cool).pixels().fetch(chrom)
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
        conserve["conservation"][i] = bins[i].get(1,1) / (bins[i].get(-1,1) +  bins[i].get(1,1))
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
