import numpy as np
import pandas as pd
import cooler

#from cooltools.lib import numutils
from cooltools.api.eigdecomp  import cis_eig
from cooltools.api.insulation import insulation
from cooler.core import put,delete,get

def _fast_oe(pixels):
    pixels["distance"] = pixels["bin2_id"] - pixels["bin1_id"]
    gr = pixels.groupby("distance")["count"]
    distance_expected = gr.mean().to_dict()
    OE = pixels.apply(lambda x : x["count"]/distance_expected[x.distance],axis=1).astype('float32')
    return OE
def fast_oe(path,method='intra',store=False):
    clr = cooler.Cooler(path)
    chroms = clr.chromnames
    OE = np.zeros(clr.info['nnz'])
    if method == 'intra':
        # Apply OE balancing only for each chromsome internal.    
        for i in chroms:
            pixels = clr.matrix(as_pixels=True,ignore_index=False,balance=False).fetch(i,i)
            p = _fast_oe(pixels)
            OE[p.index] = p.values
    elif method == "all":
        # Apply OE balancing for all chromsomes, including inter-chromsome contacts between two chromsome. Time comsuming, pay attention.
        for i in chroms:
            for j in chroms:
                pixels = clr.matrix(as_pixels=True,ignore_index=False,balance=False).fetch(i,j)
                if len(pixels) > 0: 
                    p = _fast_oe(pixels)
                    OE[p.index] = p.values
    else:
        print('Unrecognize OE balancing method, apply genomewide global OE balancing.')
        pixels = clr.pixels()[:]
        OE = _fast_oe(pixels)
    OE = pd.DataFrame(OE,columns=["OE"])
    if store:        
        with clr.open("r+") as grp:
            delete(grp["pixels"],"OE")
            put(grp["pixels"],OE)
    return pd.concat([clr.pixels()[:],OE],axis = 1)















def cal_all_strata(cell_list,chrom = "chr1",n_strata = 20):
    all_strata = [[] for i in range(n_strata)]
    for cell in tqdm(cell_list):
        mat = cooler.Cooler(cell).matrix(balance=False).fetch(chrom)[:]
        #  mat = matrix_operation(mat,['convolution'])
        for i in range(n_strata):
            all_strata[i].append(np.diag(mat,i))
    all_strata = [np.array(strata) for strata in all_strata]
    return all_strata


def cal_hicrep(chrom = "chr1",n_strata = 20,method = 'hicrep'):
    all_strata = cal_all_strata(cell_list,chrom,n_strata)
    sim = pairwise_distances(all_strata,method)
    return sim

