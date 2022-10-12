import numpy as np
import pandas as pd
import cooler
import shutil
try:
    from cooltools.lib import numutils
    from cooltools.api.eigdecomp  import eigs_cis
    from cooltools.api.insulation import insulation
    from cooler.core import put,delete,get
except:
    print("If you want call compartment or tad, please install cooltools by 'pip install cooltools'. ")

def _fast_oe(pixels):
    chr_len = pixels.bin2_id.max() - pixels.bin1_id.min()
    pixels["distance"] = pixels["bin2_id"] - pixels["bin1_id"]
    d_sum = pixels.groupby("distance")["count"].sum()
    distance_expected = d_sum / (chr_len - d_sum.index + 1)
    OE = pixels.apply(lambda x : x["count"]/distance_expected[x.distance],axis=1).astype('float32')
    return OE


def __cal_d_weight(mask,px):
    mask_bins = list(mask.index - px.bin1_id.min())
    mask_index = np.zeros((max(mask_bins)+1,1))
    mask_index[mask_bins] = 1
    mask_mat = np.dot(mask_index,mask_index.T)
    d_sum = {}
    for i in range(len(mask_index)):
        mask_sum = np.diag(mask_mat,i).sum()
    d_sum[i] = mask_sum
    d = px.groupby("distance")["count"].sum()
    d = pd.DataFrame(d)
    d["distance"] = d.index
    d["sum"] = d["distance"].apply(lambda x : d_sum.get(x,0))
    d["weight"] = d["sum"] / d["count"]
    return  d["weight"].to_dict()

def _fast_oe_with_gap(pixels,cut_off = 0.1):
    # calculate weight using matrix(chr_len,chr_len), pay attention to your memory
    chr_len = pixels.bin2_id.max() - pixels.bin1_id.min()
    # remove gap which bin sum less than chr_len * cut_off
    num_vector = pixels["bin2_id"].value_counts() + pixels["bin1_id"].value_counts()
    mask = num_vector[num_vector > chr_len * cut_off]
    pixels = pixels[pixels.bin1_id.isin(mask.index) & pixels.bin2_id.isin(mask.index)]
    # calculate OE weight base on distance considering gap
    pixels["distance"] = pixels["bin2_id"] - pixels["bin1_id"]
    d_weight = __cal_d_weight(mask,pixels)
    # calculate OE use count * distance_weight
    OE = pixels.apply(lambda x : x["count"] * d_weight.get(x["distance"],0),axis=1).astype('float32')
    return OE


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

def fast_oe(path,chrom=None,method='intra',gap=False,store=False,newpath=None,verbose = False):
    """
    fast_oe Apply OE balance directly on sparse matrix from cooler files.

    Parameters
    ----------
    path : str 
        path from cooler files
    chrom : str , optional
        A specific chromsome apply OE balance,if None calculate all chromsomes stored in cool file, by default None 
    method : str, optional from ['intra','all']
        Apply OE by different method, if not in below list, apply whole genome wide OE, by default 'intra'
        'intra': Apply OE balancing only for each chromsome internal. 
        'all':  Apply OE balancing for all chromsomes, including inter-chromsome contacts between two chromsome. Time comsuming, pay attention.
    gap : bool, optional
        Whether consider gap which bin sum less than chromsome length * cutoff, by default False
    store : bool, optional
        Whether store OE result in cooler.Cooler().pixels, by default False
    newpath : str, optional
        If newpath and store is not None, create a new .cool file with OE result, by default None

    Returns
    -------
    pd.Dataframe
        A dataframe include pixels and OE result
    """    
    clr = cooler.Cooler(path)
    chroms = clr.chromnames
    OE = np.zeros(clr.info['nnz'])
    if gap == True :
        cal_oe = _fast_oe_with_gap
    else :
        cal_oe = _fast_oe
    if chrom is not None:
        # Apply OE balancing for one chromsome.
        pixels = clr.matrix(as_pixels=True,ignore_index=False,balance=False).fetch(chrom,chrom)
        p = cal_oe(pixels)
        OE[p.index] = p.values
        print("chrom",chrom,"OE balance finished.") if verbose else 0
    elif method == 'intra':
        # Apply OE balancing only for each chromsome internal.    
        for i in chroms:
            try:
                pixels = clr.matrix(as_pixels=True,ignore_index=False,balance=False).fetch(i,i)
                p = cal_oe(pixels)
                OE[p.index] = p.values
                print("chrom",i,"OE balance finished.") if verbose else 0
            except:
                print("chrom",i,"OE balance failed, Skipping.")
    elif method == "all":
        # Apply OE balancing for all chromsomes, including inter-chromsome contacts between two chromsome. Time comsuming, pay attention.
        for i in chroms:
            for j in chroms:
                pixels = clr.matrix(as_pixels=True,ignore_index=False,balance=False).fetch(i,j)
                if len(pixels) > 0: 
                    p = cal_oe(pixels)
                    OE[p.index] = p.values
    else:
        print('Unrecognize OE balancing method, apply genome wide global OE balancing and do not distinguish chromosomes.')
        pixels = clr.pixels()[:]
        OE = cal_oe(pixels)
    OE = pd.DataFrame(OE,columns=["OE"])
    if store:        
        if newpath is not None:
            shutil.copyfile(path,newpath)
            path = newpath
            clr = cooler.Cooler(path)
        with clr.open("r+") as grp:
            delete(grp["pixels"],"OE")
            put(grp["pixels"],OE)
        print("Writing OE normalized matrix to ",path)
    return pd.concat([clr.pixels()[:],OE],axis = 1)
# read OE normalized matrix please use cooler.Cooler(path).matrix(field="OE",balance=False)

def compartment_adj(cool,chrom = None,genome = 'mm10',vec = "E1",fasta = '/store/zlwang/ref/mm10.fa',return_vector=False):
    clr = cooler.Cooler(cool)
    gc_path = f'{genome}_gc_cov_{clr.binsize}.tsv'
    try:
        gc_cov = pd.read_table(gc_path)
    except:
        import bioframe
        print(f"Can't find file {gc_path} in current path, try to generate GC coverage from fasta file by bioframe...")
        bins = clr.bins()[:]
        gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], bioframe.load_fasta(fasta))
        gc_cov.to_csv(gc_path,index=False,sep='\t')
        print(f"Generating {genome} GC covergae file {gc_path} from {fasta} successfully. ")
    if chrom is not None:
        view_df = pd.DataFrame({'chrom': clr.chromnames,'start': 0,'end': clr.chromsizes.values,'name': clr.chromnames})
        view_df = view_df[view_df.chrom==chrom]        
        comp = eigs_cis(clr,phasing_track=gc_cov,view_df=view_df)[1]
        comp = comp[comp.chrom==chrom]
    else: 
        comp = eigs_cis(clr,phasing_track=gc_cov)[1]
    E1 = comp[vec].values
    if return_vector:
        return E1
    E1[E1>1] = 1
    E1[E1<-1] = -1
    return np.outer(E1,E1)

def tad_adj(cool,chrom = None,window_size = 10,return_vector = False):
    # TODO: Multiresolution window size
    clr = cooler.Cooler(cool)
    res = clr.binsize
    if chrom is not None:
        view_df = pd.DataFrame({'chrom': clr.chromnames,'start': 0,'end': clr.chromsizes.values,'name': clr.chromnames})
        view_df = view_df[view_df.chrom==chrom]        
        ins = insulation(clr,window_bp=[res*window_size],view_df=view_df)
        ins = ins[ins.chrom==chrom]
    else :
        ins = insulation(clr,window_bp=[res*window_size])
    score = ins["log2_insulation_score_"+str(res*window_size)]
    boundary = score[ins["is_boundary_"+str(res*window_size)]]
    score = score.values
    score[score>1] = 1
    score[score<-1] = -1
    if return_vector:
        return score,boundary
    boundary_adj = np.zeros((len(score),len(score)))
    for i in range(len(boundary)-1):
        start,end = boundary.index[i],boundary.index[i+1] 
        boundary_adj[start:end,start:end] = 1
    return boundary_adj





