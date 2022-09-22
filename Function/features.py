import numpy as np
import pandas as pd
import cooler

#from cooltools.lib import numutils
from cooltools.api.eigdecomp  import eigs_cis
from cooltools.api.insulation import insulation
from cooler.core import put,delete,get


def _fast_oe(pixels):
    pixels["distance"] = pixels["bin2_id"] - pixels["bin1_id"]
    gr = pixels.groupby("distance")["count"]
    distance_expected = gr.mean().to_dict()
    OE = pixels.apply(lambda x : x["count"]/distance_expected[x.distance],axis=1).astype('float32')
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

def fast_oe(path,chrom=None,method='intra',store=False):
    clr = cooler.Cooler(path)
    chroms = clr.chromnames
    OE = np.zeros(clr.info['nnz'])
    if chrom is not None:
        # Apply OE balancing for one chromsome.
        pixels = clr.matrix(as_pixels=True,ignore_index=False,balance=False).fetch(chrom,chrom)
        p = _fast_oe(pixels)
        OE[p.index] = p.values
        print("chrom",chrom,"OE balance finished.")
    elif method == 'intra':
        # Apply OE balancing only for each chromsome internal.    
        for i in chroms:
            pixels = clr.matrix(as_pixels=True,ignore_index=False,balance=False).fetch(i,i)
            p = _fast_oe(pixels)
            OE[p.index] = p.values
            print("chrom",i,"OE balance finished.")
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
        print("Writing OE normalized matrix to ",path)
    return pd.concat([clr.pixels()[:],OE],axis = 1)
# read OE normalized matrix please use cooler.Cooler(path).matrix(field="OE",balance=False)

def compartment_adj(cool,chrom = None,genome = 'mm10',fasta = '/store/zlwang/ref/mm10.fa',return_vector=False):
    clr = cooler.Cooler(cool)
    gc_path = f'{genome}_gc_cov_{clr.binsize}.tsv'
    try:
        gc_cov = pd.read_table(gc_path)
    except:
        import bioframe
        print(f"Can't find file {gc_path} in current path, please provide correct genome or fasta files.")
        bins = clr.bins()[:]
        gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], bioframe.load_fasta(fasta))
        gc_cov.to_csv(gc_path,index=False,sep='\t')
        print(f"Generating {genome} GC covergae file {gc_path} successfully. ")
    if chrom is not None:
        view_df = pd.DataFrame({'chrom': clr.chromnames,'start': 0,'end': clr.chromsizes.values,'name': clr.chromnames})
        view_df = view_df[view_df.chrom==chrom]        
        comp = eigs_cis(clr,phasing_track=gc_cov,view_df=view_df)[1]
        comp = comp[comp.chrom==chrom]
    else: 
        comp = eigs_cis(clr,phasing_track=gc_cov)[1]
    E1 = comp["E1"].values
    E1[E1>1] = 1
    E1[E1<-1] = -1
    if return_vector:
        return E1
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





