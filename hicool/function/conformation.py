import shutil
import cooler
import numpy as np
import pandas as pd
from functools import partial

try:
    from cooler.core import delete, get, put
    from cooltools.api.eigdecomp import eigs_cis
    from cooltools.api.insulation import insulation
except:
    print("If you want call compartment or tad, please install cooltools by 'pip install cooltools' (not avaliable on Windows). ")

def calculate_feature(cool_path,
                     operation="comp",
                     chrom="chr1",
                     **kwargs) -> pd.Series:
    
    if operation in ['comp','compartment']:
        result = compartment_decomposition(cool_path,chrom=chrom,**kwargs)
    elif operation in ['tad']:
        result = tad_insulation(cool_path,chrom=chrom,**kwargs)
    else :
        print(f"Not a valid feature, now only support compartment and tad.")    
    return result

def load_gc_cov(gc_path, chrom, clr):
    try:
        gc_cov = pd.read_table(gc_path)
        if (chrom not in set(gc_cov["chrom"])) and not set(clr.chromnames).issubset(set(gc_cov["chrom"])):
            raise Exception(f"Chromsome {chrom} not in {gc_path}")
        return gc_cov
    except:
        return None

def generate_gc_cov(bins, fasta):
    import bioframe
    print(f"Generating GC coverage from fasta file by bioframe...")
    return bioframe.frac_gc(bins[['chrom', 'start', 'end']], bioframe.load_fasta(fasta))

def create_view_df(path, chrom):
    clr = cooler.Cooler(path)
    view_df = pd.DataFrame({'chrom': clr.chromnames,
                            'start': 0,
                            'end': clr.chromsizes.values,
                            'name': clr.chromnames})
    if chrom is not None:
        view_df = view_df[view_df.chrom == chrom].reset_index(drop=True)
    return view_df

def compartment_decomposition(path, chrom=None, vec="E1",
                              gc_path=None, genome='mm10',fasta='~/ref/mm10.fa',
                              return_vector=True, store=False):
    clr = cooler.Cooler(path)
    bins = clr.bins()[:]
    if chrom is not None:
        bins = bins[bins["chrom"] == chrom]
    if gc_path is None:
        gc_path = f'{genome}_gc_cov_{clr.binsize}.tsv'
    gc_cov = load_gc_cov(gc_path, chrom, clr)
    if gc_cov is None:
        gc_cov = generate_gc_cov(bins, fasta)
        gc_cov.to_csv(gc_path, index=False, sep='\t')
        print(f"Generated {genome} GC coverage file {gc_path} from {fasta} successfully.")
    view_df = create_view_df(path, chrom)
    comp = eigs_cis(clr, phasing_track=gc_cov)[1] if chrom is None else eigs_cis(clr, phasing_track=gc_cov, view_df=view_df)[1]
    if store and chrom is None:
        with clr.open("r+") as grp:
            put(grp["bins"], comp['E1'])
            put(grp["bins"], comp['E2'])
            put(grp["bins"], comp['E3'])
        print("Writing compartments score to ", path)
    Eig = comp[vec].values
    return Eig if return_vector else comp


def tad_insulation(path, chrom=None, window_size=10,
                   return_vector=True, store=False):
    clr = cooler.Cooler(path)
    res = clr.binsize
    view_df = create_view_df(clr, chrom)
    ins = insulation(clr, window_bp=[res*window_size], view_df=view_df)
    if chrom is not None:
        ins = ins[ins.chrom == chrom]
    score = ins["log2_insulation_score_"+str(res*window_size)]
    boundary_strength = ins["boundary_strength_"+str(res*window_size)]
    boundary = ins["is_boundary_"+str(res*window_size)]
    if store and (chrom is None):
        with clr.open("r+") as grp:
            put(grp["bins"], score)
            put(grp["bins"], boundary)
            put(grp["bins"], boundary_strength)
        print("Writing insulation score and boundary to ", path)
    boundary = score[ins["is_boundary_"+str(res*window_size)]]
    score = score.values
    return score if return_vector else ins



#  **************************************************old_version**********************************************  #

def _sparse_oe(pixels):
    chr_len = pixels.bin2_id.max() - pixels.bin1_id.min()
    pixels["distance"] = pixels["bin2_id"] - pixels["bin1_id"]
    d_sum = pixels.groupby("distance")["count"].sum()
    distance_expected = d_sum / (chr_len - d_sum.index + 1)
    OE = pixels.apply(
        lambda x: x["count"]/distance_expected[x.distance], axis=1).astype('float32')
    return OE


def __cal_d_weight(mask, px):
    mask_bins = list(mask.index - px.bin1_id.min())
    mask_index = np.zeros((max(mask_bins)+1, 1))
    mask_index[mask_bins] = 1
    mask_mat = np.dot(mask_index, mask_index.T)
    d_sum = {}
    for i in range(len(mask_index)):
        mask_sum = np.diag(mask_mat, i).sum()
    d_sum[i] = mask_sum
    d = px.groupby("distance")["count"].sum()
    d = pd.DataFrame(d)
    d["distance"] = d.index
    d["sum"] = d["distance"].apply(lambda x: d_sum.get(x, 0))
    d["weight"] = d["sum"] / d["count"]
    return d["weight"].to_dict()


def _sparse_oe_with_gap(pixels, cut_off=0.1):
    # calculate weight using matrix(chr_len,chr_len), pay attention to your memory
    chr_len = pixels.bin2_id.max() - pixels.bin1_id.min()
    # remove gap which bin sum less than chr_len * cut_off
    num_vector = pixels["bin2_id"].value_counts(
    ) + pixels["bin1_id"].value_counts()
    mask = num_vector[num_vector > chr_len * cut_off]
    pixels = pixels[pixels.bin1_id.isin(
        mask.index) & pixels.bin2_id.isin(mask.index)]
    # calculate OE weight base on distance considering gap
    pixels["distance"] = pixels["bin2_id"] - pixels["bin1_id"]
    d_weight = __cal_d_weight(mask, pixels)
    # calculate OE use count * distance_weight
    OE = pixels.apply(
        lambda x: x["count"] * d_weight.get(x["distance"], 0), axis=1).astype('float32')
    return OE


def fast_oe(path,
            chrom=None,
            method='intra',
            gap=False,
            store=False,
            newpath=None,
            verbose=False):
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
    if gap:
        cal_oe = _sparse_oe_with_gap
    else:
        cal_oe = _sparse_oe
    if chrom is not None:
        # Apply OE balancing for one chromsome.
        pixels = clr.matrix(as_pixels=True, ignore_index=False,
                            balance=False).fetch(chrom, chrom)
        p = cal_oe(pixels)
        OE[p.index] = p.values
        print("chrom", chrom, "OE balance finished.") if verbose else 0
    elif method == 'intra' or 'intra-only':
        # Apply OE balancing only for each chromsome internal.
        for i in chroms:
            try:
                pixels = clr.matrix(
                    as_pixels=True, ignore_index=False, balance=False).fetch(i, i)
                p = cal_oe(pixels)
                OE[p.index] = p.values
                print("chrom", i, "OE balance finished.") if verbose else 0
            except:
                print("chrom", i, "OE balance failed, Skipping.") if verbose else 0
    elif method == "all":
        # Apply OE balancing for all chromsomes, including inter-chromsome contacts between two chromsome. Time comsuming, pay attention.
        for i in chroms:
            for j in chroms:
                pixels = clr.matrix(
                    as_pixels=True, ignore_index=False, balance=False).fetch(i, j)
                if len(pixels) > 0:
                    p = cal_oe(pixels)
                    OE[p.index] = p.values
    else:
        print('Unrecognize OE balancing method, apply genome wide global OE balancing and do not distinguish chromosomes.')
        pixels = clr.pixels()[:]
        OE = cal_oe(pixels)
    OE = pd.DataFrame(OE, columns=["OE"])
    if store:
        if newpath is not None:
            shutil.copyfile(path, newpath)
            path = newpath
            clr = cooler.Cooler(path)
        with clr.open("r+") as grp:
            delete(grp["pixels"], "OE")
            put(grp["pixels"], OE)
        print("Writing OE normalized matrix to ", path)
    return pd.concat([clr.pixels()[:], OE], axis=1)
# read OE normalized matrix please use cooler.Cooler(path).matrix(field="OE",balance=False)


def __compartment_decomposition(path,
                              chrom=None,
                              gc_path=None,
                              genome='mm10',
                              vec="E1",
                              fasta='~/ref/mm10.fa',
                              return_vector=True,
                              store=False):
    """
    compartment_adj Calculate compartment adj matrix by cooler.eigs_cis

    Parameters
    ----------
    path : str
        path of cool file
    chrom : str, optional
        A chromsome name or chromsome, by default "chr1", if None calculate all the chromsomes.
    gc_path : str, optional
        Path to gc_cov file, create if not exist, by default None
    genome : str, optional
        genome prefix of generated gc_cov file, by default "mm10"
    vec : str, optional
        Choose PCA EigVec from ["E1","E2","E3"], by default "E1"
    fasta : str, optional
        if genome does not exist, create gc_cov from fasta file, by default '/store/zlwang/ref/mm10.fa'
    return_vector : bool, optional
        return adj matrix or vector, by default True

    Returns
    -------
    np.matrix or np.array or pd.DataFrame
        A matrix or array of compartment values.
    """
    clr = cooler.Cooler(path)
    if gc_path is None:
        gc_path = f'{genome}_gc_cov_{clr.binsize}.tsv'
    try:
        gc_cov = pd.read_table(gc_path)
        if (chrom not in set(gc_cov["chrom"])) and not set(clr.chromnames).issubset(set(gc_cov["chrom"])):
            raise Exception(f"Chromsome {chrom} not in {gc_path}")
    except:
        import bioframe
        print(
            f"Can't find file {gc_path} in current path, try to generate GC coverage from fasta file by bioframe...")
        bins = clr.bins()[:]
        try:
            gc_cov = bioframe.frac_gc(
                bins[['chrom', 'start', 'end']], bioframe.load_fasta(fasta))
        except:
            bins = bins[bins["chrom"] == chrom]
            gc_cov = bioframe.frac_gc(
                bins[['chrom', 'start', 'end']], bioframe.load_fasta(fasta))
        gc_cov.to_csv(gc_path, index=False, sep='\t')
        print(
            f"Generating {genome} GC covergae file {gc_path} from {fasta} successfully. ")
    if chrom is not None:
        view_df = pd.DataFrame({'chrom': clr.chromnames, 'start': 0,
                               'end': clr.chromsizes.values, 'name': clr.chromnames})
        view_df = view_df[view_df.chrom == chrom].reset_index(drop=True)
        comp = eigs_cis(clr, phasing_track=gc_cov, view_df=view_df)[1]
        comp = comp[comp.chrom == chrom]
    else:
        comp = eigs_cis(clr, phasing_track=gc_cov)[1]
    if store and (chrom is None):
        with clr.open("r+") as grp:
            put(grp["bins"], comp['E1'])
            put(grp["bins"], comp['E2'])
            put(grp["bins"], comp['E3'])
        print("Writing compartments score to ", path)
    Eig = comp[vec].values
    if return_vector:
        return Eig
    # if return_matrix:
    #     Eig[Eig > 1] = 1
    #     Eig[Eig < -1] = -1
    #     return np.outer(Eig, Eig)
    return comp


def __tad_insulation(path,
                   chrom=None,
                   window_size=10,
                   return_vector=True,
                   store=False):
    """
    tad_adj calculate insulation score for cooler file.

    Parameters
    ----------
    path : str
        path of .cool file
    chrom : str, optional
        A chromosome name, if None calculate all chroms, by default None
    window_size : int, optional
        window size to calculate insulation score, by default 10
    return_vector : bool, optional
        whether return a vector or matrix , by default True
    store : bool, optional
        whether store insulation score and boundary in .cool file, by default False

    Returns
    -------
    np.array or np.matrix or pd.DataFrame
        A matrix or array represent insulation score of chromsome in .cool file.
    """
    clr = cooler.Cooler(path)
    res = clr.binsize
    if chrom is not None:
        view_df = pd.DataFrame({'chrom': clr.chromnames, 'start': 0,
                               'end': clr.chromsizes.values, 'name': clr.chromnames})
        view_df = view_df[view_df.chrom == chrom].reset_index(drop=True)
        ins = insulation(clr, window_bp=[res*window_size], view_df=view_df)
        ins = ins[ins.chrom == chrom]
    else:
        ins = insulation(clr, window_bp=[res*window_size])
    score = ins["log2_insulation_score_"+str(res*window_size)]
    boundary_strength = ins["boundary_strength_"+str(res*window_size)]
    boundary = ins["is_boundary_"+str(res*window_size)]
    if store and (chrom is None):
        with clr.open("r+") as grp:
            put(grp["bins"], score)
            put(grp["bins"], boundary)
            put(grp["bins"], boundary_strength)
        print("Writing insulation score and boundary to ", path)
    boundary = score[ins["is_boundary_"+str(res*window_size)]]
    score = score.values
    if return_vector:
        return score
    # boundary_adj = np.zeros((len(score), len(score)))
    # if return_matrix:
    #     score[score > 1] = 1
    #     score[score < -1] = -1
    #     for i in range(len(boundary)-1):
    #         start, end = boundary.index[i], boundary.index[i+1]
    #         boundary_adj[start:end, start:end] = 1
    #     return boundary_adj
    return ins
