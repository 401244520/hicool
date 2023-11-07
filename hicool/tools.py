import cooler
from tqdm import tqdm
import pandas as pd
import numpy as np
from typing import Union
from multiprocessing import Pool
from functools import partial
from .api import HiCool


def validate_regions(regions, chromsize):
    mask = []
    for region in regions:
        chrom, start, end = cooler.util.parse_region(region)
        if (chrom in chromsize.index) :
            valid = (start >= 0) & (end <= chromsize[chrom])
        else:
            valid = False
        mask.append(valid)
    valid_regions = regions[mask]
    if sum(~mask) > 0:
        print(f"{sum(mask)} regions is valid while {sum(~mask)} not in hicool:{regions[~mask]}")
    return valid_regions

def calculate_embedding(chrom, cell_list, field, operation, nproc, kwargs):
    from .function.embedding import calculate_embedding
    cal_embs = partial(calculate_embedding, chrom=chrom, field=field, operation=operation, **kwargs)
    with Pool(processes=nproc) as pool:
        bin_emb = list(pool.imap(cal_embs, cell_list))
    df = pd.concat(bin_emb, axis=1)
    return df.T.fillna(0)

def calculate_feature(chrom, cell_list, operation, nproc, kwargs):
    from .function.conformation import calculate_feature
    cal_feat = partial(calculate_feature, chrom=chrom, operation=operation, **kwargs)
    with Pool(processes=nproc) as pool:
        bin_feat = list(pool.imap(cal_feat, cell_list))
    df = pd.concat(bin_feat, axis=1)
    return df.T.fillna(0)

def embedding_hicool(hc: HiCool,
                     regions: list = None,
                     chrom: Union[str,list] = None,
                     operation: str = "bin_strength",
                     field: str = "count",
                     nproc: int = 8,
                     embedding_name: str = None,
                     **kwargs) -> Union[np.matrix, list]:
    """
    This function compresses a Hi-C matrix into lower-dimensional embedding.
    
    Parameters
    ----------
    hc : HiCool
        A .hicool file including scool file, reference hicool.api.HiCool
    regions : list, optional
        A list including the regions to calculate. If specified, calculation will be applied to these regions only. Default is None.
    chrom : str, optional
        If regions is None, apply calculation on a single chromosome such as 'chr1'. 
        If chrom is ["all", "ALL", None, ""], apply calculation on all the chromosomes. Default is None.
    operation : str, optional
        Matrix operation for single cell Hi-C data. Default is "bin_strength".
        Reference hicool.function.embedding.pixels_embedding, or you can customize your own method. 
    field : str, optional
        Which columns to choose in cool.pixels. Default is "count".
    nproc : int, optinal
        Number of processes. Default is 8.
    embedding_name : str, optional
        Whether save embedding result in embedding dict. Default is None.
        If not None, store embedding in hicool/embedding/embedding_name and index in hicool/embedding_index/embedding_name.
        
    Returns
    -------
    numpy.ndarray
        2D matrix which rows represent each cell, columns represent compressed features,
        shape is cell*features(bins).
    list 
        Containing the feature index. 
    """

    cell_list = hc.load_cells()
    chromsize = cooler.Cooler(hc.scool_path).chromsizes
    valid_chrom = list(chromsize.index)
    # If regions is not None, calculate embedding via regions
    if regions is not None:
        valid_regions = validate_regions(regions, chromsize)
        bins_matrix = calculate_embedding(valid_regions, cell_list, field, operation, nproc, kwargs)
    # Default calculate all chroms
    elif chrom in ["all", "ALL", None, ""]:
        all_bins = [calculate_embedding(chrom, cell_list, field, operation, nproc, kwargs) for chrom in tqdm(valid_chrom)]
        bins_matrix = pd.concat(all_bins, axis=1)
    # Also support chromsome list
    elif isinstance(chrom, list):
        chrom_bins = []
        for chrom in tqdm(chrom):
            if chrom in valid_chrom:
                chrom_bins.append(calculate_embedding(chrom, cell_list, field, operation, nproc, kwargs))
            else:
                print(f"chromsome {chrom} not in list {valid_chrom}, skipping !")
        bins_matrix = pd.concat(chrom_bins, axis=1)
    # Calculate single chromsome
    elif chrom in valid_chrom:
        bins_matrix = calculate_embedding(chrom, cell_list, field, operation, nproc, kwargs)
    else:
        print(f"chromsome {chrom} not in list {valid_chrom} !")
    # Convert the bin matrix to a numpy array and get the column names
    embedding_matrix = bins_matrix.to_numpy()
    embedding_index = bins_matrix.columns.to_list()
    # Check if an embedding name is provided, save embedding to hicool if True
    if embedding_name:
        # Initialize 'embedding' and 'embedding_index' dictionaries if not already present
        if 'embedding' not in hc.groups:
            hc.groups['embedding'] = {}
        if 'embedding_index' not in hc.groups:
            hc.groups['embedding_index'] = {}
        # Add the new embedding and index to the 'embedding' and 'embedding_index' dictionaries
        hc.groups['embedding'][embedding_name] = embedding_matrix
        hc.groups['embedding_index'][embedding_name] = embedding_index
        print(f"Added embedding '{embedding_name}' to 'hicool/embedding'. Current embeddings: {list(hc.groups['embedding'].keys())}")
    return embedding_matrix, embedding_index


def features_hicool(hc: HiCool,
                    operation: str = "compartment",
                    chrom: str = None,
                    nproc: int = 8,
                    feature_name: str = None,
                    **kwargs) -> Union[np.matrix, list]:
    cell_list = hc.load_cells()
    valid_chrom = cooler.Cooler(hc.scool_path).chromnames
    if chrom in ["all", "ALL", None, ""]:
        all_bins = [calculate_feature(chrom, cell_list, operation, nproc, kwargs) for chrom in tqdm(valid_chrom)]
        bins_matrix = pd.concat(all_bins, axis=1)
    elif chrom in valid_chrom:
        bins_matrix = calculate_feature(chrom, cell_list, operation, nproc, kwargs)
    else:
        print(f"chromsome {chrom} not in list {valid_chrom} !")
    matrix = bins_matrix.to_numpy()
    index = bins_matrix.columns.to_list()
    if feature_name:
        # Initialize 'embedding' and 'index' dictionaries if not already present
        if 'feature' not in hc.groups:
            hc.groups['feature'] = {}
        if 'feature_index' not in hc.groups:
            hc.groups['feature_index'] = {}
        # Add the new embedding and index to the 'embedding' and 'index' dictionaries
        hc.groups['feature'][feature_name] = matrix
        hc.groups['feature_index'][feature_name] = index
        print(f"Added embedding '{feature_name}' to 'hicool/feature'. Current embeddings: {list(hc.groups['feature'].keys())}")
    return matrix, index