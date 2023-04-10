import cooler
from tqdm import tqdm
import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial
from .process import AutoLoad
from .api import HiCool

def compress_matrix(hc : HiCool,
                    chrom : str = None,
                    operation : str = "bin_strength",
                    field : str = "count",
                    nproc = 10,
                    regions : list = None,
                    save_embedding = False,
                    **kwargs)  -> np.array:
    """
    compress_matrix 
    

    Parameters
    ----------
    hc : HiCool
        A .hicool file including scool file, reference hicool.api.HiCool
    chrom : str, optional
        Choose one chromsome or all the chromsomes included in scool, by default None
    operation : str, optional
        matrix operation for single cell Hi-C data, by default "bin_strength"
        reference hicool.function.decomposition.pixels_operation, or you can customize your own method 
    field : str, optional
        which columns to choose in cool.pixels , by default "count"
    nproc : int, optinal
        number of process ,by default 10
    save_embedding : bool, optional
        whether save embedding result in hicool.embedding dict ,by default False

    Returns
    -------
    numpy.ndarray
        A 2D matrix which rows represent each cell, columns represent compressed features, shape is cell*features
    list 
        A list contain the feature index 
    """    
    from .function.decomposition import pixels_operation
    if hc.is_hicool:
        cell_list = AutoLoad(hc.root).load_hicool_cells()
    else:
        cell_list = AutoLoad(hc.scool_path).load_scool_cells()
    chrom_list = cooler.Cooler(hc.scool_path).chromnames
    if regions is not None:
        clr = cooler.Cooler(cell_list[0])
        mask = []
        for region in regions:
            try:
                clr.extent(region)
                mask.append(True)
            except:
                mask.append(False)
        valid_regions =  regions[mask]
        cal_bins = partial(pixels_operation,regions=valid_regions,field=field,operation=operation,**kwargs)
        with Pool(processes=nproc) as pool :
            bin_strength = list(tqdm(pool.imap(cal_bins,cell_list), total= len(cell_list)))
        bins_df = pd.concat(bin_strength,axis=1)
        bins_matrix = bins_df.T.fillna(0)
    elif chrom in ["all","ALL",None,""]:
        # calculate all bins and concat
        all_bins = []
        for chrom in tqdm(chrom_list) :
            cal_bins = partial(pixels_operation,chrom=chrom,field=field,operation=operation,**kwargs)
            with Pool(processes=nproc) as pool :
                bin_strength = list(pool.imap(cal_bins,cell_list))
            bs_df = pd.concat(bin_strength,axis=1)#.sort_index()
            bs_count = bs_df.T.fillna(0)
            all_bins.append(bs_count)
        bins_matrix = pd.concat(all_bins,axis=1)
    elif chrom in chrom_list:
        cal_bins = partial(pixels_operation,chrom=chrom,field=field,operation=operation,**kwargs)
        with Pool(processes=nproc) as pool :
            bin_strength = list(pool.imap(cal_bins,cell_list))
        bs_df = pd.concat(bin_strength,axis=1)#.sort_index()
        bins_matrix = bs_df.T.fillna(0)
    else:
        print(f"chromsome {chrom} not in list{chrom_list} !")
    if save_embedding:
        hc.embedding[operation] = bins_matrix.to_numpy()
    return bins_matrix.to_numpy(),bins_matrix.columns.to_list()




def fetch_matrix(hc:HiCool,
                 region:str,
                 normlization:str=None,
                 field:str="count",
                 nproc = 10,
                 ) -> np.array:
    """Fetch genome matrix from hicool files.

    Parameters
    ----------
    hc : HiCool
        A HiCool object.
    region : str
        A genome region such as "chr1:1000000-2000000" or whole chromsome like "chr1".
    normlization : str, optional
        Choose a normaliztion method in ["OE","KR","VC"], by default None
    field : str, optional
        which columns to choose in cool.pixels , by default "count"
    nproc : int, optinal
        number of process ,by default 10
    """
    from .function.dataloader import load_scool_region
    matrix = load_scool_region(hc.scool_path,region,as_pixels=False)
    return matrix
