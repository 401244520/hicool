from functools import partial
import os,glob,time
import pandas as pd
import numpy as np
import cooler
from tqdm import tqdm
from multiprocessing import Pool

from .loading import AutoLoad

def txt2cool(file_path,chroms,positions,resolution,chrom_list,chromsizes):
    file = AutoLoad(file_path)
    file.load_txt(chroms = chroms,positions = positions,header = None ,comment = "#", low_memory = False)
    file.resolution = resolution
    bins,pixels = file.txt2cool(chrom_list,chromsizes)
    return {file.name:pixels}

def txt2scool(file_path,chrom_list,chromsizes,prefix = 'Test',nproc = 20,chroms = [1,3],positions = [2,4],resolution = 1000000,save_dir = "."):
    """
    txt2scool Merge thousands of single txt raw conatacts files into one certain resolution scool file.

    Parameters
    ----------
    file_path : list, iterable
        A list or list like(np.array) of all txt files path.
    chrom_list : list, iterable
        A list or list like of all chroms reserved in scool file. 
        Eg.['chr1','chr2'....'chrX']
    chromsizes : str
        File path of chromsizes in txt format, download at UCSC. eg: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
    prefix : str, optional
        scool name prefix, by default 'Test'
    nproc : int, optional
        number of the process , by default 20
    chroms : list, optional
        the chroms of the contacts pairs in txt files, by default [1,3]
    positions : list, optional
        the position of the contacts pairs in txt file, by default [2,4]
        Eg. set chroms = [1,3] and position = [2,4] for this condition below:
        #columns: readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1
        .       chr1    3000195 chr1    31532740        +       -       1       .
    resolution : int, optional
        the resolution of scool file, by default 1000000

    Returns
    -------
    str
        Create a scool file named by prefix and resolution and return absolute path. 
    """    
    print(" Start processing files ...")
    file = AutoLoad(file_path[0])
    chroms = [int(i) for i in chroms]
    positions = [int(i) for i in positions]
    t1 = time.time()
    # file.load_txt(chroms = chroms,positions = positions,header = None ,comment = "#", low_memory = False)
    # file.resolution = resolution
    # bins,pixels = file.txt2cool(chrom_list,chromsizes)
    # t2 = time.time() - t1
    # l = sum([os.path.getsize(i) for i in file_path])/os.path.getsize(file_path[0])
    # eval_time = int(t2 * l / nproc ) * 10 + 20
    # print(f" Total {len(file_path)} files processing with {nproc} proc at {resolution} resolution, evaluate finish in {eval_time} seconds. ")
    bins = file.txt2cool(chrom_list,chromsizes,return_bins=True)
    
    with Pool(processes = nproc) as pool:
        result = list(tqdm(pool.imap(partial(txt2cool,chroms=chroms,positions=positions,resolution=resolution,chrom_list=chrom_list,chromsizes=chromsizes),file_path), total= len(file_path)))
    name_file_list = {}
    [name_file_list.update(i) for i in result ]
    
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    save_path = f"{save_dir}/{prefix}_{resolution}.scool"
    
    cooler.create_scool(save_path,bins,name_file_list,ordered = True,ensure_sorted=True)
    print(f' Creating {save_path} finished , including {len(result)} single cells with resolution {resolution}. ')
    
    t3 = time.time() - t1
    print(f" Totally take {t3} seconds.")
    return os.path.abspath(save_path)

def scool2cool(scool_path):
    cells = cooler.fileops.list_scool_cells(scool_path)
    cells_path = [scool_path + "::" + cell for cell in cells]
    merged_path = scool_path.replace(".scool","_merged.cool")
    print(f"Start merging scool files in {scool_path} ... " )
    cooler.merge_coolers(merged_path,cells_path,mergebuf=10e10)
    print(f"Finish create merged cool file at the same path of scool file. Start balancing ... " )
    merged_cool = cooler.Cooler(merged_path)
    cooler.balance_cooler(merged_cool,store=True)
    print(f"Finish balance and create merged cool file at {merged_path}. " )


