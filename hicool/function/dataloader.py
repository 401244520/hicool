import cooler
import numpy as np
from typing import Union
from scipy.sparse import coo_matrix
from cooler.core import delete, get, put

try : 
    import torch
    from torch.utils import data
    to_tensor = True
except:
    from functools import partial
    from multiprocessing import Pool
    print("Please install pytorch if you want to load data to_tensor. Now using multiprocessing.Pool.")
    to_tensor = False

from .features import _sparse_oe

nproc = 20

def sparse_mx_to_torch_sparse_tensor(sparse_mx): 
    """Convert a scipy sparse matrix to a torch sparse tensor."""
    sparse_mx = sparse_mx.tocoo().astype(np.float32)
    indices = torch.from_numpy(
        np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.int64))
    values = torch.from_numpy(sparse_mx.data)
    shape = torch.Size(sparse_mx.shape)
    return torch.sparse.FloatTensor(indices, values, shape)

def load_cool_adj(cool_path : str,
                  chrom: str = 'chr1',
                  sparse: bool = False,
                  field: str = 'count',
                  transform: str = None,
                  balance: bool = False,
                  to_tensor: bool = True):
    """
    load_cool_adj : Generate HiC contact adjacency matrix from cool file.

    Parameters
    ----------
    cool_path : str
        url path of cool file
    chrom : str, optional
        chromsome name, by default 'chr1'
    sparse : bool, optional
        if True return coo_matrix, by default False
    transform : str, optional
        A method in ["OE" ....] for matrix transform, by default None
    balance : bool, optional
        if True using "weight" to balance cooler, need cooler.balance_cooler to preprocess it, by default False
    field : str, optional
        fetch the columns from cooler.pixels dataframe, need preprocess cooler, by default 'count' or None

    Returns
    -------
    torch.FloatTensor or torch.sparse.FloatTensor
        A sparse or dense matrix of certain chromsome from cool files.
    """
    clr = cooler.Cooler(cool_path)
    size = int(clr.chromsizes[chrom] / clr.binsize)+1
    pixels = clr.matrix(as_pixels=True,balance=balance,field=field).fetch(chrom,chrom)
    if transform in ["OE","oe","observed/expected","O/E"]:
        oe = _sparse_oe(pixels)
        pixels["count"] = oe
    row = pixels.bin1_id - clr.offset(chrom)
    col = pixels.bin2_id - clr.offset(chrom)
    adj = coo_matrix((pixels[field], (row, col)),shape=(size,size))
    if not sparse:
        adj = adj.todense()
        adj = adj + np.triu(adj,1).T
        if to_tensor:
            adj = torch.FloatTensor(adj)
    else :
        if to_tensor:
            adj = sparse_mx_to_torch_sparse_tensor(adj)      
    return adj

"""def save_cool_adj(adj: Union[coo_matrix,np.matrix],
                  cool_path : str,
                  chrom: str = 'chr1',
                  field: str = 'count',
                  overwrite: bool = True):
    clr = cooler.Cooler(cool_path)
    clr.matrix
    if sparse :
        with clr.open("r+") as grp:
            delete(grp["pixels"],field) if overwrite else 0
            put(grp["pixels"],field)
    return cool_path
"""
if to_tensor :
    class CoolData(data.Dataset):
        def __init__(self,
                     scool_path,
                     chrom='chr1',
                     sparse=False,
                     field='count',
                     transform=None) :
            self.cell_list = [scool_path + "::" + c for c in cooler.fileops.list_scool_cells(scool_path)]
            self.chrom = chrom
            self.transform = transform
            self.sparse = sparse
            self.field = field
        def __getitem__(self, index) :
            cell = self.cell_list[index]
            data_tensor = load_cool_adj(cell,self.chrom,self.sparse,self.field,self.transform)
            return data_tensor
        def __len__(self):
            return len(self.cell_list)

    def load_scool_adj(scool_path,
                       chrom='chr1',
                       sparse=False,
                       field='count',
                       transform=None,
                       return_loader = False):
        scool_data = CoolData(scool_path,chrom,sparse,field,transform)
        dataloader = data.DataLoader(scool_data,batch_size=len(scool_data),num_workers=nproc)
        if return_loader:
            return dataloader
        return iter(dataloader).next()
else : 
    def load_scool_adj(scool_path,
                       chrom='chr1',
                       sparse=False,
                       field='count',
                       transform=None):
        cell_list = [scool_path + "::" + c for c in cooler.fileops.list_scool_cells(scool_path)]
        load_adj = partial(load_cool_adj,chrom=chrom,sparse=sparse,field=field,transform=transform,to_tensor=False)
        with Pool(processes=nproc) as pool:
            result =  list(pool.imap(load_adj,cell_list))
        return np.array(result)
