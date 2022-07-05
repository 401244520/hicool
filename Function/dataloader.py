from attr import field
import cooler
import numpy as np
import torch
from scipy.sparse import coo_matrix
from torch.utils import data
from .features import  _fast_oe

def sparse_mx_to_torch_sparse_tensor(sparse_mx): 
    """Convert a scipy sparse matrix to a torch sparse tensor."""
    sparse_mx = sparse_mx.tocoo().astype(np.float32)
    indices = torch.from_numpy(
        np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.int64))
    values = torch.from_numpy(sparse_mx.data)
    shape = torch.Size(sparse_mx.shape)
    return torch.sparse.FloatTensor(indices, values, shape)

def load_cool_adj(cool_path,chrom='chr1',sparse=False,transfrom=None,balance=False,field = 'count'):
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
    transfrom : str, optional
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
    if sparse:
        pixels = clr.matrix(as_pixels=True,balance=balance).fetch(chrom,chrom)
        if transfrom in ["OE","oe","observed/expected","O/E"]:
            oe = _fast_oe(pixels)
            pixels["count"] = oe
        row = pixels.bin1_id - clr.offset(chrom)
        col = pixels.bin2_id - clr.offset(chrom)
        adj = coo_matrix((pixels['count'], (row, col)),shape=(size,size))
        adj = sparse_mx_to_torch_sparse_tensor(adj)            
    else:
        pixels = clr.matrix(balance=balance,field=field).fetch(chrom,chrom)
        adj = torch.FloatTensor(pixels)
    return adj

class CoolData(data.Dataset):
    def __init__(self,scool_path,chrom='chr1',sparse=False,transfrom=None) :
        self.cell_list = [scool_path + "::" + c for c in cooler.fileops.list_scool_cells(scool_path)]
        self.chrom = chrom
        self.transform = transfrom
        self.sparse = sparse
    def __getitem__(self, index) :
        cell = self.cell_list[index]
        data_tensor = load_cool_adj(cell,self.chrom,self.sparse,self.transform)
        return data_tensor
    def __len__(self):
        return len(self.cell_list)

def load_scool_adj(scool_path,chrom='chr1',sparse=False,transfrom=None):
    scool_data = CoolData(scool_path,chrom,sparse,transfrom)
    dataloader = data.DataLoader(scool_data,batch_size=len(scool_data),num_workers=20)
    return iter(dataloader).next()
