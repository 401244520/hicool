import cooler
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from functools import partial
from multiprocessing import Pool

try : 
    import torch
    from torch.utils import data
    to_tensor = True
except:
    print("Please install pytorch if you want to load data to_tensor. Now using multiprocessing.Pool.")
    to_tensor = False

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
                #   transform: str = None,
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
    # if transform in ["OE","oe","observed/expected","O/E"]:
    #     oe = _sparse_oe(pixels)
    #     pixels["count"] = oe
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
                     label_list,
                     chrom='chr1',
                     sparse=False,
                     field='count',
                    #  transform=None
                     ) :
            self.cell_list = [scool_path + "::" + c for c in cooler.fileops.list_scool_cells(scool_path)]
            self.chrom = chrom
            # self.transform = transform
            self.sparse = sparse
            self.field = field
            self.label_list = label_list if label_list is not None else []
            
        def __getitem__(self, index):
            cell = self.cell_list[index]
            data_tensor = load_cool_adj(cell,self.chrom,self.sparse,self.field)#,self.transform)
            if self.label_list:
                label = self.label_list[index]
                return data_tensor, label
            else:
                return data_tensor
            
        def __len__(self):
            return len(self.cell_list)

    def load_scool_adj(scool_path,
                       chrom='chr1',
                       sparse=False,
                       field='count',
                    #    transform=None,
                       return_loader = False):
        scool_data = CoolData(scool_path,chrom,sparse,field)#,transform)
        dataloader = data.DataLoader(scool_data,batch_size=len(scool_data),num_workers=nproc)
        if return_loader:
            return dataloader
        return iter(dataloader).next()
else : 
    def load_scool_adj(scool_path,
                       chrom='chr1',
                       sparse=False,
                       field='count',
                    #    transform=None
                       ):
        cell_list = [scool_path + "::" + c for c in cooler.fileops.list_scool_cells(scool_path)]
        load_adj = partial(load_cool_adj,chrom=chrom,sparse=sparse,field=field,to_tensor=False)#,transform=transform)
        with Pool(processes=nproc) as pool:
            result =  list(pool.imap(load_adj,cell_list))
        return np.array(result)

import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np

def load_cool_region(cool_path,
                     region,
                     as_pixels=True,
                     as_graph=False,
                     field="count"):
    clr = cooler.Cooler(cool_path)
    if as_graph:
        result = clr.pixels().fetch(region)
    else:
        result = clr.matrix(as_pixels=as_pixels,balance=False,field=field).fetch(region)
    return result

def load_scool_region(scool_path:str,
                      region:str,
                      as_pixels:bool=True,
                      as_graph:bool=False,
                      field:str="count",
                      nproc:int = 10
                      ):
    """load matrix or graph from scool file.

    Parameters
    ----------
    scool_path : str
        path to scool file.
    region : str
        A genome range you want to load, such as "chr1:10000-20000".
    as_pixels : bool, optional
        Whether convert matrix to graph within the region, by default True
    as_graph : bool, optional
        Whether converse all the contacts to the region, by default False
    field : str, optional
        which column to choose , by default "count"

    Returns
    -------
    numpy.array
        If as_pixels and as_graph is True,a 3D array including all cells matrices, shape is cells*region_len*region_len.
    pandas.DataFrame
        If as_pixels or as_graph is True, return a Dataframe including contact pixels and cell id.
    """
    cell_list = [scool_path + "::" + c for c in cooler.fileops.list_scool_cells(scool_path)]
    load_region = partial(load_cool_region,region=region,as_pixels=as_pixels,field=field)
    with Pool(processes=nproc) as pool:
        results =  list(pool.imap(load_region,cell_list))
    if as_pixels or as_graph :
        for i,result in enumerate(results):
            result["cell"] = i
        hyper_pixels = pd.concat(results).reset_index(drop=True)
        return hyper_pixels
    else :
        hyper_matrix = np.array(results)
        return hyper_matrix



def train_val_test_split(dataset, 
                         test_ratio=0.1, 
                         val_ratio=0, 
                         batch_size=32, 
                         shuffle=True, 
                         random_seed=42):
    """
    train_val_test_split : Split a dataset into train, validation, and test sets.

    Parameters
    ----------
    dataset : Dataset
        The dataset to be split.
    test_ratio : float, optional
        The ratio of test set to the whole dataset, by default 0.1
    val_ratio : float, optional
        The ratio of validation set to the whole dataset, by default 0.1
    batch_size : int, optional
        The batch size for the data loaders, by default 32
    shuffle : bool, optional
        Whether to shuffle the indices before splitting, by default True
    random_seed : int, optional
        The random seed for shuffling, by default 42

    Returns
    -------
    train_loader : DataLoader
        The data loader for the validation set.
    val_loader : DataLoader
        The data loader for the test set.
    test_loader : DataLoader
        The data loader for the test set.
    indices : dict
        Dictionary containing the indices for train, val, and test sets.
    """
    
    # Determine the sizes of each dataset
    num_data = len(dataset)
    val_size = int(num_data * val_ratio)
    test_size = int(num_data * test_ratio)
    train_size = num_data - val_size - test_size

    # Split the dataset into train, validation, and test sets
    indices = np.arange(num_data)
    if shuffle:
        np.random.seed(random_seed)
        np.random.shuffle(indices)
    train_indices = indices[:train_size]
    val_indices = indices[train_size:train_size + val_size]
    test_indices = indices[train_size + val_size:]

    # Create the data loaders for each dataset
    train_loader = DataLoader(dataset, batch_size=batch_size, sampler=torch.utils.data.sampler.SubsetRandomSampler(train_indices))
    val_loader = DataLoader(dataset, batch_size=batch_size, sampler=torch.utils.data.sampler.SubsetRandomSampler(val_indices))
    test_loader = DataLoader(dataset, batch_size=batch_size, sampler=torch.utils.data.sampler.SubsetRandomSampler(test_indices))

    return train_loader, val_loader, test_loader, {'train': train_indices, 'test': test_indices, 'val': val_indices}
