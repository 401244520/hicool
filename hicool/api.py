from typing import *
import h5py
import pandas as pd
import cooler

class HiCool():
    """
    HiCool _summary_
    """
    
    def __init__(
        self,
        path :str = None,
        meta = {},
        embedding = {},
        network = {},
        imputation = {},
        
    ):
        self.root = path
        if self.is_hicool_file(path):
            print("Successful loading hicool " + path)
            self.scool_path = path + "::scool"
            self.is_hicool = True
            self.meta = self._get("meta")
            self.embedding = self._get("embedding")
            self.network = self._get("network")
            self.imputation = self._get("imputation")
        elif cooler.fileops.is_scool_file(path):
            print("Loading scool file from " + path)
            self.scool_path = path
            self.is_hicool = False
            self.meta = meta
            self.embedding = embedding
            self.network = network
            self.imputation = imputation
        else:
            raise ValueError("path has to be a hicool or scool file, see transfrom.txt2scool or cooler.create_scool to generate scool file.")

    def is_hicool_file(self, path):
        """
        Check if the input file is a HiCool file.
        """
        with h5py.File(path ,"r") as f:
            if f.attrs.get("format") == "HDF5::HICOOL":
                return True
            else : 
                return False
    
    def save_as(self, output_path=None):
        """
        Write the HiCool object to a .hicool file which is hdf5 format.
        """
        if output_path is not None:
            self.root = output_path
        if not self.root.endswith(".hicool"):
            self.root += ".hicool"
        from cooler.fileops import cp
        cp(self.scool_path,self.root+"::scool",overwrite=True)
        with h5py.File(self.root ,"r+") as f:
            import time
            f.attrs.create("format","HDF5::HICOOL")
            f.attrs.create("format-url","https://github.com/401244520/hicool")
            f.attrs.create("creation-date",time.strftime("%Y-%m-%dT%H:%M:%S"))
            f.attrs.create("bin-size",cooler.Cooler(self.scool_path).info["bin-size"])
            f.attrs.create("nbins",cooler.Cooler(self.scool_path).info["nbins"])
            f.attrs.create("ncells",cooler.Cooler(self.scool_path).info["ncells"])
            f.attrs.create("nchroms",cooler.Cooler(self.scool_path).info["nchroms"])
            f.attrs.create("embedding",",".join(self.embedding.keys()))
            f.attrs.create("network",",".join((self.network.keys())))
            f.attrs.create("imputation",",".join(self.imputation.keys()))
            #f.attrs.create("scool",self.root + "::scool")
        self._put("meta",self.meta)
        self._put("embedding",self.embedding)
        self._put("network",self.network)
        self._put("imputation",self.imputation)
        print("Hicool has been saved to " + self.root)

    def info(self):
        """
        File and user metadata dict.
        """
        d = {}
        with h5py.File(self.root ,"r") as f:
            for k, v in f.attrs.items():
                d[k] = v
            return d

    def _get(self,grp):
        d = {}
        with h5py.File(self.root ,"r") as f:
            for k,v in f[grp].items():
                d[k] = v[:]
            return d

    def _put(self,grp,d):
        with h5py.File(self.root ,"r+") as f:
            if grp not in f:
                f.create_group(grp)
            for k,v in d.items():
                try:
                    f.create_dataset(grp+"/"+k,data=v)
                except:
                    cooler.core.delete(f,grp+"/"+k)
                    f.create_dataset(grp+"/"+k,data=v)
    
    def load_hicool_cells(self) -> List:  
        with h5py.File(self.root,'r') as f :
            cells = list(f['scool/cells'].keys())
        cell_list = [self.root + "::/scool/cells/" + cell for cell in cells]
        return cell_list
    
    def to_scanpy(self,embedding_name,embedding_index=None):
        try:
            from scanpy import AnnData
        except:
            print("Please install scanpy using `pip install scanpy`.")
        mat = self.embedding[embedding_name]
        if embedding_index is None:
            embedding_index = range(mat.shape[1])
        var = pd.DataFrame(embedding_index,index=embedding_index,columns=[embedding_name])
        obs = pd.DataFrame(self.meta)
        sce = AnnData(mat, obs = obs, var = var )
        return sce
    
    def to_higashi(self):
        pass
    
    def to_fast_higashi(self):
        pass
