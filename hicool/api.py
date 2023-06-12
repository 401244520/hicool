from typing import *
import h5py
import pandas as pd
import numpy as np
import cooler

class HiCool():
    """
    HiCool _summary_
    """
    
    def __init__(self,path :str = None,**kwargs):
        self.root = path
        if self.is_hicool_file(path):
            self.scool_path = path + "::scool"
            self.is_hicool = True
            self.groups = {group:self._get(group) for group in self.index-set(['scool'])}
            # self.byte2string(self.meta)
            print(f"Loading hicool {list(self.groups.keys())} from {path} successfully.")
        elif cooler.fileops.is_scool_file(path):
            self.scool_path = path
            self.is_hicool = False
            self.groups = kwargs
            print("Loading scool file from " + path)
        else:
            raise ValueError("path has to be a hicool or scool file, see transfrom.txt2scool or cooler.create_scool to generate scool file.")
    
    def __getitem__(self, item):
        return self.groups[item]
    
    def is_hicool_file(self, path):
        """
        Check if the input file is a HiCool file.
        """
        with h5py.File(path ,"r") as f:
            if f.attrs.get("format") == "HDF5::HICOOL":
                return True
            else : 
                return False
    
    def save(self,overwrite=True):
        """
        Save HiCool object to .hicool file.
        """
        with h5py.File(self.root ,"r+") as f:
            import time
            f.attrs.create("modify-date",time.strftime("%Y-%m-%dT%H:%M:%S"))
            f
        [self._put(group,self.groups[group]) for group in self.groups]
        print("HiCool saved success to " + self.root)
    
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
            #f.attrs.create("scool",self.root + "::scool")
        [self._put(group,self.groups[group]) for group in self.groups]
        print("Hicool has been saved to " + self.root)

    @property
    def info(self):
        """
        List the attributes dict.
        """
        d = {}
        with h5py.File(self.root ,"r") as f:
            for k, v in f.attrs.items():
                d[k] = v
            return d
    
    @property
    def index(self,level=1,root="/"):
        """
        List the groups index.
        # TODO: support more than 3 level.
        """
        with h5py.File(self.root ,"r") as f:
            groups = set()
            f[root].visit(groups.add)
            ind1 = {g.split("/")[0] for g in groups}
        return ind1
            # ind2 = {g.split("/")[0] + "/" + g.split("/")[1] for g in groups-ind1}
            # ind3 = {g.split("/")[0] + "/" + g.split("/")[1]+ "/" + g.split("/")[2]  for g in groups-ind1-ind2}
        # ind = {1:ind1,2:ind2,3:ind3}
        # return ind[level]
    
    def _get(self,grp):
        try:
            d = {}
            with h5py.File(self.root ,"r") as f:
                for k,v in f[grp].items():
                    d[k] = v[:]
                return d
        except:
            print(f"Group {grp} can not slicing, check the data or del it.")

    def _put(self,grp,d):
        with h5py.File(self.root ,"r+") as f:
            if grp not in f:
                f.create_group(grp)
            for k,v in d.items():
                try:
                    f.create_dataset(grp+"/"+k,data=v)
                except:
                    del f[grp+"/"+k]
                    f.create_dataset(grp+"/"+k,data=v)
    
    def load_cells(self) -> List:  
        if self.is_hicool:
            with h5py.File(self.root,'r') as f :
                cells = list(f['scool/cells'].keys())
            cell_list = [self.root + "::/scool/cells/" + cell for cell in cells]
        else:
            with h5py.File(self.root,'r') as f :
                cells = list(f['cells'].keys())
            cell_list = [self.root + "::cells/" + cell for cell in cells]
        return cell_list
    
    def byte2string(self,d):
        c = d.copy()
        for k in d.keys():
            if d[k].dtype == object:
                c[k] = np.array([item.decode() for item in d[k]],dtype=str)
        return c
    
    def to_scanpy(self,embedding_name,embedding_index=None):
        try:
            from scanpy import AnnData
        except:
            print("Please install scanpy using `pip install scanpy`.")
        mat = self.groups["embedding"][embedding_name]
        if embedding_index is None:
            embedding_index = range(mat.shape[1])
        var = pd.DataFrame(embedding_index,index=embedding_index,columns=[embedding_name])
        obs = pd.DataFrame(self["meta"])
        sce = AnnData(mat, obs = obs, var = var )
        return sce
    
    def operation(self,method,store=True):
        cell_list = self.load_cells()
    
    def to_higashi(self):
        pass
    
    def to_fast_higashi(self):
        pass
