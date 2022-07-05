# -*- coding: utf-8 -*-
__author__ = 'zlwang'

import os
import pandas as pd
import numpy as np
import cooler,h5py


## 本模块提供Contact_Matrix这一类，包含交互矩阵，名字，总reads，各染色体的成分等属性
## 并且将三种常见的Hi-C交互数据分别定义为子类进行初始化（格式化）
## Contact_Matrix类目前仅提供了染色体作图和主成分提取的功能，更多拓展功能有待完善


class AutoLoad():
    ''' This a loading methods collection of different Hi-C contact file formats.
        Now support : txt,scool,cool...
        Input: File path
    '''
    def __init__(self, path, *args, **kwargs):
        self.path = os.path.abspath(path)
        self.name = os.path.basename(path).split(".")[0]
        self.resolution = 1000000  # TODO: set resolution and other args
        self.label = None
        self.chroms = {}
        self.contact = None
        self.matrix = None

    def load_txt(self,chroms = [0,2],positions = [1,3],header = None, *args, **kwargs):
        '''  load txt file in 
            params chroms : the columns index position of two contact chromosomes in txt file
            params positions : the columns index position of two contact position on chromosomes in txt file
            params sep,header... : details reference <read_csv>
        '''
        data = pd.read_table(self.path,header=header, *args, **kwargs)
        data = data.iloc[:,chroms + positions]
        data.columns = ['chr1','chr2','pos1','pos2']
        self.contact = data
        self.amount = len(data)
        #print(f'txt data {self.name} loaded... ')
        return data
    
    def load_scool(self):  
        # Suggest use load_scool_cells or load_scool_bins when only need cell list or bins.
        file_path = self.path
        file_list = []
        bins = []
        pixel_list = []
        if cooler.fileops.is_scool_file(file_path): 
            matrices_list = cooler.fileops.list_scool_cells(file_path)
            for cell in matrices_list:
                file_list.append(file_path + "::" + cell)
            for file_name in file_list:
                scool_file = cooler.Cooler(file_name)
                pixels = scool_file.pixels()[:]
                pixel_list.append(pixels)
            bins = scool_file.bins()[:]
            self.chroms = scool_file.chromnames
            self.resolution = scool_file.binsize
            print(f'{self.name} loaded, including {len(file_list)} single cells with resolution {self.resolution}. ')
        else :
            raise FileExistsError(f"Could not find scool file {file_path}, please check the path or format.")
        return file_list,bins,pixel_list
    
    def load_scool_cells(self):  
        with h5py.File(self.path,'r') as f :
            cells = list(f['cells'].keys())
        cell_list = [self.path + "::/cells/" + cell for cell in cells]
        return cell_list
    
    def load_scool_bins(self):  
        with h5py.File(self.path,'r') as f :
            chroms = f['chroms']['name'][:].astype(str)
            keys = [key for key in f['bins'].keys()]
            bins = pd.DataFrame([f['bins'][key][:] for key in keys],index=keys).T
            bins['chrom'] = bins['chrom'].apply(lambda x : chroms[x])
        return bins[['chrom','start','end']]
    
    def load_cool(self):
        cool_file = cooler.Cooler(self.path)
        bins = cool_file.bins()[:]
        pixels = cool_file.pixels()[:]
        return bins,pixels
    
    def txt2cool(self,chroms,chromsize_txt,save_path = None,return_bins = False):
        '''
        transform contact.txt to cool format with certain resolution
        contact_txt : include columns [chr1 chr2 pos1 pos2]. eg: self.conatact
        chroms : chromtain list you want to transform into cool. eg: ["chr1","chr2"]
        res : the resolution of cool. 
        chromsize_txt: ref genome chromsize download at UCSC. eg: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
        save_cool: Boolean. if True, save the cool file as prefix_res.cool, else only return bins and pixels.
        '''
        contact_txt = self.contact
        res = self.resolution
        prefix = self.name
        def _call_bins(chroms,chromsize_txt):
            chromsize = pd.read_table(chromsize_txt,header = None, index_col = 0).loc[chroms,1]
            bins_chrom = []
            bins_start = []
            bins_end = []
            off_set = [0] + ((chromsize/res+1).astype(int)).cumsum().tolist()[:-1]
            off_set = dict(zip(chromsize.index,off_set))
            # make bin index # TODO: assert chrom in chroms
            for chrom in chroms:
                size = int(chromsize[chrom]/res+1)
                bins_chrom += [chrom] * size
                bins_start.append(np.arange(size) * res)
                bins_end.append(np.arange(size) * res + res)
                bins_end[-1][-1] = chromsize[chrom]
            bins = pd.DataFrame(data = [bins_chrom,np.hstack(bins_start),np.hstack(bins_end)]).T
            bins.columns = ["chrom","start","end"]
            return off_set,bins

        def _call_pixels(contact_txt,off_set):
            contacts = contact_txt.copy()
            contacts['count'] = 1
            contacts[["pos1","pos2"]] = (contacts[["pos1","pos2"]] / res).astype(int)
            newcontact = contacts.groupby(["chr1","pos1","chr2","pos2"]).count()
            newcontact = newcontact.reset_index()
            newcontact = newcontact[newcontact.chr1.isin(off_set.keys()) & newcontact.chr2.isin(off_set.keys())]
            newcontact = newcontact.replace(off_set)
            newcontact["bin1_id"] = newcontact.chr1 +newcontact.pos1
            newcontact["bin2_id"] = newcontact.chr2 +newcontact.pos2
            pixels = newcontact[["bin1_id","bin2_id","count"]].reset_index(drop = True)
            # reorder bins is necessary for cooler format # TODO：Time should take into account when matrix is big.
            lowwer_loc = pixels.bin1_id > pixels.bin2_id
            if sum(lowwer_loc) > 0 :
                pixels["bin1_id"].loc[lowwer_loc],pixels["bin2_id"].loc[lowwer_loc] = \
                pixels["bin2_id"].loc[lowwer_loc],pixels["bin1_id"].loc[lowwer_loc]
            pixels = pixels.groupby(["bin1_id","bin2_id"]).sum().reset_index()
            pixels = pixels.sort_values(by = ["bin1_id","bin2_id"]).reset_index(drop = True)
            return pixels
        index,bins = _call_bins(chroms,chromsize_txt)
        if return_bins:
            return bins
        pixels = _call_pixels(contact_txt,index)
        if save_path is not None:
            if not os.path.exists(save_path) :
                os.makedirs(save_path) 
            cooler.create_cooler(save_path+'/'+prefix+'_'+str(res/1000)+'k.cool',bins,pixels,ordered=True,ensure_sorted=True)
        return bins,pixels

