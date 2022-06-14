import os,glob
import pandas as pd
import numpy as np
import cooler
import seaborn as sns
from tqdm import tqdm
from multiprocessing import Pool
from matplotlib import  pyplot as plt

from .loading import AutoLoad

class PreprocessCooler():
    def __init__(self,url):
        self.path = url
        self.clr = cooler.Cooler(url)
        self.bins = self.clr.bins()[:]
        self.px = self.clr.pixels()[:]
        self.resolution = self.clr.binsize
        
    def cal_intra_slow(self) : 
        # deprecated now cause time comsuming, using cal_intra instead.
        px = self.px
        intra = px[px.bin1_id.apply(lambda x : self.bins.chrom[x]) == px.bin2_id.apply(lambda x : self.bins.chrom[x])]["count"].sum()
        return intra
    
    def cal_intra(self) : 
        # calculate hic intra-chromosome contacts, which means conatcts within each chromsomes. inter = total - intra
        intra = [self.clr.matrix(balance= False,as_pixels = True).fetch(i)["count"].sum() for i in self.clr.chromnames]
        return np.sum(intra)
    
    def cal_amount(self) :
        # calculate hic total contacts
        return self.px["count"].sum()

    def cal_non_zero_bins(self) : 
        # calculate hic non zero bins
        return len(self.px)
    
    def cal_short_range(self):
        # calculate hic contacts less than 2 MB genomic distance, not include 2 MB
        px = self.px
        short_range = px[px.bin2_id - px.bin1_id < 2000000/self.resolution]['count'].sum()
        return short_range
    
    def cal_mitotic(self):
        # calculate hic contacts less than 12 MB and more than 2MB genomic distance
        px = self.px
        ind = (px.bin2_id - px.bin1_id >= 2000000/self.resolution) & (px.bin2_id - px.bin1_id < 12000000/self.resolution)
        mitotic = px[ind]['count'].sum()
        return mitotic
    # Add other indicator here 1/5

def preprocess(cool_url):
    pre = PreprocessCooler(cool_url)
    nonzero_bins = pre.cal_non_zero_bins()
    amount = pre.cal_amount()
    intra = pre.cal_intra()
    short_range = pre.cal_short_range()
    mitotic = pre.cal_mitotic()
    # Add other indicator here 2/5
    stat = [cool_url,amount,nonzero_bins,intra,short_range,mitotic]
    return stat


def quality_control(scool_path,meta_path,rawpath_col=0,sample_col=1,label_col=2,intra_cutoff=0.5,min_cutoff=10000,nonzero_cutoff=10000,save_pass=True,save_fig=True):
    """
    quality_control : Calculate single cell statistical indicator and set a cutoff threshold of each indicator.

    Parameters
    ----------
    scool_path : str
        scool url 
    meta_path : str
        meta information of txt format. Which need at least three columns: rawpath, sample name, labels. You can use scool_meta create it from SRA RunTable or manually.
    rawpath_col : str, optional
        raw path col name of cells store in meta files, by default col 0 
    sample_col : str, optional
        sample col index of cells in meta files, by default col 1
    label_col : str, optional
        cell type col index of cells main label used in next step, by default col 2 # TODO: support multiple labels
    intra_cutoff : float, optional
        reserve intra/total > cutoff, due to bad cells, by default 0.5
    min_cutoff : int, optional
        reserve total contacts > cutoff,due to incomplete cells, by default 10000
    nonzero_cutoff : int, optional
        reserve non zero bins > cutoff, due to very sparse cells, by default 10000
    save_pass ： Boolean, optional
        return cells pass quality control csv at the same path of scool, False only return a DataFrame including control results
    save_fig :  Boolean, optional
        save quality control figures at the same path of scool
    Returns
    -------
    DataFrame
        Save statistics DataFrame including quality control result
    str, path
        if save_pass return a csv path
    
    """    
    cell_list,bins,pixel_list =  AutoLoad(scool_path).load_scool()
    with Pool(processes = 20) as pool:
        result = list(tqdm(pool.imap(preprocess,cell_list), total= len(cell_list)))

    meta = pd.read_csv(meta_path)  # TODO: Check meta list whether correctly according scool file. 
    stats = pd.DataFrame(result,columns = ["cool_url","total_contacts","nonzero_bins","intra","short_range","mitotic"])
    stats["raw_path"] = meta.iloc[:,rawpath_col].values
    stats["sample"] = meta.iloc[:,sample_col].values
    stats["label"] = meta.iloc[:,label_col].values
    stats["qualified"] = (stats.intra/stats.total_contacts > intra_cutoff) & (stats.total_contacts > min_cutoff)\
    & (stats.nonzero_bins > nonzero_cutoff) & (-pd.isna(stats.label))
    print(f"{stats.qualified.sum()} passed quality control, with more than {min_cutoff} contacts, more than {nonzero_cutoff} nonzero bins , intra contacts percentage more than {intra_cutoff*100} %.")
    # Add other indicator here 3/5 and also cutoff in this method
    if save_fig:
        stats["short(%)"] = stats.short_range / stats.total_contacts
        stats["mitotic(%)"] = stats.mitotic / stats.total_contacts
        stats["intra(%)"] = stats.intra / stats.total_contacts
        stats_info = stats[["total_contacts","nonzero_bins","intra(%)","short(%)","mitotic(%)"]]
        qc_dict = {"total_contacts" : min_cutoff,"nonzero_bins" : nonzero_cutoff,"intra(%)" : intra_cutoff}
        # Add other indicator here 3/5 to plot figure
        cols = stats_info.columns
        fig,axes=plt.subplots(1,len(cols),figsize=(5*len(cols),4))
        for i,col in enumerate(cols):
            ax = sns.histplot(stats_info[col],bins=20,kde=True,ax = axes[i])
            if qc_dict.get(col) is not None:
                ax.axvline(qc_dict.get(col), color='red',linestyle='--')
        plt.suptitle(os.path.basename(scool_path.split("_")[0]))
        plt.gcf().subplots_adjust(bottom=0.15)
        fig_path = scool_path.replace(".scool","_qc.png")
        plt.savefig(fig_path)
        print("Saving quality control figures results to " + fig_path)            
    if save_pass:
        qual_meta = stats[stats.qualified]
        qual_path = scool_path.replace(".scool","_qc_meta.csv")
        qual_meta.to_csv(qual_path,index=0)
        print("Saving quality control cells results to " + qual_path)
        return qual_path    
    return stats


def plot_stats(qual_path):
    stats = pd.read_csv(qual_path)
    stats_info = stats[["total_contacts","nonzero_bins","short(%)","mitotic(%)","intra(%)","label"]]
    g = sns.PairGrid(stats_info,hue="label")
    g.map_offdiag(sns.scatterplot)  
    g.map_diag(sns.histplot)  
    g = g.add_legend()
    #plt.suptitle(os.path.basename(qual_path.split("_")[0]))
    plt.savefig(qual_path.replace("meta.csv","stats.png"))
    print("Saving quality controlled cells stats to " + qual_path)
    # Add other indicator here 4/5 to plot figure
    return g


def scool_meta(file_list,runtable_path,prefix = 'Test',label_colname = 'Cell_type',save_dir = "."):
    """
    scool_meta : Generate meta data for scool from NCBI RunTable informations. Which do not fit all condition, only for some dataset.

    Parameters
    ----------
    file_list : list
        A list or list like of all txt files path.
    runtable_path : str
        NCBI SRA RunTable url, which including information. 
    label_colname : str
        column name of the label col you need in SRA RunTable. Such as Cell_type, source_name ...
    Returns
    -------
    str
        Save a DataFrame have at least three columns: RawPath, Sample Name, information including cell type or others.
        return csv path
    """    
    RunTable = pd.read_csv(runtable_path).drop_duplicates("Sample Name")
    label = pd.DataFrame(columns = ["RawPath"],data = file_list)
    label["Sample Name"] = label["RawPath"].apply(lambda x : os.path.basename(x).split("_")[0])
    meta = pd.merge(label,RunTable)
    meta = meta[['RawPath',"Sample Name",label_colname]]
    meta_path = save_dir + "/" + prefix + "_raw_meta.csv"
    meta.to_csv(meta_path,index=0)
    print(f"Create raw meta data transform from {runtable_path} to match all txt files... \nStore raw meta to {meta_path}")
    return meta_path

def scool_meta_dipc2021(txt_path="/store/zlwang/Workspace/data/Dip-C/pairs/GSM*",runtable_path="/store/zlwang/Workspace/data/Dip-C/GSE146397_metadata.cells_contacts_100k.txt",prefix = 'DipC2021',save_dir = "."):
    file_list = glob.glob(txt_path)
    RunTable = pd.read_table(runtable_path,header=None)
    RunTable[1] = RunTable[1].replace('Unknown',None)
    RunTable_dict = RunTable.set_index(0).to_dict()[1]
    label = pd.DataFrame(columns = ["path"],data = file_list)
    label["SampleName"] = label.path.apply(lambda x : os.path.basename(x).split(".")[0])
    label['label'] = label["SampleName"].apply(lambda x: RunTable_dict.get(x[11:]))
    meta_path = save_dir + "/" + prefix + "_raw_meta.csv"
    label.to_csv(meta_path,index=0)
    print(f"Create raw meta data transform from {runtable_path} to match all txt files... \nStore raw meta to {meta_path}")

def scool_meta_nagano(scool_path = "/store/zlwang/Workspace/data/Nagano2017/nagano_10kb_cell_types.scool", table_path = "/store/zlwang/Workspace/data/Nagano2017/GSE94489_README.txt",prefix = 'Nagano2017',save_dir = "."):
    cell_list,bins,pixel_list = AutoLoad(scool_path).load_scool()
    rdm = pd.read_table(table_path)
    ind = pd.DataFrame(cell_list,columns = ['path'])
    ind["Five_prime_barcode"] = [i.split("/")[-1].split("_")[2] for i in cell_list]
    ind["Three_prime_barcode"] = [i.split("/")[-1].split("_")[3] for i in cell_list]
    ind["Nature_Reference_ID"] = [i.split("/")[-1].split("_")[0] + "_" + i.split("/")[-1].split("_")[1] for i in cell_list]
    meta = pd.merge(ind,rdm,how="left",on = ["Five_prime_barcode","Three_prime_barcode","Nature_Reference_ID"])
    meta["label"] = meta.Name.apply(lambda x : x.split("-")[0])
    meta_path = save_dir + "/" + prefix + "_raw_meta.csv"
    meta = meta[['path','Name','label']]    
    meta.to_csv(meta_path,index=0)
    return meta_path