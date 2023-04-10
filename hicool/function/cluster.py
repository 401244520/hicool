import pandas as pd
import cooler
from cooler.reduce import CoolerMerger

from ..process import AutoLoad

from sklearn.metrics import silhouette_score,adjusted_rand_score, pairwise_distances



def cluster_cell(scool_path,label):
    cell_list = AutoLoad(scool_path).load_scool_cells()
    # check label length
    if len(cell_list) != len(label):
        raise Exception("Label length does not match scool files.")
    scool_cluster = scool_path.replace(".scool","_cluster.scool")
    bins = cooler.Cooler(scool_path).bins()[:]
    cell_name_pixels_dict = {}
    grps = pd.DataFrame([cell_list,label]).T.groupby(1)
    inds = grps.sum().index
    for ind in inds:
        gr = grps.get_group(ind)[0]
        cluster = [cooler.Cooler(cell) for cell in gr]
        cell_name_pixels_dict[ind] = CoolerMerger(cluster,maxbuf=1e9)
    cooler.create_scool(scool_cluster,bins,cell_name_pixels_dict)
    return scool_cluster 

def merge_cell(scool_path):
    cell_list = AutoLoad(scool_path).load_scool_cells()
    merged_cool = scool_path.replace(".scool","_merged.cool")
    bins = cooler.Cooler(scool_path).bins()[:]
    merged = CoolerMerger([cooler.Cooler(cell) for cell in cell_list],maxbuf=1e9)
    cooler.create_cooler(merged_cool,bins,merged)
    return  merged_cool



