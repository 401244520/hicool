# hicool
HiCool
HiCool is a Python package that provides an object-oriented interface for working with HiCool and SCool data files. This package allows you to load, manipulate, and save HiCool and SCool files with ease.

# Installation
pip install hicool
cp hicool ./ && conda install hicool/requirements.txt -r
# Usage
from hicool.tools import compress_matrix,HiCool

from hicool.process import AutoLoad
from hicool.process.statistics import quality_control
from hicool.function.similarity import cal_similarity
from hicool.function.features import tad_insulation,compartment_decomposition

# Quality control
![image](https://user-images.githubusercontent.com/47477490/230853452-69784db2-c75e-4ad8-a883-00e08efc8747.png)
scool_qc,meta_qc = quality_control(scool_path,meta_path) \
quality_control : Calculate single cell statistical indicator and set a cutoff threshold of each indicator.
Default:
    rawpath_col: int = 0,
    sample_col: int = 1,
    label_col: int = 2,
    intra_cutoff: float = 0.5,
    min_cutoff: int = 10000,
    nonzero_cutoff: int = 10000,
    nproc: int = 20,
    save_pass: bool = True,
    save_fig: bool = True


# Load a HiCool or Scool file
hc = HiCool(scool_qc,meta_qc)
or you can directly loading .hicool file.
hc = HiCool(hicool_path)


# Cell embedding
compress_matrix is a wrapped function that compresses a HiCool object into a 2D matrix. Each row in the matrix represents a cell, and each column represents a compressed feature.

matrix, features = compress_matrix(hc, chrom="chr1", operation="bin_strength",regions=chrom_regions)
operations: \
"bin_degree" :  partial(node_degree,chrom=chrom,balance=balance,field=field), \
"bin_strength" :  partial(node_strength,chrom=chrom,balance=balance,field=field),\
"region_strength" : partial(strength_region,regions=regions,balance=balance,field=field,**kwargs),\
"region_conservation" : partial(conservation_region,regions=regions,balance=balance,field=field,**kwargs),\
"comp_strength" : partial(_node_aggregation_strength,**kwargs),\
"comp_conservation" : partial(_node_conservation,**kwargs),\
"strata_strength" : partial(strata_strength,chrom=chrom,balance=balance,field=field,**kwargs),\
![image](https://user-images.githubusercontent.com/47477490/230856698-1b2f9060-efae-4874-addc-991e4dce9c84.png)


# Cell network
from hicool.function.similarity import cal_similarity
cell_list = AutoLoad(hc_path).load_hicool_cells()\
hc.network["hicrep_chr1"] = cal_similarity(cell_list,method="hicrep",chrom="chr1")\

Or you can calculate correlation of cell embedding \
similarity = np.corrcoef(hc.embedding["bin_degree"]) \
distance_mat = np.sqrt(2 - 2 * similarity) \
hc.network["bin_degree"] = distance_mat \
![image](https://user-images.githubusercontent.com/47477490/230856733-c809fa3e-02d1-464f-b0c2-da6dcbb51653.png)


# Cell Clustering 

emb_names = list(hc.embedding.keys()) \
embs_pca,embs_tsne,embs_mds,embs_umap = [],[],[],[] \
for emb_name in emb_names: \
    embs_pca.append(PCA(2).fit_transform(hc.embedding[emb_name])) \
    embs_tsne.append(TSNE(2).fit_transform(hc.embedding[emb_name])) \
    embs_mds.append(MDS(2).fit_transform(hc.embedding[emb_name])) \
    embs_umap.append(UMAP(2).fit_transform(hc.embedding[emb_name])) \

![image](https://user-images.githubusercontent.com/47477490/230856825-78feb89b-f6fc-496b-87cb-2ae65b4a5bbb.png)

from hicool.function.estimation import cal_acroc,plt

acroc = cal_acroc(hc.embedding["bin_degree"],label)
plt.title("bin_degree")
![image](https://user-images.githubusercontent.com/47477490/230857334-a94dd90b-9e00-41c7-a4af-0b6b7c614e89.png)

# Save the HiCool object as a HiCool file
hc.save_as("path/to/output.hicool")

# Get metadata information
hc.info()

# Convert the HiCool object to a Scanpy AnnData object
sce = hc.to_scanpy(embedding_name="my_embedding")
hig = hc.to_higashi() # Pending
fast_hig = hc.to_fast_higashi()




