# hicool
**HiCool is not complete now, just a preliminary development stage, for the time being, only support copying tools to local development, we are developing a complete version that can be installed and detailed documentation, so stay tuned!**

HiCool is a Python package that provides an object-oriented interface for working with HiCool and Scool data files, which allows to load, manipulate, and save scHiC data and analysis easily.

# Installation
```
pip install hicool  # Not implement now
cp hicool ./ && conda install hicool/requirements.txt -r # Please copy or link the hicool floder to current path
```
# Usage
```
from hicool.tools import compress_matrix,HiCool 
from hicool.process import AutoLoad 
from hicool.process.statistics import quality_control 
from hicool.function.similarity import cal_similarity 
from hicool.function.features import tad_insulation,compartment_decomposition 
from hicool.function.estimation import cal_acroc 
```
# Quality control
```
scool_qc,meta_qc = quality_control(scool_path,meta_path) 
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
```

# Load a HiCool or Scool file
We provide generate HiCool instance from scool and directly loading hicool file.\
```
hc = HiCool(scool_qc,meta_qc)
hc = HiCool(hicool_path)
```

# Compartment and TAD
```
comp = compartment_decomposition("../../data/DipC2019/DipC2019_100000_qc_merged.cool") 
MOE,retina = AutoLoad("../../data/DipC2019/DipC2019_100000_qc_cluster.scool").load_scool_cells() 
```
```
tad = tad_insulation("../../data/DipC2019/DipC2019_100000_qc_merged.cool") 
```
```
cell_list = AutoLoad("../../data/DipC2019/DipC2019_100000_qc.scool").load_scool_cells() 
comp_100k = list(pool.imap(compartment_decomposition,cell_list)) 
sns.clustermap(pd.concat(comp_100k)["E1"],col_cluster=False,cmap="RdBu",figsize=(12,4),row_colors=colors,vmax = 1 ,vmin = -1) 
```
```
hc = HiCool("../../data/DipC2019/DipC2019_100000_qc.scool")
bd_chr1 = find_boundaries(comp_chr1)
c_cons,b = compress_matrix(hc,operation="comp_conservation",bin=bd_chr1,chrom="chr1")
```
```
cell_list = AutoLoad("../../data/DipC2019/DipC2019_50000_qc.scool").load_scool_cells() 
tad_50k = list(pool.imap(tad_insulation,cell_list))
```
# Cell embedding
compress_matrix is a wrapped function that compresses a HiCool object into a 2D matrix. Each row in the matrix represents a cell, and each column represents a compressed feature.
```
matrix, features = compress_matrix(hc, chrom="chr1", operation="bin_strength",regions=chrom_regions)
operations: 
"bin_degree" :  partial(node_degree,chrom=chrom,balance=balance,field=field), 
"bin_strength" :  partial(node_strength,chrom=chrom,balance=balance,field=field), 
"region_strength" : partial(strength_region,regions=regions,balance=balance,field=field,**kwargs), 
"region_conservation" : partial(conservation_region,regions=regions,balance=balance,field=field,**kwargs), 
"comp_strength" : partial(_node_aggregation_strength,**kwargs), 
"comp_conservation" : partial(_node_conservation,**kwargs), 
"strata_strength" : partial(strata_strength,chrom=chrom,balance=balance,field=field,**kwargs), 
```


# Cell network
```
from hicool.function.similarity import cal_similarity
cell_list = AutoLoad(hc_path).load_hicool_cells() 
hc.network["hicrep_chr1"] = cal_similarity(cell_list,method="hicrep",chrom="chr1") 
```
Or you can calculate correlation of cell embedding  
```
similarity = np.corrcoef(hc.embedding["bin_degree"])  
distance_mat = np.sqrt(2 - 2 * similarity) 
hc.network["bin_degree"] = distance_mat 
```


# Cell Clustering 
```
emb_names = list(hc.embedding.keys()) 
embs_pca,embs_tsne,embs_mds,embs_umap = PCA(2),TSNE(2),MDS(2),UMAP() 
for emb_name in emb_names:
    embs_pca.fit_transform(hc.embedding[emb_name]) 
acroc = cal_acroc(hc.embedding["bin_degree"],label) 
```


# Save the HiCool object as a HiCool file
```
hc.save_as("path/to/output.hicool")
```
# Get metadata information
```
hc.info()
```
# Convert the HiCool object to a Scanpy AnnData object
```
sce = hc.to_scanpy(embedding_name="my_embedding") 
hig = hc.to_higashi() # Pending 
fast_hig = hc.to_fast_higashi() 
```



