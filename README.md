# hicool
**HiCool is not complete now, just a preliminary development stage, for the time being, only support copying tools to local development, we are developing a complete version that can be installed and detailed documentation, so stay tuned!**

# Data Available
DipC2019 processed data and scripts [DipC2019](https://pan.baidu.com/s/1P1weJG0J1FdpYmGWQvfgCQ?pwd=dipc) \
More usage instance in scripts/(dataset).scool.ipynb \
All dataset at 1MB resolution and Embryo data [HiCool](https://pan.baidu.com/s/147r4oojByKrgVY61Unzr3Q?pwd=hico)
# hicool
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
![image](https://user-images.githubusercontent.com/47477490/230857501-c44798f4-0c8f-44bd-83c6-cd904eaed441.png)
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
![chr1_compartment_E1](https://user-images.githubusercontent.com/47477490/230863478-b08a8caf-45df-4e51-81a9-ab32067e5d4b.png)
```
comp = compartment_decomposition("../../data/DipC2019/DipC2019_100000_qc_merged.cool") 
MOE,retina = AutoLoad("../../data/DipC2019/DipC2019_100000_qc_cluster.scool").load_scool_cells() 
```
![chr1_tad](https://user-images.githubusercontent.com/47477490/230863583-806adc5c-c0db-47e1-bad0-9b0facf14d03.png)
```
tad = tad_insulation("../../data/DipC2019/DipC2019_100000_qc_merged.cool") 
```
![image](https://user-images.githubusercontent.com/47477490/230865140-926a9be7-e09d-41e8-8cb0-7133aad887d0.png)
```
cell_list = AutoLoad("../../data/DipC2019/DipC2019_100000_qc.scool").load_scool_cells() 
comp_100k = list(pool.imap(compartment_decomposition,cell_list)) 
sns.clustermap(pd.concat(comp_100k)["E1"],col_cluster=False,cmap="RdBu",figsize=(12,4),row_colors=colors,vmax = 1 ,vmin = -1) 
```
![image](https://user-images.githubusercontent.com/47477490/230860952-7ab4fd9a-9353-4087-b65d-215150af1bcb.png)
```
hc = HiCool("../../data/DipC2019/DipC2019_100000_qc.scool")
bd_chr1 = find_boundaries(comp_chr1)
c_cons,b = compress_matrix(hc,operation="comp_conservation",bin=bd_chr1,chrom="chr1")
```
![chr1_tad50K_cluster](https://user-images.githubusercontent.com/47477490/230863808-fc0a5833-b982-44aa-a17f-96438b2a17ba.png)
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
![image](https://user-images.githubusercontent.com/47477490/230856698-1b2f9060-efae-4874-addc-991e4dce9c84.png)


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
![image](https://user-images.githubusercontent.com/47477490/230856733-c809fa3e-02d1-464f-b0c2-da6dcbb51653.png)


# Cell Clustering 
![image](https://user-images.githubusercontent.com/47477490/230856825-78feb89b-f6fc-496b-87cb-2ae65b4a5bbb.png)
```
emb_names = list(hc.embedding.keys()) 
embs_pca,embs_tsne,embs_mds,embs_umap = PCA(2),TSNE(2),MDS(2),UMAP() 
for emb_name in emb_names:
    embs_pca.fit_transform(hc.embedding[emb_name]) 
```
![image](https://user-images.githubusercontent.com/47477490/230857441-1c4f2680-07cf-4297-9c27-15b3b50fe24b.png)
```
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
```
![image](https://user-images.githubusercontent.com/47477490/230866067-c0321bbd-bbf7-414c-990b-bd52e90b1f42.png)
```
hig = hc.to_higashi() # Pending 
fast_hig = hc.to_fast_higashi() 
```



