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
![image](https://user-images.githubusercontent.com/47477490/230857501-c44798f4-0c8f-44bd-83c6-cd904eaed441.png)

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


# Compartment and TAD
comp = compartment_decomposition("../../data/DipC2019/DipC2019_100000_qc_merged.cool")
MOE,retina = AutoLoad("../../data/DipC2019/DipC2019_100000_qc_cluster.scool").load_scool_cells()
![image](https://user-images.githubusercontent.com/47477490/230858985-0914d732-5cfb-44dd-9080-a21fa8ac07e4.png)

tad = tad_insulation("../../data/DipC2019/DipC2019_100000_qc_merged.cool")
![image](https://user-images.githubusercontent.com/47477490/230860117-dc2632ef-89ef-4239-af2f-78d3977a5d34.png)

cell_list = AutoLoad("../../data/DipC2019/DipC2019_100000_qc.scool").load_scool_cells()
with Pool(processes = 36) as pool:
    comp_100k = list(tqdm(pool.imap(compartment_decomposition,cell_list), total= len(cell_list)))
sns.clustermap(pd.concat(comp_100k)["E1"],col_cluster=False,cmap="RdBu",figsize=(12,4),row_colors=colors,vmax = 1.5 ,vmin = -1.5)
![image](https://user-images.githubusercontent.com/47477490/230857985-c610a4c2-33cc-432a-8f7e-bf49e72c0386.png)
sns.clustermap(t_conserve,row_cluster=1,cmap="RdBu_r",figsize=(12,4),row_colors=colors,vmax = 1 ,vmin = 0)
![image](https://user-images.githubusercontent.com/47477490/230860951-66c836da-b419-45f8-a28b-afa2eea4a414.png)
n_clusters = 4
cluster = AgglomerativeClustering(n_clusters=n_clusters, affinity='euclidean', linkage='ward')  
clr = cluster.fit_predict(data_scaled)
for i in range(n_clusters):
    plt.subplot(100+n_clusters*10+i+1)
    sns.heatmap(data_scaled[clr == i].T,cmap="RdBu_r",vmax = 1)
    plt.title("Cluster "+ str(i))
![image](https://user-images.githubusercontent.com/47477490/230860952-7ab4fd9a-9353-4087-b65d-215150af1bcb.png)


cell_list = AutoLoad("../../data/DipC2019/DipC2019_50000_qc.scool").load_scool_cells()
with Pool(processes = 36) as pool:
    tad_50k = list(tqdm(pool.imap(tad_insulation,cell_list), total= len(cell_list)))
![image](https://user-images.githubusercontent.com/47477490/230858044-8e0b0e07-000f-4d1e-8252-eeae9641180f.png)



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
![image](https://user-images.githubusercontent.com/47477490/230857441-1c4f2680-07cf-4297-9c27-15b3b50fe24b.png)


# Save the HiCool object as a HiCool file
hc.save_as("path/to/output.hicool")

# Get metadata information
hc.info()

# Convert the HiCool object to a Scanpy AnnData object
sce = hc.to_scanpy(embedding_name="my_embedding")
hig = hc.to_higashi() # Pending
fast_hig = hc.to_fast_higashi()




