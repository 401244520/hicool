import glob

import pandas as pd

from Module import plot_stats, quality_control, scool_meta, txt2scool

chromsizes = "./data/mm10.chrom.sizes.txt"
chrom_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",\
            "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX"]
runtable_path = "./data/GSE121791RunTable.txt"

prefix = "Test"
txt_path = "./data/Test/GSM*"
file_list = glob.glob(txt_path)

save_dir = './Test'

scool_path = txt2scool(file_list,chrom_list,chromsizes,prefix,resolution=1000000,save_dir=save_dir)

stats_path = scool_meta(file_list,runtable_path,prefix,save_dir=save_dir)

qual_path = quality_control(scool_path,stats_path)

g = plot_stats(qual_path)

qc_cells = pd.read_csv(qual_path).raw_path.values

scool_100k = txt2scool(qc_cells,chrom_list,chromsizes,prefix,resolution=100000,save_dir=save_dir)
scool_10k = txt2scool(qc_cells,chrom_list,chromsizes,prefix,nproc=48,resolution=10000,save_dir=save_dir)
# scool_1k = txt2scool(qc_cells,chrom_list,chromsizes,prefix,nproc=72,resolution=1000,save_dir=save_dir)
