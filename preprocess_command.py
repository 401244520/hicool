# -*- coding:utf-8 -*-
""" CLI """
import click
import pandas as pd

from Module import plot_stats, quality_control, txt2scool


@click.command()
#@click.option('-i','--txt_path', default="./data/GSM*",type=str, help='Support Linux regular format, use quotation mark. ')
@click.option('-m','--meta_path', default="./data/DipC2021_raw_meta.csv", help='Meta information of txt format. Which need at least three columns: rawpath, sample name, labels. You can use scool_meta create it from SRA RunTable or manually.')
@click.option('-c','--chromsizes', default="./data/hg38.chrom.sizes", help='Path to chromsizes file download from NCBI.')
@click.option('-o','--save_dir', default="Test", help='Output dir path, all the files produced in this folder.')
@click.option('-p','--prefix', default="Test", help='Prefix of your Project.')

@click.option('-n','--nproc', default=20, help='Number of cpu process.')
@click.option('-r','--resolutions', default='1000000,100000,10000', type = str, help='A list of resolutions of scool files you want, spilt using comma.')
@click.option('-l','--chrom_list', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,\
chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chrX,chrY", type = str, help='A list of chroms you want analysis, spilt using comma.')

@click.option('-cc','--chroms_col', default='1,3', help='chromsome col of pairs txt files.')
@click.option('-pc','--positions_col', default='2,4', help='position col of pairs txt files.')

@click.option('-ic','--intra_cutoff', default=0.5, help='reserve intra/total(intra+inter) contacts > cutoff due to bad cells.')
@click.option('-mc','--min_cutoff', default=10000, help='reserve total contacts > cutoff due to incomplete cells.')
@click.option('-nc','--nonzero_cutoff', default=10000, help='reserve nonzero bins contacts > cutoff due to very sparse cells, pay attention to adjust it when you only analysis several chromsomes.')
# Add other cutoff indicator here 5/5

@click.option('-resume','--create_scool_resume', is_flag=True, help=' If you have run this workflow once sucessfully and only want to crate scool file in new resolution. Please set it as True')


def create_scool(meta_path,save_dir,prefix,nproc,resolutions,chromsizes,chroms_col,positions_col,chrom_list,intra_cutoff,min_cutoff,nonzero_cutoff,create_scool_resume):
    'Preprocess txt contacts files into scool and run quality control step with certain cutoff. Then create multi resolution scool files with cells passing quality control. '
    file_list = pd.read_csv(meta_path).iloc[:,0]
    # Default run qc in 1MB resolution
    chrom_list = chrom_list.split(',')
    print(chrom_list)
    chroms = chroms_col.split(',')
    positions = positions_col.split(',')
    if not create_scool_resume :
        scool_path = txt2scool(file_list,chrom_list,chromsizes,prefix+'_raw',nproc=nproc,chroms=chroms,positions=positions,resolution=1000000,save_dir=save_dir)
        # Add other indicator here 5/5
        qual_path = quality_control(scool_path,meta_path,intra_cutoff=intra_cutoff,min_cutoff=min_cutoff,nonzero_cutoff=nonzero_cutoff)
        g = plot_stats(qual_path)
        # Create multi-resolution scool files. 
        qc_cells = pd.read_csv(qual_path).raw_path.values
    else:
        qc_cells = file_list #TODO: Add exception Handling when using -resume
    resolutions = resolutions.split(",")
    print(resolutions)
    for res in resolutions:
        txt2scool(qc_cells,chrom_list,chromsizes,prefix,nproc=nproc,chroms=chroms,positions=positions,resolution=int(res),save_dir=save_dir)
    print(f'Creating resolution {resolutions} of quality control cells sucessfully.\n')

if __name__ == '__main__':
    create_scool()
