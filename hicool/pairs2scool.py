# -*- coding:utf-8 -*-
""" CLI """
import click,glob
import pandas as pd



from process import txt2scool


@click.command()
@click.option('-i','--txt_path', default=None,type=str, help='Support Linux regular format, use quotation mark. ')
@click.option('-m','--meta_path', default="./data/test_raw_meta.csv", help='Meta information of txt format. Which need at least three columns: rawpath, sample name, labels. You can use scool_meta create it from SRA RunTable or manually.')
@click.option('-c','--chromsizes', default="./data/mm10.chrom.sizes", help='Path to chromsizes file download from NCBI.')
@click.option('-o','--save_dir', default="Test", help='Output dir path, all the files produced in this folder.')
@click.option('-p','--prefix', default="Test", help='Prefix of your Project.')

@click.option('-n','--nproc', default=20, help='Number of cpu process.')
@click.option('-r','--resolutions', default='1000000,100000,10000', type = str, help='A list of resolutions of scool files you want, spilt using comma.')
@click.option('-l','--chrom_list', default=None, type = str, help='A list of chroms you want analysis, spilt using comma.')
#"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX,chrY"

@click.option('-cc','--chroms_col', default='1,3', help='chromsome col of pairs txt files.')
@click.option('-pc','--positions_col', default='2,4', help='position col of pairs txt files.')


def create_scool(txt_path,meta_path,save_dir,prefix,nproc,resolutions,chromsizes,chroms_col,positions_col,chrom_list):
    'Preprocess valid pairs into scool with certain resolution. Then create multi resolution scool files. '
    if txt_path is not None :
        file_list = glob.glob(txt_path)
    else :
        file_list = pd.read_csv(meta_path).iloc[:,0]
    # Default run qc in 1MB resolution
    if chrom_list is not None :
        chrom_list = chrom_list.split(',')
    else :
        chrom_list = pd.read_csv(chromsizes,sep="\t",header=None)[0].tolist()
    print(chrom_list)
    chroms = chroms_col.split(',')
    positions = positions_col.split(',')
    resolutions = resolutions.split(",")
    print(resolutions)
    for res in resolutions:
        txt2scool(file_list,chrom_list,chromsizes,prefix,nproc=nproc,chroms=chroms,positions=positions,resolution=int(res),save_dir=save_dir)
    print(f'Create resolution {resolutions} of {len(file_list)} cells on {len(chrom_list)} chromosomes successful.\n')

if __name__ == '__main__':
    create_scool()