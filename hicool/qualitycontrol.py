# -*- coding:utf-8 -*-
""" CLI """
import click
import pandas as pd

from hicool.process import plot_stats, quality_control #,scool_meta


@click.command()
@click.option('-s','--scool_path', default="./Test/Test_1000000.scool",type=str, help='Support Linux regular format, use quotation mark. ')
@click.option('-m','--meta_path', default="./data/test_raw_meta.csv", help='Meta information of txt format. Which need at least three columns: rawpath, sample name, labels. You can use scool_meta create it from SRA RunTable or manually.')

@click.option('-ic','--intra_cutoff', default=0.5, help='reserve intra/total(intra+inter) contacts > cutoff due to bad cells.')
@click.option('-mc','--min_cutoff', default=10000, help='reserve total contacts > cutoff due to incomplete cells.')
@click.option('-nc','--nonzero_cutoff', default=10000, help='reserve nonzero bins contacts > cutoff due to very sparse cells, pay attention to adjust it when you only analysis several chromsomes.')
# Add other cutoff indicator here 5/5

def qc_scool(scool_path,meta_path,intra_cutoff,min_cutoff,nonzero_cutoff):
    'Run quality control step with certain cutoff. Create scool files with cells passing quality control and plot stats figures. '
    qual_scool,qual_path = quality_control(scool_path,meta_path,intra_cutoff=intra_cutoff,min_cutoff=min_cutoff,nonzero_cutoff=nonzero_cutoff)
    g = plot_stats(qual_path)

if __name__ == '__main__':
    qc_scool()
