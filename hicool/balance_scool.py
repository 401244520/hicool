# -*- coding:utf-8 -*-
""" CLI """
import click
import cooler
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

from process import AutoLoad


@click.command()
@click.option('-s','--scool_path', help='scool path you want balance and normalize.')
@click.option('-m','--max_iters', default=200, help='max iter times, default is 200.')


def balance_scool(scool_path,max_iters):
    cell_list = AutoLoad(scool_path).load_scool_cells()
    for cell in tqdm(cell_list):
        cooler.balance_cooler(cooler.Cooler(cell),store=True,max_iters=max_iters)
        
if __name__ == '__main__':
    balance_scool()