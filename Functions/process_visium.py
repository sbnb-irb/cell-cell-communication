import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io
import os
import numpy as np

from matplotlib.image import imread
import json


def visium_noh5(main_path_exp, main_path_spatial, library_id = 'library_id'): 
# create adata object to save the GEx data
    adata = sc.read_10x_mtx(
        main_path_exp,  # the directory with the `.mtx` file
        var_names='gene_symbols',# use gene symbols for the variable names (variables-axis index)
        cache=True) 

    adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`   
    
## Prepare for spatial
    files = dict(
        tissue_positions_file=main_path_spatial + 'tissue_positions_list.csv',
        scalefactors_json_file=main_path_spatial + 'scalefactors_json.json',
        hires_image=main_path_spatial + 'tissue_hires_image.png',
        lowres_image=main_path_spatial + 'tissue_lowres_image.png',
    )
    adata.uns["spatial"] = dict()
    adata.uns["spatial"][library_id] = dict()
    adata.uns["spatial"][library_id]['images'] = dict()
    adata.uns["spatial"][library_id]['images']['hires'] = imread(files['hires_image'])
    adata.uns["spatial"][library_id]['images']['lowres'] = imread(files['lowres_image'])
    
    f = open(files['scalefactors_json_file'])
    adata.uns["spatial"][library_id]['scalefactors'] = json.load(f)
    
    positions = pd.read_csv(files['tissue_positions_file'], header=None)
    positions.columns = [
        'barcode',
        'in_tissue',
        'array_row',
        'array_col',
        'pxl_col_in_fullres',
        'pxl_row_in_fullres',
    ]
    positions.index = positions['barcode']
    adata.obs = adata.obs.join(positions, how="left")
    adata.obsm['spatial'] = adata.obs[
        ['pxl_row_in_fullres', 'pxl_col_in_fullres']
    ].to_numpy()
    adata.obs.drop(
        columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'],
        inplace=True,
    )
    return(adata)