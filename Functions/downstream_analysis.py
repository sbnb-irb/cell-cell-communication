from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import liana as li
import scanpy as sc
import pandas as pd
import numpy as np
import random
from collections import Counter
from scipy.spatial import distance
import ktplotspy as kpy
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib
import seaborn as sns

title_size = 12
labelout_size = 10
label_size = 8
ticks_size = 6


def obtain_cell_counts(path):
    """
    This function reads the cyto_cell2spot file and returns the number of cells for each cell type (breast cancer)
    Input:
        - path: path to the assigned_locations file
    """
    cyto_cell2spot = pd.read_csv(path)
    cyto_cell2spot.drop_duplicates('OriginalCID', inplace=True)
    total_cancer_cells = sum(value for key, value in Counter(cyto_cell2spot.CellType).items() if "Cancer" in key)
    cell_counts = {key: value for key, value in Counter(cyto_cell2spot.CellType).items() if "Cancer" not in key}
    cell_counts["Cancer Cells"] = total_cancer_cells
    return cell_counts

def counts_major_cell_by_region(spatialfile, adatafile, cytofile, type2color, cols, region = True):
    '''
    This function returns the number of cells for each cell type in each region (heart)
    Input:
        - spatialfile: path to the spatial adata file
        - adatafile: path to the sc adata file
        - cytofile: path to the cyto_cell2spot file
        - type2color: dictionary with the colors for the different regions
        - cols: columns in the adata.obs that contains the cell type information (the first one has to be concordant with cytospace annotation)
        - region: boolean to return the counts per celltype without dividing by region
    Output:
        - plot_df: dataframe with the number of cells for each cell type in each region
        - p_c2l: AnnData object with the spatial division
    '''
    ### Load data
    p_c2l = sc.read(spatialfile)
    adata = sc.read(adatafile)
    cyto_cell2spot = pd.read_csv(cytofile)
    ### Obtain colors for the different regions
    p_c2l.obs['color'] = [type2color[i] for i in p_c2l.obs.annotation_final]
    ### Convert cell states to cell type
    info_cells = adata.obs[cols].drop_duplicates()
    state2cell = dict(zip(info_cells[cols[0]], info_cells[cols[1]]))
    cyto_cell2spot = cyto_cell2spot[cyto_cell2spot.OriginalCID.notna()]
    cyto_cell2spot['selected'] = [state2cell[i] for i in cyto_cell2spot.CellType]
    
    if region == False: 
        cell_counts = {key: value for key, value in Counter(cyto_cell2spot.selected).items()}
        return cell_counts

    spot2region = dict(zip(p_c2l.obs.index, p_c2l.obs.annotation_final))
    cyto_cell2spot['Region'] = [spot2region[i] for i in cyto_cell2spot.SpotID]

    ### Count
    data_for_plot = []
    for i in cyto_cell2spot.Region.unique(): 
        c = cyto_cell2spot[cyto_cell2spot.Region == i]
        cell_counts = Counter(c['selected'])
        for cell_type, count in cell_counts.items():
            data_for_plot.append({'Region': i, 'CellType': cell_type, 'Count': count})

    plot_df = pd.DataFrame(data_for_plot)
    return plot_df, p_c2l

def obtain_adata_embeddings(adata_file, list_loc, celltype_key = 'celltype_loc'): 
    """
    This function reads the adata file with the assigned cells to locations and returns the AnnData object with the embeddings
    Input:
        - adata_file: path to the adata file
        - list_loc: list of labels of the regions to include in the analysis
        - celltype_key: column in the adata.obs that contains the cell type information
    """
    adata = sc.read(adata_file)
    adata.uns['log1p']['base'] = None
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
    sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
    sc.tl.tsne(adata, use_rep="X_pca")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    adata = adata[adata.obs[celltype_key].isin(list_loc)]
    return adata


def pseudoreplicates (adata, condition, n = 3, layer = 'counts'): 
    """
    This function generates pseudoreplicates for the pseudobulk DEA
    Input: 
        - adata: AnnData object with the cells to compare
        - condition: column in the adata.obs that defines the condition
        - n: number of pseudoreplicates
    """
    pbs = []
    for cluster in adata.obs[condition].unique(): 
        cluster_adata = adata[adata.obs[condition] == cluster]
        if layer:
            cluster_adata.X = cluster_adata.layers[layer] ## use raw data

        indices = list(cluster_adata.obs_names)
        random.shuffle(indices)
        indices = np.array_split(np.array(indices), n)

        for i, pseudo_rep in enumerate(indices):
            cluster_adata_chunk = cluster_adata[cluster_adata.obs.index.isin(indices[i])]
            rep_adata = sc.AnnData(X = cluster_adata_chunk.X.sum(axis = 0), var = cluster_adata_chunk.var[[]])
            rep_adata.obs_names = [cluster_adata_chunk.obs[condition].iloc[0] + '_' + str(i)]
            rep_adata.obs['condition'] = cluster_adata_chunk.obs[condition].iloc[0]
            rep_adata.obs['replicate'] = i
            pbs.append(rep_adata)
    pb = sc.concat(pbs)
    return(pb)

def obtain_deg(pb, comparison, reference, design_factors = 'condition'): 
    """
    This function perform DEA using pseudobulk. You select the query and the reference condition
    Input:
        - pb: concat adata object with the different pseudobulk samples
        - comparison: the condition you want to compare with the reference
        - reference: the condition you want to use as reference
        - design_factors: the column in the metadata that contains the condition
    
    """
    counts = pd.DataFrame(pb.X, columns=pb.var_names)   
    dds = DeseqDataSet(
        counts = counts,
        metadata = pb.obs, 
        design_factors = design_factors
    )
    sc.pp.filter_genes(dds, min_cells = 1)
    dds.deseq2()
    stat_res = DeseqStats(dds, n_cpus=4, contrast = [design_factors, comparison, reference])
    stat_res.summary()
    res = stat_res.results_df
    return(res)

def select_regional_spots (cytospace_folder, adatafile, list_loc): 
    '''
    This function relate the regions in the annotated anndata object (with the space division) to the spots in cytospace. 
    It is used to obtain the input for the colocalization function. 
    Input: 
        - cytospace_folder: Folder where the cytoSPACE output is located
        - adatafile: .h5ad object with the spatial division in the target cell 
        - list_loc: list with the names of the regions to be related to the spots.
    Output: 
        - spots_loc1: list with the spots in the first region
        - spots_loc2: list with the spots in the second region
    '''
    ### output CytoSPACE
    cyto_cell2spot = pd.read_csv(f'{cytospace_folder}/assigned_locations.csv')
    cyto_spot2cell = pd.read_csv(f'{cytospace_folder}/cell_type_assignments_by_spot.csv')
    del(cyto_spot2cell['Total cells']) 

    ### select only cancer cells and assign them to the location (assign spots to locations)
    adata = sc.read(adatafile)
    cells_loc1 = adata.obs[adata.obs.celltype_loc == list_loc[0]].index
    cells_loc2 = adata.obs[adata.obs.celltype_loc == list_loc[1]].index

    spots_loc1 = cyto_cell2spot[cyto_cell2spot.OriginalCID.isin(cells_loc1)].SpotID.unique()
    spots_loc2 = cyto_cell2spot[cyto_cell2spot.OriginalCID.isin(cells_loc2)].SpotID.unique()

    return spots_loc1, spots_loc2

def compute_colocalization(folder_or_df, cell, loc1_spots, loc2_spots):
    '''
    Function to compute the binary matrix representing the spot profile for each cell type and 
    then compute the colocalization between the two groups of query cells and the rest of the cells in the dataset
    Input:
        - folder_or_df: folder where the cytoSPACE output is located
        - cell: cell type to analyze
        - loc1_spots: list of spots where the first group of query cells is located
        - loc2_spots: list of spots where the second group of query cells is located
    '''
    # Load the cytoSPACE output
    if isinstance(folder_or_df, str):
        cyto_spot2cell = pd.read_csv(f'{folder_or_df}/cell_type_assignments_by_spot.csv')
        del(cyto_spot2cell['Total cells']) ## remove total cells column
        # Remove the query cell column (we want to separate between the two groups of query cells previously defined)
        spot2cell_qcell = cyto_spot2cell[cyto_spot2cell.columns[~cyto_spot2cell.columns.str.contains(cell)]]
        spot2cell_qcell.set_index('SpotID', inplace = True)
    else: 
        spot2cell_qcell = folder_or_df[folder_or_df.columns[~folder_or_df.columns.str.contains(cell)]]

    # Binarize the matrix
    spot2cell_qcell = spot2cell_qcell.applymap(lambda x: 1 if x > 0 else 0)

    # Compute the spot profile for this two subsets of query cells
    loc2_bin, loc1_bin = [],[]
    for i in spot2cell_qcell.index: 
        if i in loc1_spots:
            loc1_bin.append(1)
        else:
            loc1_bin.append(0)        
        if i in loc2_spots:
            loc2_bin.append(1)
        else:
            loc2_bin.append(0)

    # Compute the colocalization (Jaccard similarity)
    columns = spot2cell_qcell.columns

    similarities_loc1 = {}

    for i in range(len(columns)): 
        col1 = spot2cell_qcell[columns[i]]
        col2 = loc1_bin
        jaccard_similarity = 1 - distance.jaccard(col1, col2)
        similarities_loc1[columns[i]] = jaccard_similarity           
        similarities_loc2 = {}

    for i in range(len(columns)): 
        col1 = spot2cell_qcell[columns[i]]
        col2 = loc2_bin
        jaccard_similarity = 1 - distance.jaccard(col1, col2)
        similarities_loc2[columns[i]] = jaccard_similarity

    if similarities_loc1.keys() == similarities_loc2.keys(): 
        closest_cells_df = pd.DataFrame([list(similarities_loc1.keys()), list(similarities_loc1.values()),list(similarities_loc2.values())]).T
        closest_cells_df.columns = ['cells', 'JS_loc1', 'JS_loc2']
        return closest_cells_df
    
    else:
        print('Error: the keys in the two dictionaries are different')
        return similarities_loc1, similarities_loc2
    
def create_filtered_adata_CP_subsampling (adatafile, outfolder, list_loc, celltype_key = 'celltype_loc'):
    '''
    Function to create a filtered adata and metadata with the same number of cells for each target cell
    Input:
        - adatafile: path to the annotated adata file
        - metafile: path to the metadata file
        - outfolder: path to save the output
        - list_loc: list of labels of the regions to include in the analysis
        - cell_target: the cell type to remove from the analysis
        - celltype_key: column in the adata.obs that contains the cell type information
    '''
    # Load adata and metadata
    adata = sc.read(adatafile)
    # Select the non-target cells 
    index_no_cancer = adata.obs[~adata.obs[celltype_key].isin(list_loc)].index
    # Select the minimum number of cells in the target cells
    n_ccells = adata.obs[adata.obs[celltype_key].isin(list_loc)][celltype_key].astype(str).value_counts().min()
    # Select the same number of cells for each target cell
    ind_selected = []
    for i in list_loc: 
        ind_total = adata.obs[adata.obs[celltype_key] == i].index
        ind_selected.append(np.random.choice(ind_total, n_ccells, replace=False))
    selected_ind_CP = np.concatenate([np.concatenate(ind_selected), index_no_cancer])
    # Create the new adata and metadata
    adata_cp = adata[adata.obs_names.isin(selected_ind_CP)]
    ## Create the meta dataframe for the cellphonedb
    meta_cp = pd.DataFrame([adata_cp.obs.index,adata_cp.obs[celltype_key]]).T
    meta_cp.columns = ['Cell', 'cell_type']
    adata_cp.write(outfolder + 'subsampling.h5ad')
    meta_cp.to_csv(outfolder + 'metadata_subsampling.tsv', sep = '\t')

    return adata_cp
    

def run_cellphonedb(cpdb_file_path, meta_file_path, counts_file_path, out_path, celltype_key, groups, pvalfile = ''): 
    '''
    Function to run cellphonedb statistical analysis and obtain the number of significant interactions for a given cell types
    Input:
        - cpdb_file_path: path to the cellphonedb database file (cellphonedb.zip)
        - meta_file_path: path to the metadata file 
        - counts_file_path: path to the counts file (.h5ad)
        - out_path: path to save the output
        - celltype_key: column in the adata.obs that contains the cell type information
        - groups: list of cell types to compare
        - pvalfile: path to the pval file if you already compute it(optional)
    Output: 
        - counts_qcell: dataframe with the number of significant interactions with all cells for the cell types in the groups list
    '''
    if pvalfile: 
        pvalues = pd.read_csv(pvalfile, sep = '\t', index_col = 0)
    else: 
        deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
            cpdb_file_path = cpdb_file_path,
            meta_file_path = meta_file_path,
            counts_file_path = counts_file_path,
            counts_data = 'gene_name',
            subsampling_num_cells = None,
            output_path = out_path)
    
    adata = sc.read_h5ad(counts_file_path)
    counts_nosim = kpy.plot_cpdb_heatmap(
        adata = adata,
        pvals = pvalues,
        celltype_key = celltype_key,
        figsize = (10,10),
        title = "Number of significant interactions",
        symmetrical = True,
        return_tables = True
    )
    counts_qcell = counts_nosim['count_network'][groups]
    return(counts_qcell)


def plot_CCCvsColoc (closest_cells_df, counts_ccc, subset, path_save = ''): 
    '''
    Function to show the colocalization results vs the cell-cell communication. 
    Input: 
        - closest_cells_df = output of the colocalization function. df with cells, colocalization w/ loc1 and loc2
        - counts_ccc = oytput of the run_CellPhone or run_LIANA functions. df with cells, interactions counts w/ loc1 and loc2 
        - subset of cells to repesent
        - path_save: path to save the plot
    Output: 
        - colocalization vs ccc plot
        - data_cc_coloc: df with coloc and ccc columns
    '''
        
    selected_CC = counts_ccc[counts_ccc.cell.isin(subset)]
    selected_CC.set_index('cell', inplace = True)

    selected_coloc = closest_cells_df[closest_cells_df.cells.isin(subset)]
    selected_coloc.sort_values('cells', inplace = True)
    selected_coloc.set_index('cells', inplace = True)

    data_cc_coloc = selected_CC.merge(selected_coloc, left_index = True, right_index = True)
    data_cc_coloc.reset_index(inplace = True)

    data_cc_coloc.columns = ['cells', 'Cancer_loc1', 'Cancer_loc2', 'JS_loc1', 'JS_loc2']

    data = data_cc_coloc

    # Create the plot
    fig, ax = plt.subplots(figsize=(5,6), dpi=200)


    # Plot the jaccard similarities as horizontal lolliplots
    for i, row in data.iterrows():

        ax.plot([0, row['JS_loc2']], [i, i], color='#F39426', lw=2, ls = '--', zorder=10)
        ax.scatter([row['JS_loc2']], [i], color='#F39426', s=40, marker='*', zorder=10)

        ax.plot([0, row['JS_loc1']], [i + 0.3, i + 0.3], color='#B16BA8', lw=2, ls = '--', zorder=10)
        ax.scatter([row['JS_loc1']], [i + 0.3], color='#B16BA8', s=40,  marker='*', zorder=10)
        
    # Create a second x-axis to represent the interactions with cancer cells
    ax2 = ax.twiny() 

    # Plot the number of significant LR interctions as horizontal lolliplots
    for i, row in data.iterrows():
        ax2.plot([0, row['Cancer_loc2']], [i + 0.15, i + 0.15], color='#F39426', lw=2, zorder=10)
        ax2.scatter([row['Cancer_loc2']], [i + 0.15], color='#F39426', s=20, zorder=10)

        ax2.plot([0, row['Cancer_loc1']], [i + 0.45, i + 0.45], color='#B16BA8', lw=2, zorder=10)
        ax2.scatter([row['Cancer_loc1']], [i + 0.45], color='#B16BA8', s=20, zorder=10)
        

    ax.set_yticks(np.array((range(len(data['cells'])))) + 0.2, data['cells'], fontsize = labelout_size)


    ax.set_xlabel('Jaccard Similarity - Colocalization',  fontsize=labelout_size)
    ax2.set_xlabel('# Interactions ',  fontsize=labelout_size)

    ax.tick_params(axis='both', labelsize=labelout_size)



    legend_elements = [
        Line2D([0], [0], color='#B16BA8', lw=2, label='Location 1'),
        Line2D([0], [0], color='#F39426', lw=2, label='Location 2'),
        Line2D([0], [0], color='black', lw=2, ls='--', label='Jaccard Similarity'),
        Line2D([0], [0], color='black', lw=2, label='# Interactions')
    ]
    ax.legend(handles=legend_elements, fontsize=label_size, frameon = False, loc = 'lower right')
    ax.set_title('Colocalization vs number of LR interactions', fontsize=title_size, y = 1.1)
    if path_save: 
        plt.savefig(path_save+'.svg',format='svg', dpi=300, bbox_inches='tight')
        plt.savefig(path_save+'.png',format='png', dpi=300, facecolor='white', edgecolor='none', bbox_inches='tight')

    return data_cc_coloc

    


def boxplot_subsampling (df_plot, ax, patient, list_loc, filter = [], new_labels_l = ['Location 1', 'Location 2']):
    '''
    Function to create a boxplot with the number of interactions for the different cell types (we show the results from the different subsampling iterations)
    Input:
        - df_plot: dataframe with the number of interactions (all subsampling iterations)
        - ax: axis to plot the boxplot
        - patient: patient identifier
        - list_loc: list of labels of the regions to include in the analysis (they are used as palette)
        - filter: list of cell types to include in the analysis (if empty, all cell types are included)
    '''
    if filter:
        df_plot = df_plot[df_plot['index'].isin(filter)]
        
    sns.boxplot(x='index', y='inter', hue='celltype', data=df_plot, ax=ax, palette =list_loc, saturation=1, linewidth=0.8)
    sns.stripplot(x='index', y='inter', hue='celltype', data=df_plot, ax=ax, palette = list_loc, dodge=True, size=2, jitter=True, linewidth=0.1, edgecolor='black')

    box_patches = [patch for patch in ax.patches if type(patch) == matplotlib.patches.PathPatch]
    if len(box_patches) == 0:  # in matplotlib older than 3.5, the boxes are stored in ax2.artists
        box_patches = ax.artists
    num_patches = len(box_patches)
    lines_per_boxplot = len(ax.lines) // num_patches
    for i, patch in enumerate(box_patches):
        # Set the linecolor on the patch to the facecolor, and set the facecolor to None
        col = patch.get_facecolor()
        patch.set_edgecolor(col)
        patch.set_facecolor(col[:3] + (0.2,))
        

        # Each box has associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same color as above
        for line in ax.lines[i * lines_per_boxplot: (i + 1) * lines_per_boxplot]:
            line.set_color(col)
            line.set_mfc(col)  # facecolor of fliers
            line.set_mec(col)  # edgecolor of fliers


    # ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=90,fontsize = labelout_size-2)  # Rotate x-axis labels for better readability, if needed
    # Split long x-tick labels into multiple lines
    xtick_labels = ax.get_xticklabels()
    new_labels = []
    max_label_length = 30
    for label in xtick_labels:
        text = label.get_text()
        if len(text) > max_label_length:
            text = '\n'.join([text[i:i+max_label_length] for i in range(0, len(text), max_label_length)])
        new_labels.append(text)
    ax.set_xticklabels(new_labels, rotation=90, fontsize=labelout_size-2)


    ax.set_ylabel('# LR iteractions')
    ax.set_xlabel('')
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(handles, new_labels_l)

    for handle in legend.legend_handles:
        col = handle.get_facecolor()
        handle.set_edgecolor(col)
        handle.set_facecolor(col[:3] + (0.2,))
        handle.set_linewidth(1) 

    # Set the grid with a specific zorder
    ax.set_axisbelow(True)
    ax.grid(True, zorder=-1, ls = '--')

    ax.set_title( patient, fontsize=title_size)


def run_LIANA(countsfile, list_cells_loc, celltype_key): 
    '''
    Function to run LIANA and obtain the number of significant interactions for the defined spatial locations
    Input:
        - countsfile: path to the adata with the annotated regions (.h5ad)
        - list_cells_loc: list of regions to compare
        - celltype_key: column in the adata.obs that contains the cell type information
    Output:
        - counts_target: dataframe with the number of significant interactions with all cells for the cell types in the groups list
    '''
    lr_pairs = li.resource.select_resource('consensus')
    adata = sc.read(countsfile)
    adata.obs[celltype_key] = adata.obs[celltype_key].astype('category')
    li.mt.rank_aggregate(adata,
                    groupby=celltype_key,
                    resource_name = 'consensus',
                    expr_prop=0.1, # must be expressed in expr_prop fraction of cells
                    min_cells = 5,
                    n_perms = 100,
                    use_raw = False, # run on log- and library-normalized counts
                    verbose = True,
                    inplace = True)
    # Select the significant interactions
    LIANA_consensus_sig = adata.uns['liana_res'][adata.uns['liana_res']['magnitude_rank'] <= 0.05]
    LIANA_consensus_sig = LIANA_consensus_sig[['source', 'target', 'ligand_complex', 'receptor_complex','specificity_rank', 'magnitude_rank']]
    LIANA_consensus_sig_both = LIANA_consensus_sig[LIANA_consensus_sig['specificity_rank'] <= 0.05]
    # Create the LR pairs column
    LR_pairs = [('|').join([LIANA_consensus_sig_both.iloc[i]['ligand_complex'], LIANA_consensus_sig_both.iloc[i]['receptor_complex']]) for i in range(LIANA_consensus_sig_both.shape[0])]
    LIANA_consensus_sig_both['pairs'] = LR_pairs
    # Select the interactions for the two locations
    SC1_interactions = LIANA_consensus_sig_both.loc[(LIANA_consensus_sig_both.source == list_cells_loc[0]) | (LIANA_consensus_sig_both.target == list_cells_loc[0])]
    SC2_interactions = LIANA_consensus_sig_both.loc[(LIANA_consensus_sig_both.source == list_cells_loc[1]) | (LIANA_consensus_sig_both.target == list_cells_loc[1])]
    # Count the number of interactions for each cell type and create a df
    SC1_cell_counts = pd.DataFrame(Counter(np.concatenate([SC1_interactions['source'].tolist(), SC1_interactions['target'].tolist()])).most_common(), columns = ['cell', 'counts']).iloc[1:]
    SC2_cell_counts = pd.DataFrame(Counter(np.concatenate([SC2_interactions['source'].tolist(), SC2_interactions['target'].tolist()])).most_common(), columns = ['cell', 'counts']).iloc[1:]
    SC1_cell_counts = SC1_cell_counts.sort_values(by='cell')
    SC2_cell_counts = SC2_cell_counts.sort_values(by='cell')
    counts_target = SC1_cell_counts.merge(SC2_cell_counts, left_on = 'cell', right_on = 'cell', how='outer')
    return counts_target