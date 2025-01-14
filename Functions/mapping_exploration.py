
global_path = './cell-cell-communication/' ## Path to the github downloaded repositoryimport pandas as pd
import numpy as np
import scanpy as sc
import pickle
import sys
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
sys.path.insert(1, global_path + 'Functions/')
from process_visium import visium_noh5

title_size = 12
labelout_size = 10
label_size = 8
ticks_size = 6

# Funtion for loading cell2loc data
def prepare_cell2loc(results_folder, cells): 
    # Load data
    run_name = f'{results_folder}/cell2location_map'
    adata_file = f"{run_name}/sp.h5ad"
    adata_c2l = sc.read_h5ad(adata_file)

    # Abundaces
    c2l = adata_c2l.obsm['q05_cell_abundance_w_sf']

    # Select abundances for the specific cell types
    abundances = []
    for i in cells: 
        c = adata_c2l.obsm['q05_cell_abundance_w_sf'].columns.str.contains(i)
        try: 
            sum_val = c2l[c2l.columns[c]].sum(axis = 1).values
        except:
            sum_val = c2l[c2l.columns[c]].iloc[:,0]

        abundances.append(sum_val)

    abundances = np.array(abundances)

    # Select coordinates
    coords = [adata_c2l.obs.array_row, adata_c2l.obs.array_col]

    return abundances, coords

# Function for loading cell2loc data (already processed - Heart data)
def prepare_cell2loc_processed(adatafile, cells):
    abundances = []
    adata = sc.read(adatafile)
    for i in cells: 
        abundances.append(adata.obs[i].values)
    abundances = np.array(abundances)
    coords = [adata.obs.array_row, adata.obs.array_col]
    return abundances, coords


# Function for loading CytoSPACE data
def prepare_cytospace(results_folder, cells): 
    '''
    Function to obtain the abundances and coordinates by spot of the selected cell types from the CytoSPACE results.
    Inputs:
        - results_folder: folder with the results of the CytoSPACE analysis
        - cells: list with the names of the cell types to be plotted
    '''
    # Load data
    cyto_cell2spot = pd.read_csv(f'{results_folder}/assigned_locations.csv') # Each individual cell assigned to which spot
    cyto_spot2cell = pd.read_csv(f'{results_folder}/cell_type_assignments_by_spot.csv') # For each spot, which cells types are assigned to it (quantity)
    del(cyto_spot2cell['Total cells']) # Remove total counts per spot column

   # Select abundances for the specific cell types
    abundances = []
    for i in cells: 
        sum_val = cyto_spot2cell[cyto_spot2cell.columns[cyto_spot2cell.columns.str.contains(i)]].sum(axis=1).values
        abundances.append(sum_val)
    abundances = np.array(abundances)
    
    # Select coordinates
    coordinates = cyto_cell2spot[['SpotID', 'row', 'col']].drop_duplicates()
    coords = [coordinates['row'], coordinates['col']]

    return abundances, coords


# Function for Tangram
def prepare_tangram(results_file, cells): 
    # Load data
    ad_map = sc.read(results_file)

    # Build the count matrix of number of cells per spot
    try: 
        df_spot2cell = ad_map.obs[['cell_spots', 'celltype_minor']]
        spot2cell = df_spot2cell.pivot_table(index='cell_spots', columns='celltype_minor', aggfunc=len, fill_value=0)

    except: 
        df_spot2cell = ad_map.obs[['cell_spots', 'cell_state']]
        spot2cell = df_spot2cell.pivot_table(index='cell_spots', columns='cell_state', aggfunc=len, fill_value=0)

    ## Select the coordinates
    coordinates_tg = ad_map.obs[['cell_spots', 'cell_coords_x', 'cell_coords_y']].drop_duplicates()
    coordinates_tg.sort_values(by = 'cell_spots', inplace = True)

   # Select abundances for the specific cell types
    abundances = []
    for i in cells: 
        sum_val = spot2cell[spot2cell.columns[spot2cell.columns.str.contains(i)]].sum(axis=1).values
        abundances.append(sum_val)
    abundances = np.array(abundances)

    coords = [coordinates_tg['cell_coords_x'], coordinates_tg['cell_coords_y']]

    return abundances, coords

# Funtion to convert Celltrek locations to the closest spot

def assign_loc2spot (celltrekfile, coordsfile): 
    '''
    Function to assign the location of the cells to the closest spot in the spatial data.
    Inputs:
    - celltrekfile: file with the celltrek results
    - coordsfile: file with the coordinates of the spots
    '''
    coords = pd.read_csv(coordsfile)
    coords['Unnamed: 0'] = coords['Unnamed: 0'].str.replace('.', '-')
    meta_cells_total = pd.read_csv(celltrekfile)
    spot_asig, coord_x, coord_y = [], [], []

    for i in meta_cells_total.id_raw: 
        c = meta_cells_total[meta_cells_total.id_raw == i] # select current cell
        c = c[['coord_x', 'coord_y']].values # get coordinates
        idx_min = cdist(c,coords[['imagerow', 'imagecol']].values).argsort()[0][0] # get the closest spot
        spot_asig.append(coords['Unnamed: 0'][idx_min]) # append the spot with the closest distance
        coord_x.append(coords[coords['Unnamed: 0'] == coords['Unnamed: 0'][idx_min]]['imagerow'].values[0]) # append the coordinates
        coord_y.append(coords[coords['Unnamed: 0'] == coords['Unnamed: 0'][idx_min]]['imagecol'].values[0])

    meta_cells_total['spot'] = spot_asig # add the spot to the meta_cells_total
    meta_cells_total['imagerow'] = coord_x # add the coordinates to the meta_cells_total
    meta_cells_total['imagecol'] = coord_y
    
    return meta_cells_total

def prepare_celltrek(meta_cells_total, cells):
    '''
    Function to obtain the abundances and coordinates by spot of the selected cell types from the Celltrek result mapped to the spots.
    Inputs:
        - meta_cells_total: dataframe with the celltrek results mapped to the spots
        - cells: list with the names of the cell types to be plotted
    '''
    # convert df to count matrix
    df_spot2cell = meta_cells_total[['spot', 'cell_type']]
    spot2cell = df_spot2cell.pivot_table(index='spot', columns='cell_type', aggfunc=len, fill_value=0)
    # sort the coordinates to match the order of the count matrix
    coordinates = meta_cells_total[['spot', 'imagerow', 'imagecol']].drop_duplicates()
    coordinates.sort_values(by=['spot'], inplace=True)
    # Select abundances for the specific cell types
    abundances = []
    for i in cells: 
        sum_val = spot2cell[spot2cell.columns[spot2cell.columns.str.contains(i)]].sum(axis=1).values
        abundances.append(sum_val)
    abundances = np.array(abundances)
    coords = [coordinates['imagerow'], coordinates['imagecol']]

    return abundances, coords

## Function for prepare the data for the histology figure
def prepare_patological (patient):
    ### breast cancer slide and metadata
    main_path_exp = global_path + 'Data/Breast/spatial_data/filtered_count_matrices/'+patient+'_filtered_count_matrix/'
    main_path_spatial = global_path + 'Data/Breast/spatial_data/spatial/'+patient+'_spatial/'
    library_id = 'breast_tissue'

    adata_vis = visium_noh5(main_path_exp, main_path_spatial, library_id) ## anndata object with counts, images and spatial coordinates.
    meta_patho = pd.read_csv(global_path + 'Data/Breast/spatial_data/metadata/'+patient+'_metadata.csv')

    ## Add metadata to adata_vis and asign a color to each region
    regions = meta_patho.Classification.unique()
    type2color = dict(zip(regions, sc.pl.palettes.godsnot_102))
    meta_patho['color'] = [type2color[i] for i in meta_patho.Classification.tolist()]
    adata_vis.obs = pd.merge(adata_vis.obs, meta_patho.set_index('Unnamed: 0'), left_index=True, right_index=True)
    return adata_vis, type2color

## Functions for exploring the spatial slide

## Explore the histology annotation
def histology_figure(adata_file, slide):
    adata_vis = sc.read(adata_file)
    type2color = pickle.load(open('colors/type2color_'+slide+'.pkl', 'rb'))
    adata_vis.obs['color'] = [type2color[i] for i in adata_vis.obs.annotation_final.tolist()]
    return adata_vis, type2color

## Explore the abundances per spot (Cell2loc, CytoSPACE or Tangram)
def abundances_figure(array_abundances, coords, axes, labels, title, colorbar = False, colorbar_bounds = [], save = False): 

    '''


    Function to plot the abundances of the selected cell types in the spatial coordinates.
    Inputs:
    - array_abundances: array with the abundances of the cell types. Each row is a cell type and each column is a spot (output of prepare_cell2loc or prepare_cytospace)
    - coords: list with the x and y coordinates of the spots (output of prepare_cell2loc, prepare_cytospace or prepare_tangram)
    - axes: axes of the figure
    - labels: list with the names of the cell types
    - title: title of the plot
    - colormap: if True, the plot will have a colorbar for each cell type inside the figure
    - colorbar_bounds: If colormap is True, list with the bounds for each colorbar. Each element is a list with the bounds for the colorbar of each cell type.
                       Define the for each colorbar [x, y, width, height], relative to the ax_scatter dimensions.
    - save: if filename is given, the plot will be saved in the given path

    '''

    ### Normalize the abundances of all cell types
    c1_normalized = (array_abundances[0] - min(array_abundances[0])) / (max(array_abundances[0]) - min(array_abundances[0]))
    c2_normalized = (array_abundances[1] - min(array_abundances[1])) / (max(array_abundances[1]) - min(array_abundances[1]))
    c3_normalized = (array_abundances[2] - min(array_abundances[2])) / (max(array_abundances[2]) - min(array_abundances[2]))
    c4_normalized = (array_abundances[3] - min(array_abundances[3])) / (max(array_abundances[3]) - min(array_abundances[3]))

    ### Determine values to plot (the max value for each spot)
    plot_c1 = np.maximum(c1_normalized, np.maximum(c2_normalized, np.maximum(c3_normalized, c4_normalized))) == c1_normalized
    plot_c2 = np.maximum(c1_normalized, np.maximum(c2_normalized, np.maximum(c3_normalized, c4_normalized))) == c2_normalized
    plot_c3 = np.maximum(c1_normalized, np.maximum(c2_normalized, np.maximum(c3_normalized, c4_normalized))) == c3_normalized
    plot_c4 = np.maximum(c1_normalized, np.maximum(c2_normalized, np.maximum(c3_normalized, c4_normalized))) == c4_normalized

    ################# FIGURE #################
    cm_green = plt.cm.Greens
    cm_blue = plt.cm.Blues
    cm_orange = plt.cm.Oranges
    cm_red = plt.cm.Reds

    # Plotting
    sc1 = axes.scatter(coords[1][plot_c1], -coords[0][plot_c1], s=10, 
            c=c1_normalized[plot_c1], cmap=cm_green, label=labels[0], alpha=1)

    sc2 = axes.scatter(coords[1][plot_c2], -coords[0][plot_c2], s=10,  
            c=c2_normalized[plot_c2], cmap=cm_blue, label=labels[1], alpha=1)
    
    sc3 = axes.scatter(coords[1][plot_c3], -coords[0][plot_c3], s=10,  
            c=c3_normalized[plot_c3], cmap=cm_orange, label=labels[2], alpha=1)

    sc4 = axes.scatter(coords[1][plot_c4], -coords[0][plot_c4], s=10,  
            c=c4_normalized[plot_c4], cmap=cm_red, label=labels[3], alpha=1)
    
    # Adjusting plot
    axes.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    axes.set_xlabel('')
    axes.set_ylabel('')
    axes.set_title(title, fontsize=title_size)

    if colorbar == True: 
        ############# Colorbars ################

        # Function to create an inset colorbar within the given axis
        def add_inset_colorbar(ax, mappable, label, bounds, ticks, tick_labels):
            cax = ax.inset_axes(bounds) # bounds are [x, y, width, height] relative to the ax_scatter dimensions.x and y are the coordinates, width and height are the width and height of the colorbar
            cbar = plt.colorbar(mappable, cax=cax, orientation='horizontal') # create the colorbar
            cbar.set_label(label, fontsize=label_size)
            cbar.ax.xaxis.set_ticks_position('bottom')
            cbar.ax.xaxis.set_label_position('top')
            cbar.ax.tick_params(labelsize=ticks_size, axis = 'x', length = 2, pad = 1.5)  # pad is the distance between the label and the tick and length is the length of the ticks
            cbar.set_ticks(ticks)  # Set the positions of the ticks
            cbar.set_ticklabels(tick_labels)  

        # Define the ticks of the colormap
        ticks_cancer = [0, 1]
        tick_labels_cancer = [array_abundances[0].min(), array_abundances[0].max()]
        ticks_t_cells = [0, 1]
        tick_labels_t_cells = [array_abundances[1].min(), array_abundances[1].max()]
        ticks_cafs = [0, 1]
        tick_labels_cafs = [array_abundances[2].min(),  array_abundances[2].max()]
        ticks_macrophages = [0, 1]
        tick_labels_macrophages = [array_abundances[3].min(), array_abundances[3].max()]

        # Create colorbars
        add_inset_colorbar(axes, sc1, labels[0], colorbar_bounds[0], ticks_cancer, tick_labels_cancer)
        add_inset_colorbar(axes, sc2, labels[1], colorbar_bounds[1], ticks_t_cells, tick_labels_t_cells)
        add_inset_colorbar(axes, sc3, labels[2], colorbar_bounds[2], ticks_cafs, tick_labels_cafs)
        add_inset_colorbar(axes, sc4, labels[3], colorbar_bounds[3], ticks_macrophages, tick_labels_macrophages)

    if save: 
        plt.savefig(save, dpi = 300, bbox_inches = 'tight')

## Explore the abundances per spot (CellTrek)
def celltrek_figure(file, cells, axes, labels, title): 
    meta_cells_charting = pd.read_csv(file)
    ### Define colors for each cell type
    color = []
    for i in meta_cells_charting['cell_names']: 
        if cells[0] in i: 
            color.append('green')
        elif cells[1] in i:
            color.append('blue')
        elif cells[2] in i:
            color.append('orange')
        elif cells[3] in i: 
            color.append('red')
        else:
            color.append('grey')

    meta_cells_charting['color'] = color

    ## Define coordinates for each cell type
    xs_b = meta_cells_charting[meta_cells_charting.color == 'grey'].coord_x.values
    ys_b = meta_cells_charting[meta_cells_charting.color == 'grey'].coord_y.values

    xs_c = meta_cells_charting[meta_cells_charting.color == 'green'].coord_x.values
    ys_c = meta_cells_charting[meta_cells_charting.color == 'green'].coord_y.values

    xs_t = meta_cells_charting[meta_cells_charting.color == 'blue'].coord_x.values
    ys_t = meta_cells_charting[meta_cells_charting.color == 'blue'].coord_y.values

    xs_m = meta_cells_charting[meta_cells_charting.color == 'red'].coord_x.values
    ys_m = meta_cells_charting[meta_cells_charting.color == 'red'].coord_y.values

    xs_ca = meta_cells_charting[meta_cells_charting.color == 'orange'].coord_x.values
    ys_ca = meta_cells_charting[meta_cells_charting.color == 'orange'].coord_y.values

    ## Plot them
    axes.scatter(ys_b, -xs_b, s=10, color = 'grey', alpha = 0.6, label = 'Background')
    axes.scatter(ys_t, -xs_t, s=10, color = 'blue', alpha = 0.8, label =labels[1])
    axes.scatter(ys_ca, -xs_ca, s=10, color = 'orange', alpha = 0.8, label = labels[2])
    axes.scatter(ys_m, -xs_m, s=10, color = 'red', alpha = 0.8, label = labels[3])
    axes.scatter(ys_c, -xs_c, s=10, color = 'green', alpha = 0.8, label = labels[0])

    axes.set_xlabel('')
    axes.set_ylabel('')
    axes.set_title(title, fontsize = title_size)
    axes.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)





