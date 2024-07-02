
## For adding metadata or umap coordinates. Recommended to run in command line as necessary.

# Import packages and csv files with relevant metadata
import numpy as np
import pandas as pd
umap = pd.read_csv('d14_rna_umap_embeddings_ver2.csv')
umap = umap.rename(columns = {'Unnamed: 0':'Cell ID'})

meta= pd.read_csv('/projects/p20023/Viriya/analysis/foxa2/multiome/irina/d2_metadata_small.csv')
meta = meta.rename(columns = {'Unnamed: 0':'Cell ID'})

# Import loom file containing the spliced-unspliced mRNA matricies and metadata
import loompy
ds = loompy.connect("d14_wclusters.loom")

# If necessary, for new loom files; rename the cell IDs in the loom file according to what is in the metadata.
#cell_ids =  ds.ca.CellID
#cell_ids = [i.replace('d2:','') for i in cell_ids]
#cell_ids = [i.replace('x','-1') for i in cell_ids]
#ds.ca.CellID = cell_ids


# Get cell ids from ds and from umap
cell_ids =  ds.ca.CellID
cell_ids = pd.DataFrame(cell_ids)
cell_ids = cell_ids.rename(columns = {0:'Cell ID'})

# Add cell ids to the umap dataframe
umap_ordered = cell_ids.merge(umap, on = "Cell ID")

# Filter for which cells are included in the metadat
cells_filt = pd.merge(umap_ordered['Cell ID'], cell_ids)
in_list = np.isin(ds.ca.CellID, cells_filt) # boolean array for which cell IDs in ds are left after filtering


## Add dummy empty spaces in the umap data and metadata for the cells that are not filtered
# get cell IDs for those that get filtered out
not_in_list = np.where(in_list == False)
not_in_list = np.array(not_in_list)
not_in_list = not_in_list.flatten()
not_in_list = not_in_list.tolist()
cell_ids_list = cell_ids['Cell ID']
not_cells_filt = pd.DataFrame(cell_ids_list[not_in_list])

# add 0s for the filtered out cell IDs
not_cells_filt['umap_rna_alone_30_1'] = 0 #cannot use NaNs; I think this is what causes an error in score_cv_vs_mean of ValueError: Input y contains infinity or a value too large for dtype('float64') when running velocyto
not_cells_filt['umap_rna_alone_30_2'] = 0
not_cells_filt = not_cells_filt.reset_index()

# concatenate and reorder umap with the filtered out cell IDs
umap_ordered = umap_ordered.append(not_cells_filt, ignore_index=True)
umap_ordered = cell_ids.merge(umap_ordered, on = "Cell ID")
umap_ordered = umap_ordered.drop(['index'], axis=1)

# do the same for the metadata
cell_ids_list = cell_ids['Cell ID']
not_cells_filt = pd.DataFrame(cell_ids_list[not_in_list])
not_cells_filt['klk3_status'] = 0
not_cells_filt['AR_signature1'] = 0
not_cells_filt['Luminal_signature1'] = 0
not_cells_filt['NE_signature1'] = 0
not_cells_filt = not_cells_filt.reset_index()


meta_ordered = cell_ids.merge(meta, on = "Cell ID")
meta_ordered = meta_ordered.append(not_cells_filt, ignore_index=True)
meta_ordered = cell_ids.merge(meta_ordered, on = "Cell ID")
meta_ordered = meta_ordered.drop(['index'], axis=1)

# Add column attribute, of which cells are those that are filtered
ds.ca.filtered_atac = in_list

# Add metadata and umap embeddings to filtered
umap_ordered = umap_ordered.iloc[:,1:]
umap_ordered = umap_ordered.to_numpy()
ds.ca.umap_2 = umap_ordered

#meta_ordered = meta_ordered.iloc[:,1:]
#ds.ca.klk3_status = meta_ordered.klk3_status.to_numpy()
#ds.ca.AR_signature = meta_ordered.AR_signature1.to_numpy()
#ds.ca.Luminal_signature = meta_ordered.Luminal_signature1.to_numpy()
#ds.ca.NE_signature = meta_ordered.NE_signature1.to_numpy()

#loompy.create("d14_wclusters_2.loom", ds.layers, ds.ra, ds.ca) <--- no need for this step, unless you want another file:  all additions to the loom file are saved when you close it.

ds.close()


## Similar code to above, but for isolating the heterogenous population in d21.
import numpy as np
import pandas as pd

import loompy
ds = loompy.connect("d14_wclusters.loom")

# Get cell ids
cell_ids =  ds.ca.CellID
cell_ids = pd.DataFrame(cell_ids)
cell_ids = cell_ids.rename(columns = {0:'Cell ID'})

#Get which RNA-seq umap coordinates are negative
umap = ds.ca.umap
umap = pd.DataFrame(umap)
negatives = umap[(umap < 0).all(axis=1)]
in_list = np.isin(umap, negatives)
in_list = in_list[:,0]

ds.ca.filtered_heterogeneous = in_list

ds.close()