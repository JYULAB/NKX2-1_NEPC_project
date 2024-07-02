
### Estimate RNA velocity ####
import velocyto as vcy
import matplotlib.pyplot as plt
import numpy as np

# Load data
vlm = vcy.VelocytoLoom("d14_wclusters.loom")
#vlm = vcy.VelocytoLoom("/projects/b1042/YuLab/Viriya/multiome/d14/velocyto/d14.loom")

# Normalize
vlm.normalize("S", size=True, log=True)
vlm.S_norm  # contains log normalized
 

#vlm.plot_fractions()
#plt.savefig('d21_plotfractions.png')
#plt.clf()

# Save analysis results <- this doesn't work: AttributeError: 'VelocytoLoom' object has no attribute 'dump_hdf5'
#vlm.dump_hdf5("d14_velocyto_analysis")


# Filtering
#vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5))
#instead: filter cells by the boolean array that is in the loom file
vlm.filter_cells(bool_array=vlm.ca["filtered"] == 1)

# Cluster annotations
import matplotlib.pyplot as plt
# colors to use:  #000080 (0, 0, 128) and #FFA500 (255, 165, 0)
# code for colormap:
from matplotlib import colors
new_cmap = colors.LinearSegmentedColormap.from_list('new_cmap', ["#000080","#FFA500"], N=2)
# cluster names: klk3_status, AR_signature, Luminal_signature, NE_signature. Use plt.cm.PuOr for signatures, and .Set1 for klk3
vlm.set_clusters(vlm.ca["klk3_status"], colormap = new_cmap)


# Select the genes that are expressed above a threshold of total number of molecules in any of the clusters
vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)
vlm.filter_genes(by_detection_levels=True)

# Feature selection
vlm.score_cv_vs_mean(3000, plot=True, max_expr_avg=35)
vlm.filter_genes(by_cv_vs_mean=True)
#plt.savefig('d14_PCA.png')
plt.clf()


# Normalize by size (total molecule count)
vlm._normalize_S(relative_size=vlm.S.sum(0),
             target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0),
             target_size=vlm.U.sum(0).mean())


# Prepare for gamma fit (smooth data using kNN neighbors pooling approach)
vlm.perform_PCA()

vlm.knn_imputation(n_pca_dims=20, balanced=True, n_jobs=4) # unsure about sight and maxl for this step; had to switch b_sight to be less than n_samples

# Fit gamma to every gene that survived the filtering step
vlm.fit_gammas()


# Fit visualization
#vlm.plot_phase_portraits(["AR", "KLK3"])


# Calculate velocity and predict future state of cells
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)


# Project velocities onto embeddings
## For TSNE -- not used
#from sklearn.manifold import TSNE
#bh_tsne = TSNE()
#vlm.ts = bh_tsne.fit_transform(vlm.pcs[:, :25]) # takes a long time for big dataset

## Using UMAP coordinates added to the loom file -- see scratch.py for how
plt.figure(figsize=(8, 8), dpi=80)
plt.scatter(vlm.ca["umap_2"][:,0], vlm.ca["umap_atac"][:,1])
plt.savefig('d14_UMAP_atac.png')
plt.clf()

vlm.ts = vlm.ca["umap_2"]
vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1, knn_random=True, sampled_fraction=0.5)
vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)
vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=300)

plt.figure(None,(9,9))

## Project velocities onto the UMAP
vlm.plot_grid_arrows(quiver_scale=0.2,
                    scatter_kwargs_dict={"alpha":0.7, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True}, 
		    min_mass=2.9, angles='xy', scale_units='xy',
                    headaxislength=2.75, headlength=5, headwidth=4.8,
                    plot_random = False, scale_type="absolute")
plt.savefig('d14_velocities_umap_klk3_30.png')
#plt.savefig('d14_velocities_umap_klk3.pdf', format="pdf")
plt.clf()

## Plot expression of gene as color gradient on UMAP
plt.figure(None,(9,9))
vlm.plot_expression_as_color("KLK3", s=20)
plt.savefig('d14_klk3_atc.png')
plt.clf()

