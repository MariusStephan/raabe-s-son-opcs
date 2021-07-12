############## Python script for Velocyto

#### Filter cells
import velocyto as vcy
import scanpy as sc
import re
import math
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
#%load_ext rpy2.ipython
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt

# Load dataset
OPC = anndata.read_loom("~/Velocyto/OPC_merged.loom")

# Load Seurat pipeline result metadata
cellID_obs = pd.read_csv("~/Project/OPCmerged_cellIDs.csv")
umap = pd.read_csv("~/Project/OPCmerged_UMAP.csv")
clusters = pd.read_csv("~/Project/OPCmerged_clusters.csv")
samples =  pd.read_csv("~/Project/OPCmerged_groups.csv")

# Adapt sample ID prefix to the R metadata
OPC.obs.index = OPC.obs.index.str.replace("NPC_S1:","Pre_")
OPC.obs.index = OPC.obs.index.str.replace("OPC1_S_S5:","S_")
OPC.obs.index = OPC.obs.index.str.replace("OPC2_SON_S6:","SON_")

# Filter cells
OPC = OPC[np.isin(OPC.obs.index,cellID_obs["x"])]

##### Add UMAP coordinates and cluster IDs
##########################################

# Cast index to dataframe and rename the column
OPC_index = pd.DataFrame(OPC.obs.index)
OPC_index = OPC_index.rename(columns = {0:'Cell ID'})

### Add UMAP coordinates first
# Rename ID column for UMAP as well
umap = umap.rename(columns = {'Unnamed: 0':'Cell ID'})

# Reorder umap following Cell ID
umap_ordered = OPC_index.merge(umap, on = "Cell ID")

# Strip Cell ID column from UMAP coordinates and add them to the datasets
umap_ordered = umap_ordered.iloc[:,1:]
OPC.obsm['X_umap'] = umap_ordered.values

### Add cluster IDs
# Rename, order, strip and add
clusters = clusters.rename(columns = {'Unnamed: 0':'Cell ID'})
clusters_ordered = OPC_index.merge(clusters, on = "Cell ID")
clusters_ordered = clusters_ordered.iloc[:,1:]
clusters_ordered.index = OPC.obs.index
clusters_ordered.x = clusters_ordered.x.astype('category') # Cluster IDs are stored as categorical in scVelo
OPC.obs['clusters'] = clusters_ordered.x

### Add cluster colors
clusterPalette = pd.read_csv("~/Project/OPCmerged_clusterPalette.csv")
clusters_cats = pd.DataFrame(OPC.obs.clusters.cat.categories)
clusters_cats.columns = ['Cluster ID']
clusterPalette = clusterPalette.rename(columns = {'Unnamed: 0':'Cluster ID'})
clusterPalette_ordered = clusters_cats.merge(clusterPalette, on = "Cluster ID")
OPC.uns['clusters_colors'] = clusterPalette_ordered.x.tolist()

### Add group IDs
samples = samples.rename(columns = {'Unnamed: 0':'Cell ID'})
samples_ordered = OPC_index.merge(samples, on = "Cell ID")
samples_ordered = samples_ordered.iloc[:,1:]
samples_ordered.index = OPC.obs.index
samples_ordered.x = samples_ordered.x.astype('category')
OPC.obs['samples'] = samples_ordered.x

### Add group colors
# Rename, order, strip and add
samplePalette = pd.read_csv("~/Project/OPCmerged_groupPalette.csv")
samples_cats = pd.DataFrame(OPC.obs.samples.cat.categories)
samples_cats.columns = ['sample ID']
samplePalette = samplePalette.rename(columns = {'Unnamed: 0':'sample ID'})
samplePalette_ordered = samples_cats.merge(samplePalette, on = "sample ID")
OPC.uns['samples_colors'] = samplePalette_ordered.x.tolist()

### Some gene labels are duplicated (Ensembl IDs are still unique!!)
OPC.var_names_make_unique()

#######
#### Split sample datasets
####################
cellID_obs = OPC.obs.index


cellID_obs_Pre = cellID_obs[OPC.obs.index.str.contains("Pre_")]
cellID_obs_S = cellID_obs[OPC.obs.index.str.contains("S_")]
cellID_obs_SON = cellID_obs[OPC.obs.index.str.contains("SON_")]

Pre = OPC[np.isin(OPC.obs.index, cellID_obs_Pre)]
S = OPC[np.isin(OPC.obs.index, cellID_obs_S)]
SON = OPC[np.isin(OPC.obs.index, cellID_obs_SON)]


#### Run Velocyto pipeline
##########################################
scv.pp.filter_and_normalize(Pre)
scv.pp.moments(Pre)

# Dynamic model
scv.tl.recover_dynamics(Pre)
scv.tl.velocity(Pre, mode = "dynamical")
scv.tl.velocity_graph(Pre)
scv.tl.score_genes_cell_cycle(Pre)
scv.tl.latent_time(Pre)

# Save model to file
Pre.write('Velocyto/Pre_dynamicModel.h5ad', compression = 'gzip')


scv.pp.filter_and_normalize(S)
scv.pp.moments(S)

# Dynamic model
scv.tl.recover_dynamics(S)
scv.tl.velocity(S, mode = "dynamical")
scv.tl.velocity_graph(S)
scv.tl.score_genes_cell_cycle(S)
scv.tl.latent_time(S)

# Save model to file
S.write('Velocyto/S_dynamicModel.h5ad', compression = 'gzip')

scv.pp.filter_and_normalize(SON)
scv.pp.moments(SON)

# Dynamic model
scv.tl.recover_dynamics(SON)
scv.tl.velocity(SON, mode = "dynamical")
scv.tl.velocity_graph(SON)
scv.tl.score_genes_cell_cycle(SON)
scv.tl.latent_time(SON)

# Save model to file
SON.write('Velocyto/SON_dynamicModel.h5ad', compression = 'gzip')

#### Plotting
###########################################

OPC = scv.read('~/Velocyto/OPCmerged_dynamicModel.h5ad')

pp = PdfPages('OPC_Dimplot_Vectors.pdf')
kwargs = dict(linewidth=0)
pp.savefig(scv.pl.velocity_embedding(OPC, arrow_length=3, arrow_size=3, dpi=120, legend_loc = 'on data',**kwargs))
pp.savefig(scv.pl.velocity_embedding(OPC, arrow_length=3, arrow_size=3, dpi=120, color = 'samples', legend_loc = 'on data',**kwargs))
pp.close()

pp = PdfPages('OPC_CellCycle.pdf')
pp.savefig(scv.pl.scatter(OPC, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[2, 98]),right_margin = 5, legend_loc = 'on data')
pp.close()

pp = PdfPages('OPC_Latent_Time.pdf')
pp.savefig(scv.pl.scatter(OPC, color='latent_time', color_map='gnuplot', size=80))
pp.close()

expression_feats = ["ASCL1", "BMI1", "HES1", "MSI1", "NES", "PAX6", "PROM1", "SOX1", "SOX2", "SOX3", "ST8SIA1", "PDGFRA", "CSPG4", "CD9", "CNP", "GPR17", "PTPRZ1", "MYT1", "OLIG1", "OLIG2", "ID2", "ID4", "SOX9", "SOX10", "CLDN11", "GALC", "MAG", "MAL", "MBP", "MOBP", "MOG", "MYRF", "NKX6-2", "NKX2-2", "PLP1", "SMARCA4", "SOX17", "MKI67", "KCTD13", "ZNF488", "ZNF536", "MALAT1", "STMN2", "NEFM", "NEFL", "SNAP25", "NEUROD1", "NEUROD4", "VIM", "SPARC", "SLC6A4", "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B", "CDKN1C"]
pp = PdfPages('OPC_Expression_clustercolors_v02.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.velocity(OPC, feature, add_outline=True))

pp.close()

pp = PdfPages('OPC_Expression_samplecolors_v02.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.velocity(OPC, feature, color = 'samples', add_outline=True))

pp.close()

pp = PdfPages('OPC_LatentTimeline_clustercolors_v01.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.scatter(OPC, x='latent_time', y=feature, add_outline=True, legend_loc = 'right margin'))

pp.close()

pp = PdfPages('OPC_LatentTimeline_samplecolors_v01.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.scatter(OPC, x='latent_time', y=feature,color='samples', add_outline=True, legend_loc = 'right margin'))

pp.close()

######################### Pre

Pre = scv.read('~/97.MS_NPC_Pre_RNA-Velocity/Velocyto/Premerged_dynamicModel.h5ad')

pp = PdfPages('Pre_Dimplot_Vectors.pdf')
kwargs = dict(linewidth=0)
pp.savefig(scv.pl.velocity_embedding(Pre, arrow_length=3, arrow_size=3, dpi=120, legend_loc = 'on data',**kwargs))
pp.savefig(scv.pl.velocity_embedding(Pre, arrow_length=3, arrow_size=3, dpi=120, color = 'samples', legend_loc = 'on data',**kwargs))
pp.close()

pp = PdfPages('Pre_CellCycle.pdf')
pp.savefig(scv.pl.scatter(Pre, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[2, 98]),right_margin = 5, legend_loc = 'on data')
pp.close()

pp = PdfPages('Pre_Latent_Time.pdf')
pp.savefig(scv.pl.scatter(Pre, color='latent_time', color_map='gnuplot', size=80))
pp.close()

expression_feats = ["ASCL1", "BMI1", "HES1", "MSI1", "NES", "PAX6", "PROM1", "SOX1", "SOX2", "SOX3", "ST8SIA1", "PDGFRA", "CSPG4", "CD9", "CNP", "GPR17", "PTPRZ1", "MYT1", "OLIG1", "OLIG2", "ID2", "ID4", "SOX9", "SOX10", "CLDN11", "GALC", "MAG", "MAL", "MBP", "MOBP", "MOG", "MYRF", "NKX6-2", "NKX2-2", "PLP1", "SMARCA4", "SOX17", "MKI67", "KCTD13", "ZNF488", "ZNF536", "MALAT1", "STMN2", "NEFM", "NEFL", "SNAP25", "NEUROD1", "NEUROD4", "VIM", "SPARC", "SLC6A4", "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B", "CDKN1C"]
pp = PdfPages('Pre_Expression_clustercolors_v02.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.velocity(Pre, feature, add_outline=True))

pp.close()

pp = PdfPages('Pre_Expression_samplecolors_v02.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.velocity(Pre, feature, color = 'samples', add_outline=True))

pp.close()

pp = PdfPages('Pre_LatentTimeline_clustercolors_v01.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.scatter(Pre, x='latent_time', y=feature, add_outline=True, legend_loc = 'right margin'))

pp.close()

pp = PdfPages('Pre_LatentTimeline_samplecolors_v01.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.scatter(Pre, x='latent_time', y=feature,color='samples', add_outline=True, legend_loc = 'right margin'))

pp.close()

#################################### S

S = scv.read('~/97.MS_NPC_S_RNA-Velocity/Velocyto/Smerged_dynamicModel.h5ad')

pp = PdfPages('S_Dimplot_Vectors.pdf')
kwargs = dict(linewidth=0)
pp.savefig(scv.pl.velocity_embedding(S, arrow_length=3, arrow_size=3, dpi=120, legend_loc = 'on data',**kwargs))
pp.savefig(scv.pl.velocity_embedding(S, arrow_length=3, arrow_size=3, dpi=120, color = 'samples', legend_loc = 'on data',**kwargs))
pp.close()

pp = PdfPages('S_CellCycle.pdf')
pp.savefig(scv.pl.scatter(S, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[2, 98]),right_margin = 5, legend_loc = 'on data')
pp.close()

pp = PdfPages('S_Latent_Time.pdf')
pp.savefig(scv.pl.scatter(S, color='latent_time', color_map='gnuplot', size=80))
pp.close()

expression_feats = ["ASCL1", "BMI1", "HES1", "MSI1", "NES", "PAX6", "PROM1", "SOX1", "SOX2", "SOX3", "ST8SIA1", "PDGFRA", "CSPG4", "CD9", "CNP", "GPR17", "PTPRZ1", "MYT1", "OLIG1", "OLIG2", "ID2", "ID4", "SOX9", "SOX10", "CLDN11", "GALC", "MAG", "MAL", "MBP", "MOBP", "MOG", "MYRF", "NKX6-2", "NKX2-2", "PLP1", "SMARCA4", "SOX17", "MKI67", "KCTD13", "ZNF488", "ZNF536", "MALAT1", "STMN2", "NEFM", "NEFL", "SNAP25", "NEUROD1", "NEUROD4", "VIM", "SPARC", "SLC6A4", "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B", "CDKN1C"]
pp = PdfPages('S_Expression_clustercolors_v02.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.velocity(S, feature, add_outline=True))

pp.close()

pp = PdfPages('S_Expression_samplecolors_v02.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.velocity(S, feature, color = 'samples', add_outline=True))

pp.close()

pp = PdfPages('S_LatentTimeline_clustercolors_v01.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.scatter(S, x='latent_time', y=feature, add_outline=True, legend_loc = 'right margin'))

pp.close()

pp = PdfPages('S_LatentTimeline_samplecolors_v01.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.scatter(S, x='latent_time', y=feature,color='samples', add_outline=True, legend_loc = 'right margin'))

pp.close()


##################################### SON

SON = scv.read('~/97.MS_NPC_SON_RNA-Velocity/Velocyto/SON_dynamicModel.h5ad')

pp = PdfPages('SON_Dimplot_Vectors.pdf')
kwargs = dict(linewidth=0)
pp.savefig(scv.pl.velocity_embedding(SON, arrow_length=3, arrow_size=3, dpi=120, legend_loc = 'on data',**kwargs))
pp.savefig(scv.pl.velocity_embedding(SON, arrow_length=3, arrow_size=3, dpi=120, color = 'samples', legend_loc = 'on data',**kwargs))
pp.close()

pp = PdfPages('SON_CellCycle.pdf')
pp.savefig(scv.pl.scatter(SON, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[2, 98]),right_margin = 5, legend_loc = 'on data')
pp.close()

pp = PdfPages('SON_Latent_Time.pdf')
pp.savefig(scv.pl.scatter(SON, color='latent_time', color_map='gnuplot', size=80))
pp.close()

expression_feats = ["ASCL1", "BMI1", "HES1", "MSI1", "NES", "PAX6", "PROM1", "SOX1", "SOX2", "SOX3", "ST8SIA1", "PDGFRA", "CSPG4", "CD9", "CNP", "GPR17", "PTPRZ1", "MYT1", "OLIG1", "OLIG2", "ID2", "ID4", "SOX9", "SOX10", "CLDN11", "GALC", "MAG", "MAL", "MBP", "MOBP", "MOG", "MYRF", "NKX6-2", "NKX2-2", "PLP1", "SMARCA4", "SOX17", "MKI67", "KCTD13", "ZNF488", "ZNF536", "MALAT1", "STMN2", "NEFM", "NEFL", "SNAP25", "NEUROD1", "NEUROD4", "VIM", "SPARC", "SLC6A4", "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B", "CDKN1C"]
pp = PdfPages('SON_Expression_clustercolors_v02.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.velocity(SON, feature, add_outline=True))

pp.close()

pp = PdfPages('SON_Expression_samplecolors_v02.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.velocity(SON, feature, color = 'samples', add_outline=True))

pp.close()

pp = PdfPages('SON_LatentTimeline_clustercolors_v01.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.scatter(SON, x='latent_time', y=feature, add_outline=True, legend_loc = 'right margin'))

pp.close()

pp = PdfPages('SON_LatentTimeline_samplecolors_v01.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.scatter(SON, x='latent_time', y=feature,color='samples', add_outline=True, legend_loc = 'right margin'))

pp.close()

### Without Neurons

cellID_obs = pd.read_csv("~/Project/OPCmerged_cellIDs_noNeuronalCells.csv")
noNeurons = OPC[np.isin(OPC.obs.index,cellID_obs["x"])]

pp = PdfPages('Expression_noNeurons_clustercolors_v02.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.velocity(noNeurons, feature, add_outline=True))

pp.close()

pp = PdfPages('Expression_noNeurons_samplecolors_v02.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.velocity(noNeurons, feature, color = 'samples', add_outline=True))

pp.close()

pp = PdfPages('LatentTimeline_noNeurons_clustercolors_v01.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.scatter(noNeurons, x='latent_time', y=feature, add_outline=True, legend_loc = 'right margin'))

pp.close()

pp = PdfPages('LatentTimeline_noNeurons_samplecolors_v01.pdf')
for feature in expression_feats:
 pp.savefig(scv.pl.scatter(noNeurons, x='latent_time', y=feature,color='samples', add_outline=True, legend_loc = 'right margin'))

pp.close()

### Run PAGA on individual Datasets

scv.tl.paga(Pre, groups='clusters',use_time_prior='latent_time')
scv.tl.paga(S, groups='clusters',use_time_prior='latent_time')
scv.tl.paga(SON, groups='clusters',use_time_prior='latent_time')
scv.tl.paga(OPC, groups='clusters',use_time_prior='latent_time')

pp = PdfPages('Pagaplot_v02.pdf')
pp.savefig(scv.pl.paga(OPC, basis='umap', size=50, alpha=.15, min_edge_width=2, node_size_scale=3, legend_loc='on data',title = "All samples"))
pp.savefig(scv.pl.paga(Pre, basis='umap', size=50, alpha=.15, min_edge_width=2, node_size_scale=3, legend_loc='on data',title = "NI Day+0"))
pp.savefig(scv.pl.paga(S, basis='umap', size=50, alpha=.15, min_edge_width=2, node_size_scale=3, legend_loc='on data',title = "S Day+15"))
pp.savefig(scv.pl.paga(SON, basis='umap', size=50, alpha=.15, min_edge_width=2, node_size_scale=3, legend_loc='on data',title = "SON Day+15"))
pp.close()


### Export latent time for import in Seurat
df = scv.DataFrame(OPC.obs[['velocity_pseudotime','latent_time']])
df.to_csv('OPCmerge_pseudotime.tsv', sep='\t')
