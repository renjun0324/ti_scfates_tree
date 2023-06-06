#!/usr/local/bin/python


import dynclipy
# task = dynclipy.main(definition_location="/share/data6/tmp/renjun/CellTrekResult/CellTrek/celltrek_pipline/ti_scfates_tree/definition.yml")
task = dynclipy.main(definition_location="")

#-----------------------------------------------------------------------------
#
#                                 create
#
#-----------------------------------------------------------------------------

import json
import sys
import os
import igraph
import scFates as scf
import warnings
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import palantir
import matplotlib.pyplot as plt
import time 
# warnings.filterwarnings("ignore")
# sc.settings.verbosity = 3
# sc.settings.logfile = sys.stdout

checkpoints = {}

# counts = palantir.io.from_csv('/home/renjun/文档/marrow_sample_scseq_counts.csv.gz')
# start_id = 'C658'
# parameters = {"n_top_genes": 1500,
#               "dim_used_for_output": "fr",
#               "dim_used_for_ppt": "fr",
#               "palantir_components_num": 10,
#               "knn_num": 30,
#               "n_eigs": 4,
#               "tree_method":'ppt',
#               "ppt_nodes": 200,
#               "ppt_lambda": 1,
#               "ppt_sigma": 0.1,
#               "ppt_nsteps": 50,
#               "n_map": 1,
#               "epg_lambda": 0.01,
#               "epg_mu": 0.1,
#               "epg_initnodes": 2,
#               "epg_extend_leaves": False }
# dim_used_for_output = parameters["dim_used_for_output"]
# dim_used_for_ppt = parameters["dim_used_for_ppt"]

counts = task["counts"]
parameters = task["parameters"]
start_id = task["priors"]["start_id"]

# if dim_used_for_output=='fr':
#   dim_used_for_output = 'draw_graph_fr'
# 
# if dim_used_for_ppt=='fr':
#   dim_used_for_ppt = 'draw_graph_fr'
  
if isinstance(start_id, list):
  start_id = start_id[0]

checkpoints["method_afterpreproc"] = time.time()

#-----------------------------------------------------------------------------
#
#                                 running
#
#-----------------------------------------------------------------------------

# 1. Load, normalize and log-transform count data
# norm_df = sc.pp.normalize_per_cell(counts,copy=True)
# norm_df = palantir.preprocess.log_transform(norm_df)
# adata = sc.AnnData(norm_df)
# adata.write("/public/home/renjun/x.h5")
adata = anndata.AnnData(counts)
adata = sc.pp.normalize_per_cell(adata, copy=True)
palantir.preprocess.log_transform(adata)
# adata.write("/public/home/renjun/1.h5")

# 2. Perform PCA on highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=parameters['n_top_genes'], flavor='cell_ranger')
sc.pp.pca(adata)
# print("11111")

# 3. X_Palantir
pca_projections = pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections,
                                           n_components=parameters["palantir_components_num"],
                                           knn=parameters["knn_num"])
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=parameters["n_eigs"])
adata.obsm["X_palantir"]=ms_data.values
# print("22222")

# 4. Generate embedding from the multiscale diffusion space
sc.pp.neighbors(adata,n_neighbors=parameters["knn_num"], use_rep="X_pca")
adata.obsm["X_pca2d"]=adata.obsm["X_pca"][:,:2]
sc.tl.draw_graph(adata,init_pos='X_pca2d')
# print("33333")
adata.write("/public/home/renjun/2.h5")

# 5. Tree learning with SimplePPT
t = adata.obsm['X_palantir'].shape[0]
if t<200:
  parameters["ppt_nodes"] = t-20

scf.tl.tree(adata,
            use_rep='palantir',
            device="cpu",
            seed=723,
            method=parameters["tree_method"],
            Nodes=parameters["ppt_nodes"],
            ppt_lambda=parameters["ppt_lambda"],
            ppt_sigma=parameters["ppt_sigma"],
            ppt_nsteps=parameters["ppt_nsteps"],
            epg_lambda=parameters["epg_lambda"],
            epg_mu=parameters["epg_mu"],
            epg_initnodes=parameters["epg_initnodes"],
            epg_extend_leaves =parameters["epg_extend_leaves"])
            
# 6. Selecting a root and computing pseudotime
scf.tl.root(adata, root=adata.uns['graph']['tips'][0])

if parameters["tree_method"]=='epg':
  scf.tl.convert_to_soft(adata, sigma = parameters["ppt_sigma"], lam = parameters["ppt_lambda"])
  
try:
  scf.tl.pseudotime(adata, n_jobs=1, n_map=parameters["n_map"], seed=723)
except:
  print("scf.tl.pseudotime error")
  
  
tmp = adata.obs['milestones']
i = np.where(tmp.index==start_id)
start_index = tmp[i[0][0]]
scf.tl.root(adata, int(start_index))

try:
  scf.tl.pseudotime(adata, n_jobs=1, n_map=parameters["n_map"], seed=723)
except:
  print("scf.tl.pseudotime error")
  
  
# 7. plot
# sc.pl.draw_graph(adata,color="CD34",color_map="RdBu_r",save=".png")
# scf.pl.graph(adata, basis = dim_used_for_output, save=".png")
# scf.pl.trajectory(adata, basis = dim_used_for_output, save = ".png")

checkpoints["method_aftermethod"] = time.time() 

#-----------------------------------------------------------------------------
#
#                                  output
#
#-----------------------------------------------------------------------------

# dimred
dimred = pd.DataFrame([x for x in adata.obsm["X_draw_graph_fa"].T]).T
dimred.columns = ["comp_" + str(i+1) for i in range(dimred.shape[1])]
dimred["cell_id"] = adata.obs.index

df = adata.uns['graph']['pp_seg']
milestone_network = pd.DataFrame({
  "from": df["from"].astype(str),
  "to": df["to"].astype(str),
  "length": df["d"].astype(float),
  "directed": True
})
ind = np.where(milestone_network['from']!=milestone_network['to'])
milestone_network = milestone_network.iloc[ind]

# grouping
grouping = pd.DataFrame({
  "cell_id": adata.obs.index,
  "group_id": adata.obs.milestones.astype(str)
})

# add nan judge
ind = np.where(grouping['group_id']!="nan")[0]
grouping = grouping.iloc[ind]
cell_ids = adata.obs.index[ind]
dimred = dimred.iloc[ind]


# output
dataset = dynclipy.wrap_data(cell_ids = cell_ids)
dataset.add_cluster_graph(
  grouping = grouping,
  milestone_network = milestone_network
)

dataset.add_dimred(dimred=dimred)
dataset.add_timings(timings = checkpoints)
dataset.write_output(task["output"])

# dimred_milestones
# emb = adata.obsm["X_"+dim_used_for_output]
# emb = emb[:, :2]
# R = adata.obsm["X_R"]
# proj = (np.dot(emb.T, R) / R.sum(axis=0)).T
# dimred_milestones = pd.DataFrame(proj)
# dimred_milestones = pd.DataFrame(adata.uns["graph"]['F']).transpose()
# dimred_milestones.columns = ["comp_" + str(i+1) for i in range(dimred_milestones.shape[1])]

# 200个采样点的邻接矩阵
# A = adata.uns['graph']['B']
# g = igraph.Graph.Adjacency((A > 0).tolist())
# g = igraph.Graph.Adjacency(A.tolist())
# df = g.get_edge_dataframe()
# df.index.name=None
# df.columns = ["from", "to"]

# .uns['ppt'], dictionnary containing inferred tree.
# adata.obsm['X_R'] soft assignment of cells to principal points.
# adata.uns["graph"]['B'] adjacency matrix of the principal points.
# adata.uns["graph"]['F'] coordinates of principal points in representation space.
# adata.write("/public/home/renjun/x2.h5")

# basis = dim_used_for_output
# graph = adata.uns["graph"]
# basis = get_basis(adata, basis)
# emb = adata.obsm[f"X_{basis}"]
# emb = emb[:, :2]
# R = adata.obsm["X_R"]
# proj = (np.dot(emb.T, R) / R.sum(axis=0)).T
# B = graph["B"]
# ax = sc.pl.embedding( adata, basis=basis, show=False)
# al = np.array( igraph.Graph.Adjacency((B > 0).tolist(), mode="undirected").get_edgelist() )
# segs = al.tolist()
# vertices = proj.tolist()
# lines = [[tuple(vertices[j]) for j in i] for i in segs]
# lc = matplotlib.collections.LineCollection( lines, colors="k", linewidths=1, alpha=1, rasterized=True )
# ax.add_collection(lc)
# ax.scatter( proj[:, 0], proj[:, 1], s=1, c="k", alpha=1, rasterized=True )
# savefig_or_show("graph", show=False, save="1.png")

# import matplotlib.pyplot as plt
# import numpy as np
# plt.plot(dimred["comp_1"], dimred["comp_2"], 'o')
# plt.plot(dimred_milestones["comp_1"], dimred_milestones["comp_2"], 'o')
# plt.plot(dimred_milestones[:,0],dimred_milestones[:,1],'o')
# plt.show()
# plt.savefig(dim_used_for_output+'_mil.png')
# plt.close()

# import os
# import rpy2.robjects as robjects
# from rpy2.robjects import pandas2ri
# from rpy2.robjects.packages import importr
# pandas2ri.activate()
# readRDS = robjects.r['readRDS']
# data = readRDS('~/x.rds')
# rdsKeys = data.dimnames
