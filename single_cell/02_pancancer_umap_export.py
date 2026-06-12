"""Export UMAP coordinates and cell-type labels from the pan-cancer AnnData.

The R companion script (02_pancancer_umap_plot.R) reads the resulting CSV to draw a
rasterised UMAP. Run this after the AnnData `adata` has been pre-processed
(QC -> Harmony -> neighbours -> UMAP) and annotated.

Output: data/umap_metadata.csv  (columns: UMAP_1, UMAP_2, major_celltypes_tier2)
"""

import pandas as pd

# `adata` is the pre-processed pan-cancer AnnData, already in the session.
umap = adata.obsm["X_umap"]
umap_df = pd.DataFrame(umap, columns=["UMAP_1", "UMAP_2"])
umap_df["major_celltypes_tier2"] = adata.obs["major_celltypes_tier2"].values

umap_df.to_csv("../data/umap_metadata.csv", index=False)
print("Exported: ../data/umap_metadata.csv")


# ---------------------------------------------------------------------------
# Alternative: draw the UMAP directly in scanpy (same palette as the R script).
# ---------------------------------------------------------------------------
# import scanpy as sc
# import matplotlib.pyplot as plt
#
# celltype_colors = {
#     "T": "#FDBF6F", "Tumor": "#A6CEE3", "Macrophage": "#B2DF8A", "EC": "#FB9A99",
#     "Mesenchymal": "#1F78B4", "B": "#33A02C", "Neutrophil": "#FF7F00", "NK": "#FF7F00",
#     "Monocyte": "#CAB2D6", "DC": "#6A3D9A", "Mast": "#B15928", "ILC": "#E31A1C",
#     "N": "#FDBF6F", "GMP": "#B2DF8A", "Mesothelial": "#999999",
# }
# with plt.rc_context({"figure.figsize": (6, 6)}):
#     sc.pl.umap(adata, color=["major_celltypes_tier2"], neighbors_key="neighbors_harmony",
#                frameon=False, legend_fontsize=5, legend_fontoutline=1, palette=celltype_colors)
