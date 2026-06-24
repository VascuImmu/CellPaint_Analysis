"""
Cell Painting Pipeline — generated from Cell_Paint_3.ipynb
All parameters are defined in the CONFIG block at the top.
Edit them here or use the GUI to regenerate.
"""

import os
import re
import sys
from pathlib import Path

from metadata_utils import plot_quality_control

sys.path.append(str(Path(__file__).resolve().parent))

from tkinter import FALSE
import numpy as np
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.cluster import AgglomerativeClustering
from pycytominer import normalize, feature_select
from metadata_utils import parse_value, compute_alpha_per_group, assign_colors, plot_well_map



# =============================================================================
# CONFIG — edit all parameters here
# =============================================================================

if len(sys.argv) < 2:
    print("Usage: python cell_paint_pipeline.py <experiment_folder> <bool_new>")
    print("  e.g. python cell_paint_pipeline.py /Volumes/KINGSTON/Nico_data/Cellpaint/20260319 True")
    sys.exit(1)
 
EXPERIMENT_DIR = Path(sys.argv[1]).resolve()
ROOT_PATH      = EXPERIMENT_DIR.parent          # e.g. …/Cellpaint/
FILE_NAME      = EXPERIMENT_DIR.name            # e.g. 20260319
 
SQLITE_FILE      = EXPERIMENT_DIR / "cell_profiler_out" / "CellPaint.db"
OUTPUT_CHUNK_DIR = EXPERIMENT_DIR / "output_chunks"
OUTPUT_PLT_DIR   = EXPERIMENT_DIR / "plt_analysis"
OUTPUT_RSLT_DIR  = EXPERIMENT_DIR / "rslt"
FINAL_OUTPUT     = OUTPUT_RSLT_DIR / "final_profiles.csv"
FILTERED_OUTPUT  = OUTPUT_RSLT_DIR / "filtered_profiles.csv"
METADATA_FILE    = Path(__file__).resolve().parent / "00_Metadata_Experiments.xlsx"
 
# SELECT_NAMES: derived from the folder name — everything up to the first
# purely-numeric segment is kept, then rejoined.
# e.g. "20260319_HDMVEC" → "20260319_HDMVEC"  (folder name is already the Name)
SELECT_NAMES = FILE_NAME if not FILE_NAME.startswith("20") else FILE_NAME[2:]   # override here if your folder name differs from the Name column
 
CHUNKSIZE  = 100_000
AGG_METHOD = "median"   # "median" or "mean"
BOOL_NEW   = sys.argv[2] == "True" # True = recompute from SQLite even if CSV already exists
 
#only create output directories if the SQLITE file exists, otherwise we might be pointing to the wrong experiment folder
if not SQLITE_FILE.exists():
    raise FileNotFoundError(f"SQLite file not found at {SQLITE_FILE}. Please check the Path you provided: {EXPERIMENT_DIR}.")

os.makedirs(OUTPUT_CHUNK_DIR, exist_ok=True)
os.makedirs(OUTPUT_PLT_DIR,   exist_ok=True)
os.makedirs(OUTPUT_RSLT_DIR,  exist_ok=True)
 
print(f"Experiment : {EXPERIMENT_DIR}")
print(f"Root       : {ROOT_PATH}")
print(f"SQLite     : {SQLITE_FILE}")
print(f"Metadata   : {METADATA_FILE}")
print(f"Plots  →   : {OUTPUT_PLT_DIR}")
print(f"Results →  : {OUTPUT_RSLT_DIR}")


# --- Metadata & experiment --- Must match "Name" column in Excel
MAX_IMAGES_WARNING   = 864                 # Warn if Per_Image rows exceed this
CELL_COUNT_THRESH    = 2.0                 # Drop images < mean - N*std cell count
CELL_AREA_THRESH     = 2.0                 # Drop images > mean + N*std cell area
QC_FOCUS_THRESH      = 0#.5                 # Drop images < median * N (focus score)

# --- Data loading & aggregation ---
DROP_KEYWORDS        = ["COSTES", "PARENT", "CHILD"]  # Feature name exclusions
DROP_OBJECT_NUMBER   = True
DROP_WELL            = True

# --- Normalization ---
NORM_METHOD          = "standardize"     # "mad_robustize", "standardize", "robustize"
MAD_EPS              = 0                  # Epsilon for MAD robustize

# --- Feature selection ---
FEATURE_OPS          = ["variance_threshold", "correlation_threshold", "drop_na_columns"]
CORR_THRESHOLD       = 0.8
OUTLIER_CUTOFF       = 10

# --- Heatmap & treatment clustering ---
SCALER               = "StandardScaler"   # "StandardScaler" or "RobustScaler"
TOP_FEATURES         = 60
HEATMAP_CMAP         = "PiYG"
N_CLUSTERS           = 8
CLUST_METRIC         = "cosine"
CLUST_LINKAGE        = "average"

# --- Dimensionality reduction ---
N_SAMPLE             = 20_000
PCA_N_COMPONENTS     = 60
TSNE_PERPLEXITY      = 100
TSNE_MAX_ITER        = 1000
TSNE_INIT            = "pca"
TSNE_LR              = "auto"
UMAP_N_NEIGHBORS     = 10
UMAP_MIN_DIST        = 0.1
UMAP_METRIC          = "cosine"
RANDOM_SEED          = 42

DO_PCA               = True
DO_TSNE              = True
DO_SC_TSNE           = True
DO_UMAP              = False

# --- Plot output ---
PLT_DPI              = 450
TRANSPARENT          = True
SCATTER_SIZE_SC      = 5
SCATTER_SIZE         = 20
LEGEND_FONTSIZE      = 10
LEGEND_MARKERSCALE   = 3
ANNOT_FONTSIZE       = 5

COLOR_BY_TREATMENT         = False
PLOT_WELLMAP               = True
PLOT_CLUSTERMAP_V          = True
PLOT_CLUSTERMAP_H          = True
PLOT_SIMILARITY            = True
PLOT_SIMILARITY_CLUSTERED  = False
PLOT_TSNE                  = True
PLOT_TSNE_CLUSTERED        = True
PLOT_TSNE_GRID             = True
PLOT_UMAP                  = True


# =============================================================================
# METADATA
# =============================================================================
print("\n[1/8] Loading metadata...")
metadata = pd.read_excel(METADATA_FILE)
metadata = metadata.dropna(
    subset=["Well_Row", "Well_Column", "Treatment", "Incub time [h]", "Date"]
)
metadata["Image_Metadata_Well"] = (
    metadata["Well_Row"].astype(str) + metadata["Well_Column"].astype(int).astype(str)
)
metadata["treatment_full"] = (
    metadata["Treatment"].astype(str).str.strip() + "_" +
    metadata["Concentration (ng)"].astype(str).str.strip() + "_" +
    metadata["Incub time [h]"].astype(str).str.strip()
)
# remove space from treatment_full

# Parse date from FILE_NAME
date_str = [s for s in str(FILE_NAME).split("_") if s.isdigit()]
if not date_str:
    raise ValueError(f"No numeric date found in FILE_NAME='{FILE_NAME}'. Expected e.g. '20260319' or '20260319_HAEC'.")
if date_str[0].startswith("20"):
    date = int(date_str[0][2:])  # last 6 digits → YYMMDD
else:
    date = int(date_str[0])       # already YYMMDD
metadata = metadata[metadata["Date"] == date]
# metadata = metadata.loc[metadata["Concentration (ng)"].isna()].fillna(0).copy()
metadata["concentration_numeric"] = metadata["Concentration (ng)"].apply(parse_value)
if COLOR_BY_TREATMENT:
    comp_col = "Treatment"
else:
    comp_col = "treatment_full"
metadata = compute_alpha_per_group(metadata, comp_col=comp_col, conc_col="concentration_numeric").copy()
metadata, color_dict = assign_colors(metadata, comp_col=comp_col, conc_col="concentration_numeric")
metadata = metadata.loc[:, ~metadata.columns.str.startswith("Unnamed")]
# remove "20" prefix from Name if present, to match SELECT_NAMES
metadata["Name"] = metadata["Name"].apply(lambda x: x[2:] if isinstance(x, str) and x.startswith("20") else x)
metadata = metadata[metadata["Name"] == SELECT_NAMES]

# Plot Wellmap
if PLOT_WELLMAP:
    plot_well_map(metadata,save_path=Path(OUTPUT_PLT_DIR) / 'PlateMap.png')


if metadata.empty:
    raise ValueError(f"No metadata rows match SELECT_NAMES='{SELECT_NAMES}' for date {date}.")
print(metadata["Treatment"].unique())
print(f"  {len(metadata)} metadata rows loaded.")
print(metadata.head())

ctrl_rows = metadata[metadata.Treatment.str.contains(r"control|ctrl|dmso|dms0", regex=True, case=False, na=False)]
if ctrl_rows.empty:
    raise ValueError("No control treatment found in metadata (looked for 'control' in Treatment column).")
print(f"  Control treatments found: {ctrl_rows['treatment_full'].tolist()}")
CONTROL_LABEL = ctrl_rows["treatment_full"].values[0]
for arg in sys.argv:
    if arg.startswith('CONTROL_LABE='):
        CONTROL_LABEL = arg.split('=')[-1]
        
print(f"  Control label: {CONTROL_LABEL}")

metadata_cols = [
    "Image_Metadata_Well", "treatment_full", "Treatment",
    "Concentration (ng)", "Incub time [h]"
]

# =============================================================================
# CONNECT TO SQLITE
# =============================================================================
print("\n[2/8] Connecting to SQLite...")
conn = sqlite3.connect(SQLITE_FILE)
table_names = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table';", conn)
print(f"  Tables: {table_names['name'].tolist()}")

# =============================================================================
# PER-IMAGE METADATA + OUTLIER DETECTION
# =============================================================================
print("\n[3/8] Loading Per_Image metadata and detecting outliers...")
per_image = pd.read_sql(
    "SELECT ImageNumber, Image_Metadata_Well, Image_Count_Cells, Mean_Cells_AreaShape_Area , Image_ImageQuality_FocusScore_OrigDNA "
    "FROM Per_Image",
    conn
)
print(per_image.head())
if per_image.shape[0] > MAX_IMAGES_WARNING:
    print(f"  WARNING: Per_Image has {per_image.shape[0]} rows (> {MAX_IMAGES_WARNING}). "
          "Truncating to first {MAX_IMAGES_WARNING}.")
    per_image = per_image.head(MAX_IMAGES_WARNING)

per_image = per_image.merge(metadata, on="Image_Metadata_Well", how="left")
cell_count = per_image["Image_Count_Cells"]
cell_area  = per_image["Mean_Cells_AreaShape_Area"]
qc_focus   = per_image['Image_ImageQuality_FocusScore_OrigDNA']

plot_quality_control(OUTPUT_PLT_DIR,cell_area,cell_count,qc_focus,CELL_AREA_THRESH,CELL_COUNT_THRESH, QC_FOCUS_THRESH)

drop_images       = (
    (cell_count < cell_count.mean() - CELL_COUNT_THRESH * cell_count.std()) |
    (cell_area  > cell_area.mean()  + CELL_AREA_THRESH  * cell_area.std()) |
    (qc_focus   < qc_focus.median() * QC_FOCUS_THRESH)
)

drop_image_numbers = per_image["ImageNumber"][drop_images]
print(f"  Dropping {drop_images.sum()} outlier images out of {len(per_image)}.")

# =============================================================================
# STREAM + AGGREGATE PER_OBJECT
# =============================================================================
print("\n[4/8] Streaming Per_Object and aggregating...")
if not FINAL_OUTPUT.exists() or BOOL_NEW:
    if BOOL_NEW:
        for f in os.listdir(OUTPUT_CHUNK_DIR):
            if f.startswith("chunk"):
                os.remove(OUTPUT_CHUNK_DIR / f)

    chunk_iter = pd.read_sql("SELECT * FROM Per_Object", conn, chunksize=CHUNKSIZE)
    chunk_files = []

    for i, chunk in tqdm(enumerate(chunk_iter)):
        chunk = chunk.merge(per_image, on="ImageNumber", how="left")
        chunk = chunk[~chunk["ImageNumber"].isin(drop_image_numbers)]

        numeric_cols = chunk.select_dtypes(include="number").columns.tolist()
        numeric_cols = [c for c in numeric_cols if c != "ImageNumber"]
        if DROP_OBJECT_NUMBER:
            numeric_cols = [c for c in numeric_cols if "ObjectNumber" not in c]
        if DROP_WELL:
            numeric_cols = [c for c in numeric_cols if "Well" not in c]

        if AGG_METHOD == "median":
            agg = chunk.groupby(["ImageNumber", "treatment_full"])[numeric_cols].median()
        else:
            agg = chunk.groupby(["ImageNumber", "treatment_full"])[numeric_cols].mean()

        chunk_file = OUTPUT_CHUNK_DIR / f"chunk_{i}.parquet"
        agg.to_parquet(chunk_file)
        chunk_files.append(chunk_file)
        del chunk, agg

    print("  Combining chunks...")
    dfs     = [pd.read_parquet(f) for f in chunk_files]
    combined = pd.concat(dfs)

    if AGG_METHOD == "median":
        final_profiles = combined.groupby(level=[0, 1]).median().reset_index()
    else:
        final_profiles = combined.groupby(level=[0, 1]).mean().reset_index()
    del dfs, combined

   
    if "Name" not in final_profiles.columns:
        final_profiles = final_profiles.merge(
            metadata[["Name", "treatment_full"]].drop_duplicates(),
            on="treatment_full", how="left"
        )
    print(final_profiles)
    final_profiles = final_profiles.merge(
    per_image.loc[
        ~per_image.ImageNumber.isin(drop_image_numbers),
        ["ImageNumber", "Image_Metadata_Well"]
    ],
    on="ImageNumber",
    how="left"
    )
    keep_keys = [c for c in final_profiles.columns if "unnamed" not in c.lower()]
    final_profiles[keep_keys].to_csv(FINAL_OUTPUT, index=False)
    print(f"  Saved → {FINAL_OUTPUT}")
else:
    print(f"  Loading existing {FINAL_OUTPUT}")
    final_profiles = pd.read_csv(FINAL_OUTPUT)
    keep_keys = [c for c in final_profiles.columns if "unnamed" not in c.lower()]
    final_profiles = final_profiles[keep_keys]

# =============================================================================
# NORMALIZE
# =============================================================================
print("\n[5/8] Normalizing to control...")
meta_cols    = ["ImageNumber", "Cells_Number_Object_Number", "Name", "treatment_full",
                "Date", "concentration_numeric", "alpha","Image_Metadata_Well"]
feature_cols = [c for c in final_profiles.columns if c not in meta_cols]

normalized = normalize(
    profiles=final_profiles,
    features=feature_cols,
    meta_features=[c for c in meta_cols if c in final_profiles.columns],
    method=NORM_METHOD,
    mad_robustize_epsilon=MAD_EPS,
    samples=f"treatment_full == '{CONTROL_LABEL}'"
)

# =============================================================================
# FEATURE SELECTION
# =============================================================================
print("\n[6/8] Feature selection...")
feature_cols = [c for c in final_profiles.columns if c not in meta_cols]

selected = feature_select(
    profiles=normalized,
    features=feature_cols,
    operation=FEATURE_OPS,
    corr_threshold=CORR_THRESHOLD,
    outlier_cutoff=OUTLIER_CUTOFF,
)

selected = selected[[
    c for c in selected.columns
    if all(kw not in c.upper() for kw in DROP_KEYWORDS)
]]
print(f"  {selected.shape[1]} features retained, {selected.shape[0]} images.")
selected.to_csv(FILTERED_OUTPUT, index=False)
print(f"  Saved → {FILTERED_OUTPUT}")

# =============================================================================
# CONDITION-LEVEL AGGREGATION + SCALING
# =============================================================================
print("\n[7/8] Condition-level aggregation and scaling...")
meta_cols_ext = meta_cols + ["Nuclei_Number_Object_Number"]
feature_cols  = [c for c in selected.columns
                 if c not in meta_cols_ext and c != "Mean_Cells_AreaShape_Area"]

condition_profiles = selected.groupby("treatment_full")[feature_cols].median()

if SCALER == "RobustScaler":
    scaler = RobustScaler()
else:
    scaler = StandardScaler()

scaled_data = scaler.fit_transform(condition_profiles)
scaled_df   = pd.DataFrame(scaled_data, index=condition_profiles.index,
                            columns=condition_profiles.columns)

variances    = condition_profiles.var(axis=0)
top_features = variances.sort_values(ascending=False).head(TOP_FEATURES).index
heatmap_data = scaled_df[top_features]

# =============================================================================
# PLOTS — HEATMAPS
# =============================================================================
print("\n[8/8] Generating plots...")

if PLOT_CLUSTERMAP_V:
    g = sns.clustermap(
        heatmap_data.T, cmap=HEATMAP_CMAP, center=0,
        cbar_pos=(0.6, 0.999, 0.3, 0.01),
        xticklabels=heatmap_data.index,
        yticklabels=heatmap_data.columns,
        figsize=(9, 12), dendrogram_ratio=(0.1, 0.03),
        cbar_kws={"label": "Z-score", "orientation": "horizontal"},
    )
    g.ax_heatmap.tick_params(axis="x", labelsize=6)
    g.ax_heatmap.tick_params(axis="y", labelsize=6)
    g.ax_heatmap.set_xlabel(""); g.ax_heatmap.set_ylabel("")
    plt.savefig(OUTPUT_PLT_DIR / "clustermap.png", dpi=PLT_DPI,
                bbox_inches="tight", transparent=TRANSPARENT)
    plt.close()
    print("  Saved clustermap.png")

if PLOT_CLUSTERMAP_H:
    g = sns.clustermap(
        heatmap_data, cmap=HEATMAP_CMAP, center=0,
        cbar_pos=(0.85, 0.2, 0.1, 0.05),
        yticklabels=heatmap_data.index,
        xticklabels=heatmap_data.columns,
        figsize=(12, 8), dendrogram_ratio=(0.03, 0.1),
        cbar_kws={"label": "Z-score", "orientation": "horizontal"},
    )
    g.ax_heatmap.tick_params(axis="x", labelsize=8)
    g.ax_heatmap.tick_params(axis="y", labelsize=8)
    g.ax_heatmap.set_xlabel(""); g.ax_heatmap.set_ylabel("")
    plt.savefig(OUTPUT_PLT_DIR / "clustermap_horizontal.png", dpi=PLT_DPI,
                bbox_inches="tight", transparent=TRANSPARENT)
    plt.close()
    print("  Saved clustermap_horizontal.png")

# --- Treatment similarity ---
similarity    = cosine_similarity(condition_profiles[top_features])
similarity_df = pd.DataFrame(similarity, index=condition_profiles.index,
                              columns=condition_profiles.index)

if PLOT_SIMILARITY:
    plt.figure(figsize=(6, 5))
    g = sns.heatmap(similarity_df, cmap="viridis", xticklabels=True, yticklabels=True,
                    annot=True, linewidth=0.5, annot_kws={"size": ANNOT_FONTSIZE})
    g.set_xticklabels(g.get_xticklabels(), rotation=90, fontsize=6)
    g.set_yticklabels(g.get_yticklabels(), rotation=0,  fontsize=6)
    g.set_xlabel(""); g.set_ylabel("")
    g.xaxis.tick_top()
    plt.title("Treatment Similarity (Cosine)")
    plt.savefig(OUTPUT_PLT_DIR / "treatment_similarity.png", dpi=PLT_DPI,
                bbox_inches="tight", transparent=TRANSPARENT)
    plt.close()
    print("  Saved treatment_similarity.png")

# --- Clustered similarity ---
clustering = AgglomerativeClustering(
    n_clusters=N_CLUSTERS, linkage=CLUST_LINKAGE, metric=CLUST_METRIC
)
clusters = clustering.fit_predict(condition_profiles[top_features])
similarity_df["Cluster"] = clusters
similarity_df = similarity_df.sort_values("Cluster")
treatment_cluster_mapping = dict(zip(similarity_df.index, similarity_df["Cluster"]))
similarity_df = similarity_df.drop(columns=["Cluster"])

if PLOT_SIMILARITY_CLUSTERED:
    plt.figure(figsize=(6, 5))
    g = sns.heatmap(similarity_df, cmap="viridis", xticklabels=True, yticklabels=True,
                    annot=True, linewidth=0.5, annot_kws={"size": ANNOT_FONTSIZE})
    g.set_xticklabels(g.get_xticklabels(), rotation=90, fontsize=6)
    g.set_yticklabels(g.get_yticklabels(), rotation=0,  fontsize=6)
    g.set_xlabel(""); g.set_ylabel("")
    g.xaxis.tick_top()
    plt.title("Treatment Similarity (Cosine) — Clustered")
    plt.savefig(OUTPUT_PLT_DIR / "treatment_similarity_clustered.png", dpi=PLT_DPI,
                bbox_inches="tight", transparent=TRANSPARENT)
    plt.close()
    print("  Saved treatment_similarity_clustered.png")

# =============================================================================
# PER-CELL DIM REDUCTION
# =============================================================================
if DO_PCA or DO_TSNE or DO_UMAP:
    print("\n  Sampling per-cell data from SQLite...")
    meta_cols_sample = ["ImageNumber", "Image_Count_Cells", "Cells_Number_Object_Number", 'Image_ImageQuality_FocusScore_OrigDNA',
                        "Name", "treatment_full", "Date", "concentration_numeric", "alpha","Image_Metadata_Well"]
    feature_cols_sample = [c for c in selected.columns
                           if c not in meta_cols_sample and c != "Mean_Cells_AreaShape_Area"]
    all_cols = ["ImageNumber"] + feature_cols_sample

    n_total = pd.read_sql_query("SELECT COUNT(*) as cnt FROM Per_Object", conn).iloc[0, 0]
    print(f"  Total Per_Object rows: {n_total}")

    placeholders = ",".join(["?"] * len(drop_image_numbers.values))
    query = f"""
        SELECT {', '.join(all_cols)}
        FROM Per_Object
        WHERE ImageNumber NOT IN ({placeholders})
        AND ImageNumber < {MAX_IMAGES_WARNING}
        ORDER BY RANDOM()
        LIMIT {N_SAMPLE}
    """
    df_sample = pd.read_sql_query(query, conn, params=drop_image_numbers.values)
    df_sample = df_sample.merge(final_profiles[["ImageNumber", "treatment_full"]],
                                on="ImageNumber", how="left")
    df_sample = df_sample.merge(per_image[["ImageNumber", "Image_Metadata_Well"]],
                                on="ImageNumber", how="left")

    X = df_sample[feature_cols_sample].copy()
    y = df_sample[["treatment_full"]].copy()
    y = y.merge(
        metadata[["treatment_full", "Name", "concentration_numeric",
                  "plot_color", "alpha"]].drop_duplicates(),
        on="treatment_full", how="left"
    )
    X = X.replace([np.inf, -np.inf], np.nan).dropna()
    X = X.astype(np.float32)
    y = y.loc[X.index]

    cell_scaler = StandardScaler()
    X_scaled    = cell_scaler.fit_transform(X)

# --- PCA ---
if DO_PCA:
    print("  Running PCA...")
    pca   = PCA(n_components=PCA_N_COMPONENTS, random_state=RANDOM_SEED)
    X_pca = pca.fit_transform(X_scaled)

    cmap_pca = plt.cm.magma
    norm_pca  = plt.Normalize(
        vmin=np.percentile(X_pca[:, 2], 0.1),
        vmax=np.percentile(X_pca[:, 2], 99.9)
    )
    fig, ax = plt.subplots(ncols=2, figsize=(8, 4))
    ax[1].scatter(X_pca[:, 0], X_pca[:, 1], c=cmap_pca(norm_pca(X_pca[:, 2])), alpha=0.2)
    ax[1].set_xlim(np.percentile(X_pca[:, 0], 0.01), np.percentile(X_pca[:, 0], 99.99))
    ax[1].set_ylim(np.percentile(X_pca[:, 1], 0.01), np.percentile(X_pca[:, 1], 99.99))
    ax[0].scatter(np.arange(len(pca.explained_variance_)), pca.explained_variance_)
    ax[0].set_xlabel("PCA component"); ax[0].set_ylabel("Explained variance")
    plt.tight_layout()
    plt.savefig(OUTPUT_PLT_DIR / "pca.png", dpi=PLT_DPI,
                bbox_inches="tight", transparent=TRANSPARENT)
    plt.close()
    print("  Saved pca.png")

# --- t-SNE ---
if DO_SC_TSNE:
    print(f"  Running Single cell t-SNE (perplexity={TSNE_PERPLEXITY}, max_iter={TSNE_MAX_ITER})...")
    tsne = TSNE(
        n_components=2, perplexity=TSNE_PERPLEXITY,
        learning_rate=TSNE_LR, init=TSNE_INIT,
        random_state=RANDOM_SEED, max_iter=TSNE_MAX_ITER
    )
    X_tsne = tsne.fit_transform(X_scaled)

    treatments = metadata["treatment_full"].unique()

    if PLOT_TSNE:
        plt.figure(figsize=(10, 8))
        for t in treatments:
            color = metadata[metadata["treatment_full"] == t]["plot_color"].unique()[0]
            for conc in metadata[metadata["treatment_full"] == t]["concentration_numeric"].unique():
                mask = (y["treatment_full"] == t) & (y["concentration_numeric"] == conc)
                if t == CONTROL_LABEL:
                    clr, alpha = "k", 0.5
                else:
                    clr   = y[mask]["plot_color"].unique()[0] if mask.any() else color
                    alpha = y[mask]["alpha"].unique()[0]      if mask.any() else 0.5
                plt.scatter(X_tsne[mask, 0], X_tsne[mask, 1],
                            s=SCATTER_SIZE_SC, label=t, color=clr, alpha=alpha)
        plt.xlabel("t-SNE 1"); plt.ylabel("t-SNE 2")
        plt.title("Per-cell morphology t-SNE")
        plt.legend(markerscale=LEGEND_MARKERSCALE, fontsize=LEGEND_FONTSIZE,
                   bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)
        plt.savefig(OUTPUT_PLT_DIR / "tsne_per_cell.png", dpi=PLT_DPI,
                    bbox_inches="tight", transparent=TRANSPARENT)
        plt.close()
        print("  Saved tsne_per_cell.png")

    if PLOT_TSNE_CLUSTERED:
        clusters_unique = np.unique(list(treatment_cluster_mapping.values()))
        cmap_c = plt.cm.tab20
        norm_c = plt.Normalize(vmin=clusters_unique.min(), vmax=clusters_unique.max())
        plt.figure(figsize=(10, 8))
        for t in treatment_cluster_mapping:
            mask = y["treatment_full"] == t
            plt.scatter(X_tsne[mask, 0], X_tsne[mask, 1],
                        s=SCATTER_SIZE_SC, label=t,
                        color=cmap_c(norm_c(treatment_cluster_mapping[t])))
        plt.xlabel("t-SNE 1"); plt.ylabel("t-SNE 2")
        plt.title("Per-cell morphology t-SNE — by cluster")
        plt.legend(markerscale=LEGEND_MARKERSCALE, fontsize=LEGEND_FONTSIZE,
                   bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)
        plt.savefig(OUTPUT_PLT_DIR / "tsne_per_cell_clustered.png", dpi=PLT_DPI,
                    bbox_inches="tight", transparent=TRANSPARENT)
        plt.close()
        print("  Saved tsne_per_cell_clustered.png")

    if PLOT_TSNE_GRID:
        n_side = int(np.ceil(np.sqrt(len(treatments))))
        fig, ax = plt.subplots(n_side, n_side, figsize=(10, 10))
        for i, t in enumerate(treatments):
            for conc in metadata[metadata["treatment_full"] == t]["concentration_numeric"].unique():
                mask = (y["treatment_full"] == t) & (y["concentration_numeric"] == conc)
                ax.ravel()[i].scatter(
                    X_tsne[mask, 0], X_tsne[mask, 1],
                    s=SCATTER_SIZE_SC, label=t,
                    color=y[mask]["plot_color"].unique()[0] if mask.any() else "gray",
                    alpha=y[mask]["alpha"].unique()[0]      if mask.any() else 0.5
                )
                ax.ravel()[i].axis("off")
                ax.ravel()[i].set_title(t, fontsize=ANNOT_FONTSIZE)
                ax.ravel()[i].scatter(X_tsne[~mask, 0], X_tsne[~mask, 1],
                                      s=SCATTER_SIZE_SC, color="#AAAAAA", alpha=0.1, zorder=-1)
        plt.suptitle("Per-cell t-SNE — small multiples", fontsize=10)
        plt.tight_layout()
        plt.savefig(OUTPUT_PLT_DIR / "tsne_grid.png", dpi=PLT_DPI,
                    bbox_inches="tight", transparent=TRANSPARENT)
        plt.close()
        print("  Saved tsne_grid.png")

if DO_TSNE:
    print(f"  Running image-level t-SNE (perplexity={TSNE_PERPLEXITY}, max_iter={TSNE_MAX_ITER})...")

    # --- Use final_profiles: one row per ImageNumber ---
    # Re-use the same feature cols as condition-level, but from final_profiles
    meta_cols_ext = meta_cols + ["Nuclei_Number_Object_Number"]
    feature_cols_img = [c for c in selected.columns
                        if c not in meta_cols_ext and c != "Mean_Cells_AreaShape_Area"]

    # Align: only keep images that survived feature selection
    img_profiles = selected.copy()
    img_profiles = img_profiles.merge(
        metadata[["treatment_full", "plot_color", "alpha", "concentration_numeric"]].drop_duplicates(),
        on="treatment_full", how="left", suffixes=('','_drop')
    ).filter(regex='^(?!.*_drop$)')  # remove duplicate columns from merge

    X_img = img_profiles[feature_cols_img].replace([np.inf, -np.inf], np.nan).dropna()
    img_profiles = img_profiles.loc[X_img.index]

    img_scaler = StandardScaler()
    X_img_scaled = img_scaler.fit_transform(X_img)

    perp = min(TSNE_PERPLEXITY, len(X_img_scaled) - 1)
    tsne_img = TSNE(
        n_components=2, perplexity=perp,
        learning_rate=TSNE_LR, init=TSNE_INIT,
        random_state=RANDOM_SEED, max_iter=TSNE_MAX_ITER
    )
    X_tsne_img = tsne_img.fit_transform(X_img_scaled)

    img_profiles = img_profiles.copy()
    img_profiles["tSNE1"] = X_tsne_img[:, 0]
    img_profiles["tSNE2"] = X_tsne_img[:, 1]
    treatments = img_profiles["treatment_full"].unique()

    if PLOT_TSNE:
        fig, ax = plt.subplots(figsize=(10, 8))
        for t in treatments:
            mask = img_profiles["treatment_full"] == t
            rows = img_profiles[mask]
            color = "k" if t == CONTROL_LABEL else rows["plot_color"].iloc[0]
            alpha = 0.5 if t == CONTROL_LABEL else rows["alpha"].iloc[0]
            ax.scatter(rows["tSNE1"], rows["tSNE2"],
                       s=SCATTER_SIZE, color=color, alpha=alpha, label=t, zorder=2)
        ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
        ax.set_title("Per-image morphology t-SNE")
        ax.legend(markerscale=LEGEND_MARKERSCALE, fontsize=LEGEND_FONTSIZE,
                  bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(OUTPUT_PLT_DIR / "tsne_per_image.png", dpi=PLT_DPI,
                    bbox_inches="tight", transparent=TRANSPARENT)
        plt.close()
        print("  Saved tsne_per_image.png")

    if PLOT_TSNE_CLUSTERED:
        clusters_unique = np.unique(list(treatment_cluster_mapping.values()))
        cmap_c = plt.cm.tab20
        norm_c = plt.Normalize(vmin=clusters_unique.min(), vmax=clusters_unique.max())
        fig, ax = plt.subplots(figsize=(10, 8))
        for t in treatments:
            mask = img_profiles["treatment_full"] == t
            rows = img_profiles[mask]
            cluster_id = treatment_cluster_mapping.get(t, 0)
            ax.scatter(rows["tSNE1"], rows["tSNE2"],
                       s=SCATTER_SIZE, label=t,
                       color=cmap_c(norm_c(cluster_id)), zorder=2)
        ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
        ax.set_title("Per-image t-SNE — by cluster")
        ax.legend(markerscale=LEGEND_MARKERSCALE, fontsize=LEGEND_FONTSIZE,
                  bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(OUTPUT_PLT_DIR / "tsne_per_image_clustered.png", dpi=PLT_DPI,
                    bbox_inches="tight", transparent=TRANSPARENT)
        plt.close()
        print("  Saved tsne_per_image_clustered.png")

    if PLOT_TSNE_GRID:
        n_side = int(np.ceil(np.sqrt(len(treatments))))
        fig, axes = plt.subplots(n_side, n_side, figsize=(10, 10))
        for i, t in enumerate(treatments):
            ax_i = axes.ravel()[i]
            mask = img_profiles["treatment_full"] == t
            rows = img_profiles[mask]
            ax_i.scatter(img_profiles["tSNE1"], img_profiles["tSNE2"],
                         s=SCATTER_SIZE, color="#AAAAAA", alpha=0.1, zorder=1)
            color = "k" if t == CONTROL_LABEL else rows["plot_color"].iloc[0]
            alpha = 0.5 if t == CONTROL_LABEL else rows["alpha"].iloc[0]
            ax_i.scatter(rows["tSNE1"], rows["tSNE2"],
                         s=SCATTER_SIZE, color=color, alpha=alpha, zorder=2)
            ax_i.set_title(t, fontsize=ANNOT_FONTSIZE)
            ax_i.axis("off")
        for j in range(len(treatments), n_side ** 2):
            axes.ravel()[j].axis("off")
        plt.suptitle("Per-image t-SNE — small multiples", fontsize=10)
        plt.tight_layout()
        plt.savefig(OUTPUT_PLT_DIR / "tsne_per_image_grid.png", dpi=PLT_DPI,
                    bbox_inches="tight", transparent=TRANSPARENT)
        plt.close()
        print("  Saved tsne_per_image_grid.png")


# --- UMAP ---
if DO_UMAP:
    print("  Running UMAP (condition-level)...")
    import umap
    reducer = umap.UMAP(
        n_neighbors=UMAP_N_NEIGHBORS, min_dist=UMAP_MIN_DIST,
        metric=UMAP_METRIC, random_state=RANDOM_SEED
    )
    embedding = reducer.fit_transform(condition_profiles)
    umap_df   = pd.DataFrame(embedding, columns=["UMAP1", "UMAP2"],
                              index=condition_profiles.index)

    if PLOT_UMAP:
        plt.figure(figsize=(8, 6))
        plt.scatter(umap_df["UMAP1"], umap_df["UMAP2"])
        for i, label in enumerate(umap_df.index):
            plt.text(umap_df.iloc[i, 0], umap_df.iloc[i, 1], label, fontsize=8)
        plt.title("UMAP of treatments")
        plt.xlabel("UMAP1"); plt.ylabel("UMAP2")
        plt.savefig(OUTPUT_PLT_DIR / "umap_conditions.png", dpi=PLT_DPI,
                    bbox_inches="tight", transparent=TRANSPARENT)
        plt.close()
        print("  Saved umap_conditions.png")

conn.close()
print("\n✅ Pipeline complete.")
print(f"   Plots  → {OUTPUT_PLT_DIR}")
print(f"   Results → {OUTPUT_RSLT_DIR}")
