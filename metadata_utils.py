import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def parse_value(x):
    # Case 1: already numeric
    if isinstance(x, (int, float, np.number)):
        return float(x)
    
    if isinstance(x, str):
        x = x.strip()
        
        # Replace comma decimal with dot (European format)
        # BUT avoid messing up lists like "1,000" (optional improvement later)
        x = re.sub(r'(\d),(\d)', r'\1.\2', x)
        
        # Case 2: expression like "50 + 100 + 10"
        if "+" in x:
            try:
                parts = re.split(r"\+", x)
                numbers = []
                for p in parts:
                    match = re.search(r"[-+]?\d*\.?\d+", p)
                    if match:
                        numbers.append(float(match.group()))
                return sum(numbers) if numbers else np.nan
            except:
                pass
        
        # Case 3: value with unit like "1 µM", "0.5 mM", "0,3 µM"
        match = re.search(r"[-+]?\d*\.?\d+", x)
        if match:
            return float(match.group())
    
    return np.nan

def compute_alpha(conc, min_c, max_c):
    if max_c == min_c:
        return 1.0  # avoid division by zero
    
    norm = (conc - min_c) / (max_c - min_c)
    return 0.5 + 0.5 * norm

def get_color_map(compositions,cmap_name="tab20"):
    unique = sorted(compositions.unique())
    cmap = plt.get_cmap(cmap_name)  # good categorical palette
    
    color_dict = {
        comp: cmap(i % 20)
        for i, comp in enumerate(unique)
    }
    return color_dict

def compute_alpha_per_group(df, comp_col, conc_col):
    df = df.copy()
    
    df["alpha"] = df.groupby(comp_col)[conc_col].transform(
        lambda x: 0.5 + 0.5 * (x - x.min()) / (x.max() - x.min() + 1e-9)
    )
    
    return df

def assign_colors(df, comp_col="composition", conc_col="concentration", cmap_name="tab20b"):
    df = df.copy()
    # Get color per composition
    color_dict = get_color_map(df[comp_col], cmap_name=cmap_name)
    
    
    colors = []
    
    for _, row in df.iterrows():
        base_color = color_dict[row[comp_col]]
        colors.append(base_color)
    
    df["plot_color"] = colors
    return df, color_dict

def plot_well_map(metadata, well_col="Image_Metadata_Well", treatment_col="treatment_full",save_path=None):
    
    # extract row and column from well names (A01 → A, 1)
    rows = metadata[well_col].str[0]
    cols = metadata[well_col].str[1:].astype(int)
    
    plot_df = pd.DataFrame({
        "row": rows,
        "col": cols,
        "treatment": metadata[treatment_col],
        "color": metadata["plot_color"],
        "alpha": metadata["alpha"]  
    })
    
    # convert rows to numeric for plotting
    row_order = sorted(plot_df["row"].unique())
    row_map = {r:i for i,r in enumerate(row_order)}
    plot_df["row_num"] = plot_df["row"].map(row_map)
    # color treatments
    treatments = plot_df["treatment"].unique()
    colors = plt.cm.tab20(range(len(treatments)))
    color_map = dict(zip(treatments, colors))

    plt.figure(figsize=(10,5))

    for _, f in plot_df.sort_values('row_num').groupby('treatment'):
        plt.scatter(
            f.col,
            f.row_num,
            s=250,
            color=f.color,
            edgecolor="black", label=f.treatment.iloc[0], alpha=f.alpha.iloc[0]
        )

    plt.yticks(range(len(row_order)), row_order)
    plt.xticks(np.arange(1,1+len(np.unique(cols))))
    plt.xlabel("Column")
    plt.ylabel("Row")
    plt.title("Plate well map")
    plt.gca().invert_yaxis()
    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left",fontsize=10, markerscale=.7, borderaxespad=0.2)

    plt.tight_layout()
    if save_path == None:
        plt.show()
    else:
        plt.savefig(save_path)

def plot_quality_control(save_dir, cell_area, cell_count, qc_focus, area_thresh, count_thresh, qc_focus_thresh):
    plt.figure(figsize=(15,5))
    
    plt.subplot(1,3,1)
    cell_area.hist(bins=100)
    plt.axvline(cell_area.mean(), color='k', linestyle='dashed', linewidth=1)
    plt.axvline(cell_area.mean() - area_thresh * cell_area.std(), color='g', linestyle='dashed', linewidth=1)
    plt.title("Cell Area Distribution")
    
    plt.subplot(1,3,2)
    cell_count.hist(bins=100)
    plt.axvline(cell_count.mean(), color='k', linestyle='dashed', linewidth=1)
    plt.axvline(cell_count.mean() - count_thresh * cell_count.std(), color='g', linestyle='dashed', linewidth=1)
    plt.title("Cell Count Distribution")
    
    plt.subplot(1,3,3)
    qc_focus.hist(bins=100)
    plt.axvline(qc_focus.median(), color='k', linestyle='dashed', linewidth=1)
    plt.axvline(qc_focus.median()*qc_focus_thresh, color='g', linestyle='dashed', linewidth=1)
    plt.title("Focus Score Distribution")
    
    plt.tight_layout()
    
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(os.path.join(save_dir, "quality_control.png"))
    else:
        plt.show()
