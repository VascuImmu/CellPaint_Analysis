import os
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from aicsimageio import AICSImage
from aicsimageio.writers import OmeTiffWriter
import dask.array as da
import numpy as np
import matplotlib.pyplot as plt
import cmap
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib
matplotlib.use("Agg")  # critical for multiprocessing


def five_color_plot(
    img,
    channel_names=None,
    merge=False,
    colored_subplots = True,
    title=None,
    pixel_size=None,
    scalebar_length=100,  # in microns
    save_path=None,
    dpi = 300
):
    if isinstance(img, da.Array):
        img = img.compute().squeeze()
    # perceptually balanced colormaps
    cmaps = [
        cmap.Colormap("gray"),     # grayscale
        cmap.Colormap("green"),    # green
        cmap.Colormap("cyan"),     # cyan
        cmap.Colormap("yellow"),   # yellow
        cmap.Colormap("magenta"),  # magenta
    ]

    n_channels = 5

    if channel_names is None:
        channel_names = [f"Ch{i}" for i in range(n_channels)]

    n_panels = n_channels + 1 if merge else n_channels
    fig, ax = plt.subplots(1, n_panels, figsize=(4 * n_panels, 4))

    if n_panels == 1:
        ax = [ax]

    # normalize channels for plotting
    norm_channels = []
    for i in range(n_channels):
        vmax = np.percentile(img[i], 99.5)
        norm = np.clip(img[i] / vmax, 0, 1)
        norm_channels.append(norm)
    del img  # free memory
    # -------------------
    # individual channels
    # -------------------
    for i in range(n_channels):

        cmap_use = cmaps[i] if colored_subplots else cmaps[0]

        ax[i].imshow(norm_channels[i], cmap=cmap_use.to_mpl())
        ax[i].set_title(channel_names[i], fontsize=12)
        ax[i].axis("off")

    # -------------------
    # merged image
    # -------------------
    if merge:

        rgb_stack = []

        for i in range(n_channels):

            # convert normalized channel -> RGBA
            rgba = cmaps[i].to_mpl()(norm_channels[i])

            # drop alpha channel
            rgb_channel = rgba[..., :3]

            rgb_stack.append(rgb_channel)

        # shape → (C, Y, X, 3)
        rgb_stack = np.stack(rgb_stack, axis=0)

        # napari-style additive blending
        rgb = 1 - np.prod(1 - rgb_stack, axis=0)

        rgb = np.clip(rgb, 0, 1)

        ax[-1].imshow(rgb)
        ax[-1].set_title("Merge", fontsize=12)
        ax[-1].axis("off")

    # -------------------
    # title
    # -------------------
    if title is not None:
        plt.suptitle(title, fontsize=16)

    plt.subplots_adjust(wspace=0.02)

    # -------------------
    # scalebar
    # -------------------
    if pixel_size is not None:

        scalebar_um = scalebar_length
        scalebar_px = scalebar_um / pixel_size
        for a in ax:
            scalebar = AnchoredSizeBar(
                a.transData,
                scalebar_px,
                "",
                "lower right",
                pad=0.1,
                borderpad=0.4,
                color="white",
                frameon=False,
                size_vertical=30,
                label_top=True
            )
            a.add_artist(scalebar)

    # -------------------
    # save or show
    # -------------------
    if save_path is not None:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close()
        print(f'✓ Saved plot to {save_path}')
    else:
        plt.show()

# -------------------------------------------------
# Worker function (must be top-level for Windows)
# -------------------------------------------------
def process_scene(
    lif_path_str,
    scene,
    output_dir_str,
    chunk_z,
    bool_imsave,
    plt_dir,
    overwrite=False
):
    lif_path = Path(lif_path_str)
    output_dir = Path(output_dir_str)
    print(scene)

    safe_scene = scene.replace("/", "_").replace(":", "_")
    out_path = output_dir / f"{lif_path.stem}_{safe_scene}.ome.tif"
    plt_path = plt_dir / f"{lif_path.stem}_{safe_scene}.png"

    is_p5 = scene.endswith('P5')

    # -------------------------
    # Decide what is needed
    # -------------------------
    need_tif = overwrite or not out_path.exists()
    need_plot = (
        bool_imsave
        and is_p5
        and (overwrite or not plt_path.exists())
    )

    # -------------------------
    # Skip entirely if nothing needed
    # -------------------------
    if not need_tif and not need_plot:
        return f"↷ Skipped (exists): {out_path.name}"

    try:
        img = AICSImage(lif_path)
        img.set_scene(scene)
        pixelsize = img.physical_pixel_sizes[1]

        # -------------------------
        # Load data ONLY if needed
        # -------------------------
        if chunk_z:
            data = img.get_image_dask_data("TCZYX", chunks={"Z": 1})
        else:
            data = img.get_image_dask_data("TCZYX")

        # -------------------------
        # Save TIF
        # -------------------------
        if need_tif:
            OmeTiffWriter.save(
                data,
                uri=str(out_path),
                dim_order="TCZYX"
            )

        # -------------------------
        # Save plot (P5 only)
        # -------------------------
        if need_plot:
            five_color_plot(
                data,
                channel_names=['DNA','ER','RNA','AGP','Mito'],
                pixel_size=pixelsize,
                save_path=plt_path
            )

        # -------------------------
        # Return status
        # -------------------------
        if need_tif and need_plot:
            del data
            del img
            return f"✓ Saved TIF + plot: {out_path.name}"
        elif need_tif:
            del data
            del img
            return f"✓ Saved TIF: {out_path.name}"
        elif need_plot:
            del data
            del img
            return f"✓ Saved plot: {plt_path.name}"
        

    except Exception as e:
        return f"✗ Error processing {scene}: {e}"

# -------------------------------------------------
# Process one LIF file (parallel over scenes)
# -------------------------------------------------
def process_lif_parallel(lif_path, output_dir, n_workers, chunk_z,bool_imsave=False,plt_path='',overwrite=False):
    print(f"\nProcessing {lif_path.name}")

    img = AICSImage(lif_path)
    scenes = img.scenes

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [
            executor.submit(
                process_scene,
                str(lif_path),
                scene,
                str(output_dir),
                chunk_z,
                bool_imsave,
                Path(plt_path),
                overwrite=overwrite
            )
            for scene in scenes
        ]

        for future in as_completed(futures):
            print(future.result())


# -------------------------------------------------
# Main logic
# -------------------------------------------------
def main(input_path, output_path, other_paths, n_workers=None, chunk_z=True, bool_imsave=False,overwrite=False):

    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
   
    for p in other_paths:
        p.mkdir(parents=True, exist_ok=True)
        if str(p).endswith('plt'):
            plt_dir = p

    if n_workers is None:
        n_workers = max(1, multiprocessing.cpu_count() - 1)

    print(f"Using {n_workers} worker processes")

    if input_path.is_file() and input_path.suffix.lower() == ".lif":

        process_lif_parallel(
            input_path,
            output_path,
            n_workers,
            chunk_z,
            bool_imsave,
            plt_dir,
            overwrite
        )

    elif input_path.is_dir():

        lif_files = list(input_path.glob("*.lif")) + list(input_path.glob("*.LIF"))

        for lif in lif_files:
            lif_output = output_path / lif.stem
            lif_output.mkdir(exist_ok=True)

            process_lif_parallel(
                lif,
                lif_output,
                n_workers,
                chunk_z,
                bool_imsave,
                plt_dir,
                overwrite
            )

    elif input_path.suffix.lower() == ".lif":
        print("File not found: {}".format(input_path))
    else:
        print("Input must be a .lif file or directory.")


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage:")
        print("python Lif_to_Tif.py <input> [n_workers]")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = Path(input_path).parent / 'tifs_out'
    other_paths = [Path(input_path).parent/'cell_profiler_out', Path(input_path).parent/'plt']
        
    workers=None
    bool_imsave=False
    overwrite = False
    if len(sys.argv)>2:
        for i in range(len(sys.argv)-2):
            if sys.argv[2+i].lower().startswith('workers=') or sys.argv[2+i].lower().startswith('w='):
                workers = int(sys.argv[2+i].split('=')[1])
            if [sys.argv[2+i].lower() == s for s in ['imsave','saveimg','image_save','save','bool_imsave=true','imsave=true']]:
                bool_imsave = True
            if sys.argv[2+i].lower() in ['overwrite','overwrite=true','ow=true']:
                overwrite = True

    main(input_path, output_path, other_paths, n_workers=workers, bool_imsave=bool_imsave, overwrite=overwrite)