# How to run the Cell Painting pipeline
### No prior Python experience needed — follow every step in order.

---

## Step 1 — Install Python (once, ever)

1. Go to **https://www.python.org/downloads**
2. Click the big yellow **"Download Python"** button
3. Run the installer
   - ⚠️ On the first screen, tick **"Add Python to PATH"** before clicking Install — this is easy to miss
4. Click **Install Now** and wait until it finishes

To check it worked: open Terminal (Mac) or Command Prompt (Windows) and type `python --version`. You should see something like `Python 3.12.x`.

> **What is Python?** A free programming language. You are installing it the same way you would install any other app.

---

## Step 2 — Our package

The package consists of 3 python (.py) files in the **same folder**:

| File | What it does |
|---|---|
| `cell_paint_pipeline.py` | The analysis script |
| `install_packages.py` | Installs the required tools (run only once before the first use) |
| `metadata_utils.py` | Helper functions (already part of your project) |
| `Lif_to_Tif.py` | Lif to Tif converter used before the Cell_Profiler setup

---

## Step 3 — Install the required packages (once, ever)

Packages are small add-ons that the script needs. You only do this once.

1. Open **Terminal** (Mac) or **Command Prompt** (Windows)
   - **Mac:** press `Cmd + Space`, type `Terminal`, press Enter
   - **Windows:** press `Windows key`, type `cmd`, press Enter
2. Drag and drop **`install_packages.py`** from your file explorer into the Terminal window — the path to the file appears automatically
3. Press **Enter**
4. Wait. You will see a lot of text scrolling past — that is normal. It is done when you see `✅ All packages installed!`

> **What are packages?** Extra libraries of code, like plugins. `pip` is Python's built-in package manager — it downloads them automatically from the internet.

---

## Step 4 — Find the path to your experiment folder

The cell paint analysis script needs to know *where* your data is. This is pre-processed data from CellProfiler (i.e. a database, but not the raw imaging files).
You give it the **folder path** — the address of the folder on your hard drive.

**Mac:**
1. Find your experiment folder in Finder (e.g. `20260319`)
2. Right-click the folder → **"Get Info"**
3. Under "Where" you see the parent path. Your full path is that + `/` + the folder name.
   - Example: `/Volumes/KINGSTON/Nico_data/Cellpaint/20260319`
4. Alternatively: drag the folder into a Terminal window — the path appears automatically ✨

**Windows:**
1. Open the folder in File Explorer
2. Click the address bar at the top — the path is highlighted in blue
3. Copy it (`Ctrl+C`)
   - Example: `C:\Users\Nico\Cellpaint\20260319`

> **What is a path?** Just the address of a folder — like a postal address but for your computer.

---

## Step 5 — Run the pipeline

1. Open Terminal / Command Prompt (same as Step 3)
2. Type `python ` (with a space after), then drag your **`cell_paint_pipeline.py`** file into the Terminal window — the path fills in automatically
3. Add another space, then drag your **experiment folder** into the Terminal window
4. Your command should look like one of these:

**Mac:**
```
python /Users/USERNAME/.../CellPaint/cell_paint_pipeline.py /Volumes/KINGSTON/Nico_data/Cellpaint/20260319
```

**Windows:**
```
python C:\Users\USERNAME\...\CellPaint\cell_paint_pipeline.py C:\Users\Nico\Cellpaint\20260319
```

5. Press **Enter**
6. The script prints its progress. A full run takes roughly 10–30 minutes depending on dataset size.

---

## Step 6 — Find your results

When the script finishes you will see `✅ Pipeline complete.` and the output locations printed. All results are saved **inside your experiment folder**:

```
20260319/
├── plt/          ← all plots (.png files), open with any image viewer
│   ├── clustermap.png
│   ├── treatment_similarity.png
│   ├── tsne_per_cell.png
│   └── ...
└── rslt/         ← data tables (.csv files), open with Excel
    ├── final_profiles.csv
    └── filtered_profiles.csv
```

---

## Running a different experiment

Just repeat Step 5 with a different experiment folder path. Everything else stays the same.

---

## Something went wrong?

| Error message | What it means | Fix |
|---|---|---|
| `python: command not found` | Python is not installed or not on PATH | Redo Step 1, make sure to tick "Add to PATH" |
| `No module named '...'` | A package is missing | Redo Step 3 |
| `No metadata rows found for Name=...` | The folder name doesn't match the `Name` column in your Excel | Open `cell_paint_pipeline.py` in a text editor, find `SELECT_NAMES` near the top, and set it to the exact name from your Excel file |
| `No control treatment found` | The control well label in your Excel doesn't contain "control", "ctrl", or "dmso" | Open the script, find `CONTROL_KEYWORDS`, add the exact word used in your Excel |
| `No such file or directory: CellPaint.db` | The SQLite file is not where expected | Check that `cell_profiler_out/CellPaint.db` exists inside your experiment folder |

---

*Script written for Python 3.10+. Tested on macOS and Windows 11.*
