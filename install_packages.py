"""
Run this once to install all required packages.
Double-click it in your file explorer, or run: python install_packages.py
"""
import subprocess
import sys

packages = [
    # data & analysis
    "numpy",
    "pandas",
    "scikit-learn",
    "pycytominer",
    "umap-learn",
    "pyarrow",
    "openpyxl",
    # plotting
    "matplotlib",
    "seaborn",
    "cmap",
    # image I/O
    "aicsimageio",
    "aicspylibczi",
    "dask",
    # utilities
    "tqdm",
]

print("Installing packages... this may take a few minutes.\n")
for pkg in packages:
    print(f"  Installing {pkg}...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])

print("\n✅ All packages installed! You are ready to run the pipeline.")
input("\nPress Enter to close this window.")
