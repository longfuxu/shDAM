# Single-molecule high-resolution dynamics and activity mapping (shDAM) of DNA-binding proteins

## 🚀 Quick Start Guide

**New to this analysis?** Follow these steps:

1. **Mac/Linux users**: Use **pip installation** (see [Method 1](#method-1-pip-installation-maclinux) below)
2. **Windows users**: Use **Conda installation** (see [Method 2](#method-2-conda-installation-recommended) below) for the easiest setup
3. **Test your installation** with the [verification commands](#verification)
4. **Run the example**: start with `1_CalculatingDNApTrace_OT.ipynb`

---

## About

This repository provides a comprehensive suite of tools for analyzing DNA polymerase dynamics at the replication fork, using data from fluorescence microscopy and optical tweezers. It offers multiple ways to process and visualize data, including Jupyter notebooks for exploration, Python scripts for automation, and Graphical User Interfaces (GUIs) for interactive analysis. All analysis scripts process mechanical data and associated fluorescence kymographs from `.tdms` files.

An example dataset is provided in the `example_dataset` folder for a complete walkthrough of the analysis pipeline.

## Analysis Workflow

The core analysis pipeline consists of three main steps, which correlate mechanical data from optical tweezers with fluorescence microscopy data:

1.  **Optical Tweezer (OT) Data Analysis**: This step processes raw mechanical data (from `.tdms` files) to determine the enzyme's activity. It involves calculating the corrected end-to-end distance (EED) of the DNA tether using polymer elasticity models (e.g., eWLC and FJC). This EED is then used to derive the high-resolution position of the enzyme's active site over time, creating a "DNA kymograph" that tracks activity in nanometers. Enzymatic rates (e.g., polymerization, exonuclease) and pausing events are identified from this mechanical data.

2.  **Kymograph Analysis and Correlation**: This step analyzes fluorescence kymographs to detect fluorescently labeled protein binding/unbinding. The mechanically-derived "DNA kymograph" from Step 1 is temporally and spatially overlaid onto the fluorescence kymograph. This alignment is critical for precisely correlating enzymatic activity with the presence of a protein.

3.  **Correlated Segment Analysis**: This final step performs an in-depth analysis of the synchronized data. The fluorescence signal at the activity site is extracted and binarized into "bound" (on) and "unbound" (off) states. Similarly, the mechanical data is segmented into "active" and "paused" states. By combining these, each time point is categorized into one of four states (e.g., fluorescent and active, fluorescent and paused, etc.), enabling detailed statistical analysis of binding lifetimes, activity bursts, and processivity.

## Installation and Setup

### Prerequisites
- Python 3.9 or 3.10 (recommended)
- Git (for cloning with submodules)

---

## Installation Methods

### Method 1: Pip Installation (Mac/Linux)

1.  **Clone the repository with submodules:**
    ```bash
    git clone --recurse-submodules https://github.com/longfuxu/shDAM.git
    cd shDAM
    ```

2.  **Create and activate a virtual environment:**
    ```bash
    # On macOS/Linux
    python3 -m venv shdam
    source shdam/bin/activate
    ```

3.  **Install packages:**
    ```bash
    pip install --upgrade pip
    pip install -r requirements.txt
    ```

---

### Method 2: Conda Installation (Recommended)

**✅ This method works on all platforms (Windows, Mac, Linux), especially for Windows users!**

1.  **Install Miniconda:**
    - Download from: https://docs.conda.io/en/latest/miniconda.html
    - Follow the installation instructions for your operating system (ensure to add to your system PATH)

2.  **Clone the repository with submodules:**
    ```bash
    git clone --recurse-submodules https://github.com/longfuxu/shDAM.git
    cd shDAM
    ```

3.  **Initialize conda for your shell (first-time setup):**
    ```bash
    # For PowerShell (Windows)
    conda init powershell
    
    # For Command Prompt (Windows)
    conda init cmd.exe
    
    # For Mac/Linux (bash/zsh)
    conda init
    ```
    
    **⚠️ Important**: After running `conda init`, **close and reopen your terminal** for changes to take effect.

4.  **Create and activate conda environment:**
    ```bash
    # Create environment with Python and essential packages
    conda create -n shdam python=3.9 numpy scipy pandas matplotlib opencv numba
    conda activate shdam
    ```
    
    **✅ Verify activation**: You should see `(shdam)` at the beginning of your prompt when the environment is active.

5.  **Install additional packages with conda:**
    ```bash
    # Install packages available in conda
    conda install jupyterlab seaborn sympy openpyxl
    ```

6.  **Install remaining packages with pip:**
    ```bash
    # Install packages only available via pip or requiring specific versions
    pip install plotly>=5.18.0 kaleido==0.2.1
    pip install nptdms tifffile more-itertools pwlf tabulate ipywidgets
    
    # For GUI applications
    pip install customtkinter ttkbootstrap
    ```

---

## Verification

After installation, test your setup:

```bash
# Run the comprehensive test script
python test_installation.py
```

This script will check all packages and provide a detailed report of what's working.

---

## Troubleshooting

### Common Issues

| Problem | Solution |
|---------|----------|
| **Missing submodules** | Use `git clone --recurse-submodules` or run `git submodule update --init --recursive` |
| **Package conflicts** | Use conda environment or virtual environment |
| **Compilation errors** | Use **Method 2 (Conda)** which provides pre-compiled packages |
| **Permission errors (Windows)** | Run terminal as Administrator |

### Conda Environment Issues

#### Problem: `conda` command not recognized
```
conda: The term 'conda' is not recognized...
```

**Solutions:**
1. **Initialize conda for your shell:**
   ```bash
   # Find your conda installation path and run:
   C:\Users\YourUsername\miniconda3\Scripts\conda.exe init powershell
   # Then close and reopen your terminal
   ```

2. **Use Anaconda Prompt instead:**
   - Search for "Anaconda Prompt" in Windows Start menu
   - Run all conda commands from there

#### Problem: Environment not properly activated
Even after running `conda activate shdam`, you don't see `(shdam)` in your prompt.

**Check activation status:**
```bash
# This should show the conda environment path, not system Python
python -c "import sys; print('Python path:', sys.executable)"

# Should show: C:\Users\YourUsername\miniconda3\envs\shdam\python.exe
# NOT: C:\Users\YourUsername\AppData\Local\Programs\Python\...
```

**Solutions:**
1. **Re-initialize conda:** `conda init powershell` then restart terminal
2. **Use full path:** `C:\Users\YourUsername\miniconda3\Scripts\activate shdam`
3. **Use Anaconda Prompt:** Always works with conda commands

### Package-Specific Issues

| Package | Issue | Solution |
|---------|-------|----------|
| `numba` | Compilation errors | Use `conda install numba` |
| `jupyter` | Build failures | Use `conda install jupyter` |
| `matplotlib` | Missing dependencies | Use `conda install matplotlib` |
| `plotly` + `kaleido` | Version compatibility errors | Use `pip install kaleido==0.2.1` (not conda) |

### Jupyter Notebook Matplotlib Backend Issues

#### Problem: `%matplotlib widget` not working
```python
ValueError: ... matplotlib ... widget backend error
```

**Solutions - Try these matplotlib backends in order:**

```python
# Option 1: Interactive widget backend (best for zooming/panning)
%matplotlib widget

# Option 2: Inline static plots (most reliable)
%matplotlib inline

# Option 3: Interactive plots in separate window
%matplotlib qt

# Option 4: Interactive plots inline (newer method)
%matplotlib ipympl
```

**Choose the one that works for your system:**
- **`%matplotlib inline`**: Most reliable, works everywhere (static plots)
- **`%matplotlib widget`**: Best interactivity but requires `ipywidgets`
- **`%matplotlib qt`**: Good for detailed analysis but opens separate windows
- **`%matplotlib ipympl`**: Modern alternative to widget backend

**Note**: Always restart your Jupyter kernel after changing backends!

### Kaleido Path Issues (Windows/Mac with spaces in paths)

#### Problem: Plotly image export fails with path errors
```python
ValueError: Failed to start Kaleido subprocess. Error stream:
cd: /path/with spaces/: No such file or directory
```

**Solution**: The notebooks and scripts have been updated to use matplotlib for image saving instead of Plotly's Kaleido engine. This ensures:
- ✅ **Cross-platform compatibility** (Windows, Mac, Linux)
- ✅ **Works with spaces in directory names**
- ✅ **Same quality output** (300 DPI)
- ✅ **Interactive Plotly displays preserved** for exploration

**Manual Fix**: If you encounter this error, replace `fig.write_image()` calls with matplotlib equivalents.

### General Tips
- **Use Conda**: It handles compiled packages automatically across all platforms
- **Virtual environments**: Always use them to avoid conflicts
- **Python version**: Stick to Python 3.9 or 3.10 for best compatibility

---

## Getting Help

If you're still having issues:
1. Check the error message carefully
2. Try the **Conda installation** (Method 2)
3. Contact [longfu2.xu@gmail.com](mailto:longfu2.xu@gmail.com)
4. Include your operating system and Python version in any support request

---

## Next Steps After Installation

Once your installation is verified, you can:

### 1. 🎯 Start with the Example Dataset
```bash
# Navigate to the example folder
cd example_dataset

# Open Jupyter Lab (if installed)
jupyter lab

# Or run Python scripts directly
python ../python_scripts/1_CalculatingDNApTrace_OT.py
```

### 2. 📚 Follow the Analysis Pipeline
The complete analysis consists of three main steps:
1. **Optical Tweezer Data Analysis** (`1_CalculatingDNApTrace_OT.ipynb`)
2. **Kymograph Correlation** (`2_Correlation_image_force.ipynb`) 
3. **Correlated Segment Analysis** (`3_Correlated_segement_analysis.ipynb`)

### 3. 🖥️ Try the GUI Tools
For interactive analysis, explore the GUI applications:
- `GUI_OTdataAnalyzer/` - For optical tweezer data
- `GUI_KymographAnalyzer/` - For kymograph analysis
- `GUI_SegmentsAnalyzer/` - For segment analysis

### 4. 📖 Learn More
- Read the detailed methodology in the [paper](https://doi.org/10.1038/s41467-024-49612-3)
- Check out the advanced algorithms in `autostepfinder_GUI/` and `ChangePointDetection_slope_GUI/`
- Explore the `property/` folder for example outputs

## Usage

You can run the analysis using Jupyter Notebooks, Python scripts, or interactive GUIs. For a first-time user, we recommend starting with the Jupyter Notebooks and the provided example dataset.

### 1. Jupyter Notebooks
Ideal for a step-by-step walkthrough of the ana lysis pipeline.

1.  Start JupyterLab:
    ```bash
    jupyter lab
    ```
2.  Navigate to and run the notebooks in order:
    - `1_CalculatingDNApTrace_OT.ipynb`
    - `2_Correlation_image_force.ipynb`
    - `3_Correlated_segement_analysis.ipynb`

### 2. Python Scripts
Located in the `python_scripts/` directory, these are suitable for batch processing and automation. To run the analysis, execute the scripts in order:
```bash
python3 python_scripts/1_CalculatingDNApTrace_OT.py
python3 python_scripts/2_Correlation_image_force.py
python3 python_scripts/3_Correlated_segement_analysis.py
```

### 3. Graphical User Interfaces (GUIs)
For users who prefer interactive controls for data analysis. We provide GUIs for the key analysis steps:
- **GUI_OTdataAnalyzer**: For processing optical tweezer data.
- **GUI_KymographAnalyzer**: For analyzing kymographs and finding steps.
- **GUI_SegmentsAnalyzer**: For analyzing correlated segments.

![OT Data Analyzer](property/OTAnalyzer_screenshot.png)
_The GUI for Optical Tweezer data analysis._

![Kymograph Analyzer](property/kymoanalzyer.png)
_The GUI for Kymograph analysis._

![Segments Analyzer](property/SegmentsAnalyzer_screenshot.png)
_The GUI for Correlated Segments analysis._

## Reproducibility

### Parameter Tuning
The parameters in the notebooks, scripts, and GUIs are tuned for the provided example dataset. For your own data, these parameters will likely need to be adjusted to fit your specific experimental conditions and research system.

### Code Ocean
For fully reproducible results, a Code Ocean capsule is available:
- [Link to Code Ocean Capsule] (To be added upon publication)

## Advanced Algorithms for Data Analysis

This project integrates several powerful, standalone tools for specialized single-molecule data analysis. Each tool is available in its own directory and includes a graphical user interface (GUI) for ease of use.

### `autostepfinder_GUI`: Automated Step-Finding for Discrete Events
This tool provides an implementation of the `AutoStepfinder` algorithm, designed to rapidly and automatically detect discrete steps in time-series data. It is particularly well-suited for analyzing single-molecule traces where proteins or enzymes exhibit step-like movements or changes in state. The algorithm uses a dual-pass method to determine the optimal number of steps, providing a robust fit even for noisy data.

**Key Features**:
- Automated, dual-pass step detection.
- Interactive GUI to tune parameters and visualize fits in real-time.
- Exports step properties, including size, dwell time, and error.

**Original Algorithm Credit**:
This tool is based on the MATLAB code from the 2021 publication by Loeff et al. Please cite the original paper when using this algorithm:
> Loeff, L., Kerssemakers, J. W. J., Joo, C., & Dekker, C. (2021). AutoStepfinder: A fast and automated step detection method for single-molecule analysis. *Patterns*, Volume 2, Issue 5, 100256

### `ChangePointDetection_slope_GUI`: Detecting Changes in Gradual Traces
This tool is designed to identify significant changes in the slope of a data trace, making it ideal for analyzing processes with gradual transitions rather than abrupt steps. It implements the change-point detection algorithm described by Kerssemakers et al. (2006), which is effective for parsing phenomena such as the assembly/disassembly dynamics of biomolecules.

**Key Features**:
- Detects change-points in gradually varying data.
- Allows for fine-tuning of sensitivity via moving window size and noise parameters.
- Visualizes the original data, fitted segments, and the derivative trace.

**Original Algorithm Credit**:
The underlying method was first described in the 2006 Nature paper by Kerssemakers et al. Please cite this work when using the change-point detection algorithm:
> Kerssemakers, J., Munteanu, E., Laan, L., Noetzel, T. L., Janson, M. E., & Dogterom, M. (2006). Assembly dynamics of microtubules at molecular resolution. *Nature*, 442(7103), 709–712.

### `fast_pwl_fit_GUI`: Fast Piece-Wise Linear Fitting
This tool offers a fast algorithm for fitting time-series data with a series of connected straight lines (piece-wise linear). It is a general-purpose utility for segmenting traces into linear portions, which can be useful for identifying periods of constant velocity or rate.

## Example Walkthrough

The following images illustrate the key outputs from each stage of the analysis pipeline, run on the provided `example_dataset`.

#### Step 1: Optical Tweezer Data Analysis (`1_CalculatingDNApTrace_OT.ipynb`)
The initial step processes the mechanical data. It involves fitting the force-extension data to determine the tether's properties, calculating the percentage of ssDNA over time, and deriving the high-resolution DNA junction position, which represents the enzyme's activity trace.
![Cycle Fit](property/1_cycle-fit.png)
![DNAp Traces](property/1_DNAp-trace.png)

#### Step 2: Kymograph Correlation (`2_Correlation_image_force.ipynb`)
Here, the mechanically-derived activity trace is correlated with the fluorescence kymograph. The fluorescence intensity along the enzyme's path is extracted, filtered, and binarized using a step-finding algorithm to identify protein binding events. All traces can then be visualized together.
![Protein Trajectory](property/2_protein_trace.jpg)
![Correlated Traces](property/2_correlated_trace.png)
![Filtered and Binarized Intensity](property/2_DNAp-intensity.png)
![Step Intensity along DNAp](property/2_DNAp-intensity_stepwise.png)
![Final Correlated Segemnts](property/2_correlated_segments.png)

#### Step 3: Segment Analysis and Visualization (`3_Correlated_segement_analysis.ipynb`)
The final step is to segment the continuous data into discrete states (e.g., active vs. paused, fluorescent vs. non-fluorescent). This allows for detailed statistical analysis, such as identifying activity bursts and their correlation with fluorescence, and generating summary visualizations like heatmaps.
![Activity Bursts](property/3_activity_burst.png)
![Exo Segments](property/3_exo-segments.png)
![Correlation Heatmap](property/3_heatmap.png)

## Contributing
We welcome contributions to enhance and expand this project. Please fork the repository, make your changes, and submit a pull request. You can also contact Dr.Longfu Xu or Prof. Gijs Wuite for contributions.

## Support and Contact
Please note that the code in this repository is custom-written for internal lab use and may still contain bugs. For questions, support, or feedback, please contact Dr. Longfu Xu at [longfu2.xu@gmail.com](mailto:longfu2.xu@gmail.com).

## Citation
Xu, L., Halma, M.T.J. & Wuite, G.J.L. Mapping fast DNA polymerase exchange during replication. *Nat Commun* **15**, 5328 (2024). https://doi.org/10.1038/s41467-024-49612-3

## License
This project is licensed under the MPL-2.0 license. See the `LICENSE` file for more details.

## Acknowledgments
All code in this repository was developed by Dr. Longfu Xu (longfuxu.com) during his PhD at the [Gijs Wuite Lab](http://www.gijswuite.com/). 