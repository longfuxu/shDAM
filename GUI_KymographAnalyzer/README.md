# DNA Polymerase Trace & Force Correlation GUI

This application provides a comprehensive suite of tools for analyzing and correlating kymograph images with corresponding force-trace data for single-molecule DNA polymerase studies. The GUI allows users to load raw kymograph and trace data, perform semi-automated alignment, extract intensity profiles, detect DNAp binding events, identify change-points in synthesis activity, and generate publication-quality plots and summary data.

![Kymograph Analyzer GUI](../property/kymoanalzyer.png)

## Table of Contents
- [Features](#features)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [How to Run](#how-to-run)
- [GUI Usage Guide](#gui-usage-guide)
  - [1. Setup Tab](#1-setup-tab)
  - [2. Alignment Tab](#2-alignment-tab)
  - [3. Intensity Tab](#3-intensity-tab)
  - [4. Detection Tab](#4-detection-tab)
  - [Main Controls](#main-controls)
- [Example Workflow](#example-workflow)
- [Understanding the Output Files](#understanding-the-output-files)
- [Troubleshooting](#troubleshooting)

## Features

-   Load kymographs from `.tdms` files and force-trace data from `.xlsx` files.
-   Visualize the kymograph and overlay the force-trace trajectory.
-   **Auto-Alignment:** Automatically find the optimal pixel offsets to align the force-trace with the kymograph features.
-   Extract intensity along the aligned DNA polymerase trace.
-   Define a background region of interest (ROI) interactively from the plot.
-   **Intensity Step-Finding:** Use the `AutoStepFinder` algorithm to detect discrete binding/unbinding events from the intensity trace.
-   **Basepair Change-Point Detection:** Analyze the basepairs-vs-time data to find segments of constant synthesis velocity.
-   Generate and save multiple analysis plots, including a final correlated plot showing basepairs, intensity, and binding state over time.
-   Save all analysis parameters and results for reproducibility.

## Project Structure

The application is designed with a modular structure, where the main GUI script relies on separate modules for specific analysis tasks. For the application to function correctly, the following directory layout must be maintained at the root of your project. The main script (`KymographAnalyzer_GUI.py`) depends on relative paths to find its analysis modules.

### Explanation:

-   **`your_project_root_folder/`**: The main project folder that contains all other components. The relative paths in the script are based on this root.
-   **`GUI_KymographAnalyzer/`**: Contains the main `KymographAnalyzer_GUI.py` script. You will execute the application by running this file.
-   **`autostepfinder_GUI/`**: Contains the `AutoStepFinder` algorithm, which is imported by the main GUI for detecting binding/unbinding events from the intensity trace.
-   **`ChangePointDetection_slope_GUI/`**: Contains the algorithm for detecting change-points in the basepair synthesis trace (i.e., finding segments with constant velocity).

When you run an analysis, an output directory (e.g., `analysis_results/`) will be created by the GUI to store all generated files. By default, this directory is created in the location from where you launch the script, but you can specify any location using the "Output Directory" field in the GUI.

## Installation

1.  **Clone the repository** or download the project files and arrange them according to the structure above.

2.  **Install required Python packages.** You can install them all using pip:

    ```sh
    pip install customtkinter pandas numpy opencv-python nptdms scipy tifffile matplotlib plotly kaleido openpyxl
    ```

    -   `customtkinter`: For the modern graphical user interface.
    -   `nptdms`: For reading National Instruments `.tdms` files.
    -   `pandas` & `openpyxl`: For handling Excel trace files.
    -   `numpy` & `scipy`: For numerical and scientific computations.
    -   `opencv-python`: For image processing.
    -   `matplotlib`: For plotting.
    -   `tifffile`: For saving kymographs as TIFF images.
    -   `plotly` & `kaleido`: Used for generating some advanced plots (though the main interface uses Matplotlib).

## How to Run

Navigate to the `GUI_KymographAnalyzer` directory and run the Python script from your terminal:

```sh
python KymographAnalyzer_GUI.py
```

## GUI Usage Guide

The GUI is organized into a control panel on the left (with multiple tabs), a plot panel on the right, and a log panel at the bottom right.

### 1. Setup Tab

This is the starting point for any analysis.

-   **Input Files:**
    -   `Kymo (.tdms)`: Click to select the kymograph data file. This is required.
    -   `Trace (.xlsx)`: Click to select the corresponding force trace data. This is optional, but required for alignment, intensity extraction, and basepair analysis.
-   **Output Directory:** Choose a folder where all results (plots, data files, logs) will be saved. It defaults to `analysis_results` in the current directory.
-   **Kymo Cycle:** An identifier (e.g., `#1`) for the current experiment, which will be included in the output filenames.
-   **Parameter Management:**
    -   `Save Parameters`: Saves all current settings from all tabs into a `.json` file.
    -   `Load Parameters`: Loads settings from a previously saved `.json` file.
-   **Kymograph Display ROI (px):** Sets the initial zoom level for the kymograph plot. Adjust these pixel values to focus on the relevant area.
-   **Load Data & Visualize Kymo Button:** After selecting your files, click this button. It loads the data, saves a `.tiff` version of the kymograph in your output folder, and displays the kymograph in the plot panel. This enables the analysis buttons in other tabs.

### 2. Alignment Tab

This tab is for aligning the force trace trajectory with the kymograph image.

-   **Initial Offsets (px):** If the trace overlay is misaligned, you can provide initial `X` and `Y` pixel offsets here to roughly correct it.
-   **Auto-Alignment Search Range (px):** Defines the search area (in pixels) around the initial offset for the automatic alignment algorithm.
-   **Find Optimal Alignment Button:** Runs an algorithm that searches for the offset that maximizes the kymograph intensity along the trace path. The `Initial Offsets` will be updated with the optimal values, and the plot will show the new alignment in cyan.

### 3. Intensity Tab

Here, you extract and analyze the intensity profile along the aligned trace.

-   **Signal Filter Window:** A window size for the Savitzky-Golay filter to smooth the extracted intensity signal. Must be an odd integer > 3.
-   **Background ROI (px):** Defines a rectangle (`rec_x`, `rec_y`, `rec_w`, `rec_h`) on the kymograph to be used for background noise calculation.
    -   You can set these values manually or click `Select BG ROI on Plot`. This will let you draw a rectangle directly on the kymograph.
-   **Background Analysis:**
    -   `BG Filter Window`: Smoothing window for the background intensity.
    -   `Threshold (σ)`: Sets the threshold for identifying a "bound" state. The threshold is calculated as `mean(background) + σ * stdev(background)`.
-   **Analyze Intensity & Background Button:** Performs the intensity extraction and background analysis. It plots the resulting intensity trace over time, with the calculated threshold shown as a red dashed line.

### 4. Detection Tab

This tab contains tools for event detection in both the intensity and basepair data.

-   **Intensity Step-Finding:** Parameters for the `AutoStepFinder` algorithm.
    -   `S_max Thresh`, `Error Tolerance`: Core parameters controlling the sensitivity of step detection.
    -   `Min Event (pts)`: A post-processing filter to remove detected on/off states that are shorter than this number of data points.
-   **Run Intensity Step-Finding Button:** Runs the algorithm and overlays the detected steps (a "staircase" fit) on the intensity plot.
-   **Basepair Change-Point Detection:** Parameters for finding points where the DNA synthesis rate changes.
    -   `Window Size`: The size of the sliding window used to calculate velocity.
    -   `Sigma Noise`: An estimate of the noise in the velocity signal.
-   **Run BP Change-Point Detection Button:** Runs the analysis and plots the basepair data with linear fits over the detected segments.

### Main Controls

-   **Run Full Analysis & Generate Report:** The primary "one-click" button. It sequentially runs the entire pipeline (Intensity Analysis -> Step Finding -> BP Change-Point Detection) using the current parameters. It then generates and saves the final summary plots and data files.
-   **Save Current Plot:** Saves the plot currently displayed in the plot panel to a `.png` or `.pdf` file.

## Example Workflow

1.  **Launch the GUI:** `python KymographAnalyzer_GUI.py`.
2.  **Setup:**
    -   In the **Setup** tab, select your `kymo.tdms` and `trace.xlsx` files.
    -   Choose an output directory.
    -   Click `Load Data & Visualize Kymo`. The kymograph appears on the right.
3.  **Align:**
    -   Go to the **Alignment** tab.
    -   Click `Find Optimal Alignment`. Wait for the process to finish. The GUI will log the optimal offsets and update the plot with a cyan trace overlay.
4.  **Analyze Intensity:**
    -   Go to the **Intensity** tab.
    -   Click `Select BG ROI on Plot`. On the kymograph, click and drag to draw a rectangle in an area with no signal.
    -   Click `Analyze Intensity & Background`. The plot view will switch to the intensity trace.
5.  **Run Detections (Optional):**
    -   To see the step-fitting, go to the **Detection** tab and click `Run Intensity Step-Finding`.
    -   To see the basepair segmentation, click `Run BP Change-Point Detection`.
6.  **Generate Final Report:**
    -   Click the main `Run Full Analysis & Generate Report` button at the bottom of the control panel.
    -   The entire pipeline will run, and all output files will be saved to your specified output directory. A confirmation message will appear when done. Check the log panel for details of the process.

## Understanding the Output Files

When you run a full analysis, the following files will be generated in your output directory (using `my_kymo` and cycle `#1` as an example):

-   `my_kymo.tiff`: A standard TIFF image of your kymograph.
-   `stepfinding_results/my_kymo-cycle#1-Intensity_along_DNAp-filtered_fit_trace.csv`: A CSV file with the time, raw intensity data, and the step-fitted intensity data.
-   `stepfinding_results/my_kymo-cycle#1-Intensity_along_DNAp-filtered_step_properties.csv`: A CSV with detailed properties for each detected step (start/end time, amplitude, etc.).
-   `ChangePoints_Results/my_trace-correlated_analysis.xlsx`: An Excel file containing the raw trace data and a sheet with the start/end times and basepair positions for each detected synthesis segment.
-   `my_kymo-cycle#1-final_overlay.png`: The kymograph image with the final aligned trace and background ROI overlaid.
-   `my_kymo-cycle#1-final_correlated_plot.png`: A comprehensive plot showing basepair synthesis, intensity (raw and fitted), and the binarized "bound" state all on a shared time axis. This is the key summary plot.
-   `my_kymo-cycle#1-analysis_log.txt`: A text file containing all the parameters used for the analysis, for reproducibility.

## Troubleshooting

-   **Parameter Error:** If you see a "Parameter Error" popup, it usually means a non-numeric value was entered, or a filter window was not an odd number. Please check the values in all tabs.
-   **Auto-Alignment Failed:** This can happen if the initial trace position is too far from the real signal. Try adjusting the `Initial Offsets` in the **Alignment** tab to be closer to the bright line in the kymograph before running `Find Optimal Alignment`. Increasing the `Search Range` may also help.
-   **Script Crashes or Freezes:** Check the terminal where you ran the script for any error messages. Ensure all dependencies are correctly installed.