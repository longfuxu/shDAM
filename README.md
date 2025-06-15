# Single-molecule high-resolution dynamics and activity mapping (shDAM) of DNA-binding proteins

This repository provides a comprehensive suite of tools for analyzing DNA polymerase dynamics at the replication fork, using data from fluorescence microscopy and optical tweezers. It offers multiple ways to process and visualize data, including Jupyter notebooks for exploration, Python scripts for automation, and Graphical User Interfaces (GUIs) for interactive analysis. All analysis scripts process mechanical data and associated fluorescence kymographs from `.tdms` files.

An example dataset is provided in the `example_dataset` folder for a complete walkthrough of the analysis pipeline.

## Analysis Workflow

The core analysis pipeline consists of three main steps, which correlate mechanical data from optical tweezers with fluorescence microscopy data:

1.  **Optical Tweezer (OT) Data Analysis**: This step processes raw mechanical data (from `.tdms` files) to determine the enzyme's activity. It involves calculating the corrected end-to-end distance (EED) of the DNA tether using polymer elasticity models (e.g., eWLC and FJC). This EED is then used to derive the high-resolution position of the enzyme's active site over time, creating a "DNA kymograph" that tracks activity in nanometers. Enzymatic rates (e.g., polymerization, exonuclease) and pausing events are identified from this mechanical data.

2.  **Kymograph Analysis and Correlation**: This step analyzes fluorescence kymographs to detect fluorescently labeled protein binding/unbinding. The mechanically-derived "DNA kymograph" from Step 1 is temporally and spatially overlaid onto the fluorescence kymograph. This alignment is critical for precisely correlating enzymatic activity with the presence of a protein.

3.  **Correlated Segment Analysis**: This final step performs an in-depth analysis of the synchronized data. The fluorescence signal at the activity site is extracted and binarized into "bound" (on) and "unbound" (off) states. Similarly, the mechanical data is segmented into "active" and "paused" states. By combining these, each time point is categorized into one of four states (e.g., fluorescent and active, fluorescent and paused, etc.), enabling detailed statistical analysis of binding lifetimes, activity bursts, and processivity.

## Installation and Setup

### Prerequisites
- Python 3.9
- Git

### Installation Steps

1.  Clone the repository:
    ```bash
    git clone https://github.com/longfuxu/shDAM.git
    cd shDAM
    ```

2.  Create and activate a virtual environment:
    ```bash
    # On macOS/Linux
    python -m venv venv
    source venv/bin/activate

    # On Windows
    python -m venv venv
    .\venv\Scripts\activate
    ```

3.  Install the required packages:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

You can run the analysis using Jupyter Notebooks, Python scripts, or interactive GUIs. For a first-time user, we recommend starting with the Jupyter Notebooks and the provided example dataset.

### 1. Jupyter Notebooks
Ideal for a step-by-step walkthrough of the analysis pipeline.

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

## Reproducibility

### Parameter Tuning
The parameters in the notebooks, scripts, and GUIs are tuned for the provided example dataset. For your own data, these parameters will likely need to be adjusted to fit your specific experimental conditions and research system.

### Code Ocean
For fully reproducible results, a Code Ocean capsule is available:
- [Link to Code Ocean Capsule] (To be added upon publication)

## Additional Algorithms
This project utilizes and provides fast, custom-built algorithms for single-molecule data analysis. These tools are included as submodules or directories within this repository. For more details, please see the specific directories:
- **`fast_pwl_fit_GUI`**: A tool for fast piece-wise linear fitting of data traces.
- **`autostepfinder_GUI`**: A robust step-finding algorithm used to identify discrete steps.
- **`ChangePointDetection_slope_GUI`**: A efficient changepoint detection algorithm used to identify change points for gradual changing traces.


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