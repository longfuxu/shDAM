# Interplay Between DNA Polymerase and SSB Analysis Pipeline

This repository contains a pipeline for analyzing optical tweezers (OT) data and kymograph data related to the interplay between DNA Polymerase and SSB proteins. 

## Quick Start: Click the 'Reproduce Run' in the right side to check all the analyzed results.

## Scripts

- **`main.py`**: The entry point script that orchestrates the analysis pipeline. It runs the following scripts in sequence:
  - **`OTdataAnalyzer.py`**: Processes optical tweezers data (e.g., force and distance measurements) to calculate basepairs synthesized over time. Results are saved to the `/results` folder.
  - **`KymographAnalyzer.py`**: Analyzes kymograph imaging data to track DNA polymerase fluorescence. It correlates this with the force data results and performs step-finding on the intensity trace. Results are saved to the `/results` folder.
  - **`SegementsAnalyzer.py`**: Performs the final analysis by segmenting the correlated data into binding/unbinding and exo/pol events. It generates summary plots, heatmaps, and a final segmented data table in the `/results` folder.

## Instructions for Use in Code Ocean

1. **Directory Structure**:
- The `code` folder contains all the necessary Python scripts: `main.py`, `OTdataAnalyzer.py`, `KymographAnalyzer.py`, and `SegementsAnalyzer.py`.
- The `data` folder should contain your input files. The expected structure is:
  - `data/force data/your_force_data.tdms`
  - `data/image data/your_image_data.tdms`
- The scripts will automatically create a `results` folder at the root level to store all output files (plots, spreadsheets, etc.).

2. **Run the Pipeline**:
- Click the "Reproduce Run" button in Code Ocean. This will execute `main.py`, which runs the full analysis pipeline:
1. `OTdataAnalyzer.py` processes the force data.
2. `KymographAnalyzer.py` analyzes the kymograph data, correlating it with the force data results.
3. `SegementsAnalyzer.py` performs final segmentation and correlation analysis on the combined results.
- All intermediate and final results will be saved in the `/results` folder.

3. **Verify Results**:
- Check the `/results` folder for output files, including `.png` plots, `.tiff` images, and `.xlsx` data files.
- Review the console output in Code Ocean for any error messages or completion confirmation.


## Open-Source Code
The open-sourced code for this project is available on GitHub at:  
[https://github.com/longfuxu/Interplay_Between_DNAPol_and_SSB](https://github.com/longfuxu/Interplay_Between_DNAPol_and_SSB)

## Contact
For questions or contributions, please open an issue on the GitHub repository or contact the author directly (longfuxu@berkeley.edu).
