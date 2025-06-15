#!/usr/bin/env python
# coding: utf-8
"""
Correlation Image Force Analysis GUI

This script provides a graphical user interface for the correlation analysis 
between DNA polymerase trace from images and force measurements, including 
change-point detection and visualization.

Required packages:
pip install customtkinter pandas numpy opencv-python nptdms scipy tifffile matplotlib kaleido
"""
# --- Core Libraries ---
import os
import sys
import json
import numpy as np
import cv2
import pandas as pd
from nptdms import TdmsFile
from scipy.signal import savgol_filter
from scipy import interpolate
from scipy.interpolate import interp1d
import tifffile as tif
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import cm
from matplotlib.widgets import RectangleSelector

# --- GUI Libraries ---
# Ensure you have customtkinter installed: pip install customtkinter
import customtkinter as ctk
from tkinter import filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import plotly.graph_objects as go
import plotly.express as px

# --- Analysis Modules (Bundled with the script) ---
# Make the module paths robust to where the script is run from
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
autostepfinder_path = os.path.join(project_root, 'autostepfinder_GUI')
changepoint_path = os.path.join(project_root, 'ChangePointDetection_slope_GUI')
sys.path.append(autostepfinder_path)
sys.path.append(changepoint_path)

from autostepfinder import AutoStepFinder, StepFinderParameters
from bp_detection.bp_batch_segments import bp_batch_segments_for


class AnalysisPipeline:
    """Handles the data loading, processing, and analysis logic."""
    def __init__(self, logger_func=print):
        self.data = {}
        self.results = {}
        self.log = logger_func

    def load_data(self, kymo_file, trace_file, output_dir):
        """Loads TDMS and trace data."""
        self.log(f"Loading kymograph: {kymo_file}")
        
        metadata = TdmsFile.read_metadata(kymo_file)
        width = metadata.properties['Pixels per line']
        tdms_file = TdmsFile(kymo_file)
        height = len(tdms_file['Data']['Time (ms)'][:]) / width
        
        chn_r = np.array([int(i) for i in tdms_file['Data']['Pixel ch 1'][:]])
        chn_g = np.array([int(i) for i in tdms_file['Data']['Pixel ch 2'][:]])
        chn_b = np.array([int(i) for i in tdms_file['Data']['Pixel ch 3'][:]])

        chn_rgb = np.vstack((chn_r, chn_g, chn_b)).T
        img = chn_rgb.reshape((int(height), int(width), 3)).transpose((1, 0, 2)).astype(np.uint16)
        
        self.data['img'] = img
        self.data['g'] = cv2.split(img)[1]
        self.data['metadata'] = metadata
        self.data['time_per_line'] = tdms_file['Data']['Time (ms)'][-1] / height
        self.data['px_size'] = float(metadata.properties['Scan Command.Scan Command.scanning_axes.0.pix_size_nm'])
        self.log("Kymograph loaded successfully.")

        if os.path.exists(trace_file):
            self.log(f"Loading trace file: {trace_file}")
            self.data['trace'] = pd.read_excel(trace_file, engine='openpyxl')
            self.log("Trace file loaded.")
        else:
            self.data['trace'] = None
            self.log("No trace file provided or found. Some features will be disabled.")
            
        results_dir = os.path.join(output_dir, "results")
        os.makedirs(results_dir, exist_ok=True)
        file_stem = os.path.splitext(os.path.basename(kymo_file))[0]
        tif_path = os.path.join(results_dir, f"{file_stem}.tiff")
        tif.imwrite(tif_path, img)
        return tif_path

    def run_auto_alignment(self, params):
        self.log("Starting automatic alignment...")
        g = self.data['g']
        trace = self.data['trace']
        x_cali = 1000 / self.data['time_per_line']
        y_cali = 1000 / self.data['px_size']

        x_offset_ls, y_offset_ls, intens_sum_ls = [], [], []
        x_range = np.arange(-params['x_search'], params['x_search'] + 1)
        y_range = np.arange(-params['y_search'], params['y_search'] + 1)

        for i_x in x_range:
            for i_y in y_range:
                x_offset = params['x_offset'] + i_x
                y_offset = params['y_offset'] + i_y
                
                trace_time = trace['time'] / 1000 * x_cali + x_offset
                position = pd.to_numeric(trace['junction_position_all'], errors='coerce') * y_cali + y_offset
                trace_time, position = trace_time.dropna(), position.dropna()

                if len(trace_time) < 101: continue
                
                func = interpolate.interp1d(trace_time, position, kind='slinear', fill_value="extrapolate")
                x = np.arange(round(trace_time.iloc[100]), round(trace_time.iloc[-1]))
                y = np.rint(func(x)).astype(int)
                
                y = np.clip(y, 0, g.shape[0] - 1)
                x_valid_indices = (x >= 0) & (x < g.shape[1])
                x = x[x_valid_indices]
                y = y[x_valid_indices]

                intens_sum = np.sum(g[y, x])
                
                x_offset_ls.append(x_offset)
                y_offset_ls.append(y_offset)
                intens_sum_ls.append(intens_sum)

        if not intens_sum_ls:
            raise RuntimeError("Auto-alignment failed. Could not find valid trace segments. Try adjusting initial offset or search range.")

        intens_ind = np.argmax(intens_sum_ls)
        x_offset_optimal = x_offset_ls[intens_ind]
        y_offset_optimal = y_offset_ls[intens_ind]
        self.log(f"Optimal alignment found: x_offset={x_offset_optimal}, y_offset={y_offset_optimal}")
        return x_offset_optimal, y_offset_optimal

    def calculate_intensity_trace(self, params):
        g = self.data['g']
        trace = self.data['trace']
        x_cali = 1000 / self.data['time_per_line']
        y_cali = 1000 / self.data['px_size']
        x_offset, y_offset = params['x_offset'], params['y_offset']

        trace_time = trace['time'] / 1000 * x_cali + x_offset
        position = pd.to_numeric(trace['junction_position_all'], errors='coerce') * y_cali + y_offset
        trace_time, position = trace_time.dropna(), position.dropna()

        func = interpolate.interp1d(trace_time, position, kind='slinear', fill_value="extrapolate")
        x = np.arange(round(trace_time.iloc[0]), round(trace_time.iloc[-1]))
        x1 = (x - x_offset) / x_cali
        y = np.rint(func(x)).astype(int)
        
        y_range = [y - 2, y - 1, y, y + 1, y + 2]
        all_intensity = np.zeros_like(x, dtype=float)
        for y_coord in y_range:
            y_clipped = np.clip(y_coord, 0, g.shape[0] - 1)
            x_clipped = np.clip(x, 0, g.shape[1] - 1)
            all_intensity += g[y_clipped, x_clipped]
        
        signal_filter = savgol_filter(all_intensity, params['signal_filter_window'], 3)
        
        self.results.update({
            'x_px': x, 'y_px': y, 'time_s': x1, 
            'raw_intensity': all_intensity, 'filtered_intensity': signal_filter,
            'x_cali': x_cali, 'y_cali': y_cali, 'trace_time_px': trace_time, 'position_px': position
        })
        self.log("Intensity trace extracted and filtered.")

    def analyze_background(self, params):
        g = self.data['g']
        rec_x_list = np.arange(params['rec_x'], params['rec_x'] + params['rec_w'])
        
        if rec_x_list[-1] >= g.shape[1] or params['rec_y'] + params['rec_h'] > g.shape[0]:
            self.log("Warning: Background ROI is partially outside the image. It will be clipped.")
        
        y_coords = np.arange(params['rec_y'], params['rec_y'] + params['rec_h'])
        y_coords = np.clip(y_coords, 0, g.shape[0]-1)
        rec_x_list = np.clip(rec_x_list, 0, g.shape[1]-1)

        bagrnd_all_intensity = np.sum([g[y_val, rec_x_list] for y_val in y_coords], axis=0)
        bagrnd_signal_filter = savgol_filter(bagrnd_all_intensity, params['bg_filter_window'], 3)
        self.results['bagrnd_signal_filter'] = bagrnd_signal_filter
        self.results['bagrnd_raw_intensity'] = bagrnd_all_intensity
        self.results['bagrnd_time_px'] = rec_x_list
        
        avg, std = np.average(bagrnd_signal_filter), np.std(bagrnd_signal_filter)
        self.log(f"Background analysis done. Mean: {avg:.2f}, StdDev: {std:.2f}")

    def run_intensity_stepfinder(self, params, output_dir, file_stem):
        time_data = self.results['time_s']
        signal_data = self.results['filtered_intensity']

        # Save the filtered intensity data to a .txt file, matching the original script
        results_dir = os.path.join(output_dir, "results")
        os.makedirs(results_dir, exist_ok=True)
        df_to_save = pd.DataFrame({'time serials': time_data, 'filtered intensity': signal_data})
        txt_path = os.path.join(results_dir, f"{file_stem}-cycle#{params['kymo_cycle']}-Intensity_along_DNAp-filtered.txt")
        np.savetxt(txt_path, df_to_save.values, fmt='%1.3f')
        self.log(f"Saved filtered intensity data to: {txt_path}")
        
        output_dir_steps = os.path.join(output_dir, 'stepfinding_results')
        os.makedirs(output_dir_steps, exist_ok=True)
        
        time_res = np.mean(np.diff(time_data))
        asf_params = StepFinderParameters(s_max_threshold=params['s_max_thresh'], resolution=time_res, fit_range=10000, fit_mode='mean', local_step_merge=True, error_tolerance=params['error_tol'], overshoot=1)
        finder = AutoStepFinder(params=asf_params)
        final_fit, final_steps, _, _ = finder.run(signal_data)
        self.log(f"AutoStepFinder found {len(final_steps)} steps.")

        self.results['final_fit'] = final_fit
        self.results['final_steps'] = final_steps

        threshold = np.mean(self.results['bagrnd_signal_filter']) + (params['threshold_sigma'] * np.std(self.results['bagrnd_signal_filter']))
        signal_binarized = np.where(final_fit > threshold, 1, 0)
        
        def filter_short_events(binary_trace, min_length_points):
            filtered_trace = np.copy(binary_trace)
            for state_to_filter in [1, 0]:
                padded = np.pad(filtered_trace, 1, mode='constant', constant_values=1 - state_to_filter)
                diff = np.diff(padded)
                starts = np.where(diff == (1 if state_to_filter == 1 else -1))[0]
                ends = np.where(diff == (-1 if state_to_filter == 1 else 1))[0]
                for start, end in zip(starts, ends):
                    if end - start < min_length_points:
                        filtered_trace[start:end] = 1 - state_to_filter
            return filtered_trace

        self.results['signal_binarized_filtered'] = filter_short_events(signal_binarized, params['min_event_pts'])
        self.log("Binarized intensity trace based on threshold.")
        
        file_stem_step = f"{file_stem}-cycle#{params['kymo_cycle']}"
        trace_df = pd.DataFrame({'time': time_data, 'data': signal_data, 'fit': final_fit})
        trace_df.to_csv(os.path.join(output_dir_steps, f"{file_stem_step}-Intensity_along_DNAp-filtered_fit_trace.csv"), index=False)
        if not final_steps.empty:
            final_steps.to_csv(os.path.join(output_dir_steps, f"{file_stem_step}-Intensity_along_DNAp-filtered_step_properties.csv"), index=False)
        self.log(f"Step-finding results saved in '{output_dir_steps}'")

    def run_bp_changepoint(self, params, output_dir, trace_filename):
        df = self.data['trace']
        if df is None:
            self.log("No trace data loaded for BP changepoint detection.")
            self.results['bp_segments'] = None
            return

        folder_save = os.path.join(output_dir, 'ChangePoints_Results/')
        os.makedirs(folder_save, exist_ok=True)
        base_filename = os.path.splitext(os.path.basename(trace_filename))[0]

        data = np.zeros((len(df), 5))
        data[:, 0] = np.arange(1, len(df) + 1)
        data[:, 1] = df['time'].values
        data[:, 4] = df['basepairs'].values

        velocity_par, time_s = [], []
        nbpts = params['cpd_win_size']
        if len(data) > 2 * nbpts:
            for it in range(nbpts, len(data[:, 0]) - nbpts):
                win = np.arange(it - nbpts, it + nbpts + 1).astype(int)
                p = np.polyfit(data[win, 1], data[win, 4], 1)
                velocity_par.append(p[0])
                time_s.append(data[it, 1])
        
        if not velocity_par:
            self.log("Warning: Data too short for basepair change-point detection with current window size.")
            self.results['bp_segments'] = None
            return

        velocity_par_all = np.zeros((len(data[:, 0]), 2))
        velocity_par_all[:nbpts, 1] = velocity_par[0]
        velocity_par_all[nbpts:nbpts + len(time_s), 1] = velocity_par
        velocity_par_all[nbpts + len(time_s):, 1] = velocity_par[-1]
        
        segments, _ = bp_batch_segments_for([data[:, 0]], [velocity_par_all[:, 1]], sigma=params['cpd_sigma'], linear=False)
        
        self.results['velocity_par_all'] = velocity_par_all
        self.results['cp_data_array'] = data

        if segments and segments[0].size > 0:
            segments_array = segments[0]
            self.results['bp_segments'] = segments_array
            self.log(f"Basepair change-point detection found {len(segments_array)} segments.")
            
            bp_time, bp_signal = df['time'].values, df['basepairs'].values
            start_indices = np.clip((segments_array[:, 0] - 1).astype(int), 0, len(bp_time) - 1)
            end_indices = np.clip((segments_array[:, 2] - 1).astype(int), 0, len(bp_time) - 1)
            
            df_to_save = pd.DataFrame({
                'cp_startTime_s': bp_time[start_indices], 'cp_startBasepair': bp_signal[start_indices],
                'cp_endTime_s': bp_time[end_indices], 'cp_endBasepair': bp_signal[end_indices],
            })

            excel_filename = os.path.join(folder_save, f"{base_filename}-correlated_analysis.xlsx")
            with pd.ExcelWriter(excel_filename) as writer:
                df.to_excel(writer, sheet_name='raw_data', index=False)
                df_to_save.to_excel(writer, sheet_name='change_points', index=False)
            self.log(f"BP changepoint data saved to {excel_filename}")
        else:
            self.log("No basepair change-points were detected.")
            self.results['bp_segments'] = None


class AnalysisGUI(ctk.CTk):
    """
    An interactive GUI for correlating kymograph and force trace data.
    """
    def __init__(self):
        super().__init__()
        
        self.title("DNA Polymerase Trace & Force Correlation GUI")
        self.geometry("1600x900")
        ctk.set_appearance_mode("Dark")
        ctk.set_default_color_theme("blue")

        # --- Initialize state variables ---
        self.kymo_filename = ctk.StringVar()
        self.trace_filename = ctk.StringVar()
        default_output_dir = os.path.join(project_root, "example_dataset", "image data")
        self.output_dir = ctk.StringVar(value=default_output_dir)
        
        self.pipeline = AnalysisPipeline(logger_func=self.log)
        self.param_widgets = {}

        # --- Main layout ---
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)
        
        self.create_control_panel()
        self.create_plot_panel()
        self.create_log_panel()
    
    def create_control_panel(self):
        """Creates the left-side panel with all user controls and parameters."""
        control_frame = ctk.CTkFrame(self, width=400)
        control_frame.grid(row=0, column=0, rowspan=2, padx=10, pady=10, sticky="nswe")
        
        title_label = ctk.CTkLabel(control_frame, text="Analysis Controls", font=ctk.CTkFont(size=20, weight="bold"))
        title_label.pack(pady=10)

        self.tab_view = ctk.CTkTabview(control_frame, width=380)
        self.tab_view.pack(fill="both", expand=True, padx=5, pady=5)
        
        self.create_setup_tab()
        self.create_alignment_tab()
        self.create_intensity_tab()
        self.create_detection_tab()
        
        self.run_full_analysis_btn = ctk.CTkButton(control_frame, text="Run Full Analysis & Generate Report", command=self.run_full_analysis, state="disabled")
        self.run_full_analysis_btn.pack(side="bottom", fill="x", padx=10, pady=10)

    def create_setup_tab(self):
        tab = self.tab_view.add("Setup")
        
        # File Inputs
        ctk.CTkLabel(tab, text="Input Files").pack(pady=5, anchor="w")
        
        kymo_frame = ctk.CTkFrame(tab)
        ctk.CTkEntry(kymo_frame, textvariable=self.kymo_filename, width=250).pack(side="left", padx=5, expand=True, fill="x")
        ctk.CTkButton(kymo_frame, text="Kymo (.tdms)", command=lambda: self.select_file(self.kymo_filename)).pack(side="right")
        kymo_frame.pack(fill="x", padx=5, pady=2)
        
        trace_frame = ctk.CTkFrame(tab)
        ctk.CTkEntry(trace_frame, textvariable=self.trace_filename, width=250).pack(side="left", padx=5, expand=True, fill="x")
        ctk.CTkButton(trace_frame, text="Trace (.xlsx)", command=lambda: self.select_file(self.trace_filename, [("Excel files", "*.xlsx")])).pack(side="right")
        trace_frame.pack(fill="x", padx=5, pady=2)
        
        # Output Directory
        ctk.CTkLabel(tab, text="Output Directory").pack(pady=5, anchor="w")
        output_frame = ctk.CTkFrame(tab)
        ctk.CTkEntry(output_frame, textvariable=self.output_dir, width=250).pack(side="left", padx=5, expand=True, fill="x")
        ctk.CTkButton(output_frame, text="Select...", command=lambda: self.select_directory(self.output_dir)).pack(side="right")
        output_frame.pack(fill="x", padx=5, pady=2)

        # Kymo Cycle
        ctk.CTkLabel(tab, text="Kymo Cycle (e.g., #1)").pack(pady=5, anchor="w")
        self.kymo_cycle_entry = ctk.CTkEntry(tab)
        self.kymo_cycle_entry.insert(0, "#1")
        self.kymo_cycle_entry.pack(fill="x", padx=5, pady=2)
        self.param_widgets['kymo_cycle'] = self.kymo_cycle_entry
        
        # Parameter Management
        ctk.CTkLabel(tab, text="Parameter Management").pack(pady=(10,5), anchor="w")
        param_frame = ctk.CTkFrame(tab)
        param_frame.pack(fill="x", padx=5, pady=2)
        ctk.CTkButton(param_frame, text="Load Parameters", command=self.load_parameters).pack(side="left", padx=5, pady=5, expand=True, fill="x")
        ctk.CTkButton(param_frame, text="Save Parameters", command=self.save_parameters).pack(side="right", padx=5, pady=5, expand=True, fill="x")
        
        # Kymograph ROI
        ctk.CTkLabel(tab, text="Kymograph Display ROI (px)").pack(pady=(10,5), anchor="w")
        roi_frame = ctk.CTkFrame(tab)
        self.kymo_xlim_left = self.create_param_entry(roi_frame, "X Left", "5", 0, 0, 'kymo_xlim_left')
        self.kymo_xlim_right = self.create_param_entry(roi_frame, "X Right", "550", 0, 1, 'kymo_xlim_right')
        self.kymo_ylim_top = self.create_param_entry(roi_frame, "Y Top", "25", 1, 0, 'kymo_ylim_top')
        self.kymo_ylim_bottom = self.create_param_entry(roi_frame, "Y Bottom", "83", 1, 1, 'kymo_ylim_bottom')
        roi_frame.pack(fill="x", padx=5, pady=2)
        
        load_button = ctk.CTkButton(tab, text="Load Data & Visualize Kymo", command=self.load_and_visualize)
        load_button.pack(fill="x", padx=5, pady=10)

    def create_alignment_tab(self):
        tab = self.tab_view.add("Alignment")
        
        ctk.CTkLabel(tab, text="Initial Offsets (px)").pack(pady=5, anchor="w")
        offset_frame = ctk.CTkFrame(tab)
        self.x_offset_entry = self.create_param_entry(offset_frame, "X Offset", "-28", 0, 0, 'x_offset')
        self.y_offset_entry = self.create_param_entry(offset_frame, "Y Offset", "25", 0, 1, 'y_offset')
        offset_frame.pack(fill="x", padx=5, pady=2)

        ctk.CTkLabel(tab, text="Auto-Alignment Search Range (px)").pack(pady=5, anchor="w")
        search_frame = ctk.CTkFrame(tab)
        self.x_search_range = self.create_param_entry(search_frame, "X Range", "5", 0, 0, 'x_search')
        self.y_search_range = self.create_param_entry(search_frame, "Y Range", "2", 0, 1, 'y_search')
        search_frame.pack(fill="x", padx=5, pady=2)
        
        auto_align_button = ctk.CTkButton(tab, text="Find Optimal Alignment", command=self.run_auto_alignment)
        auto_align_button.pack(fill="x", padx=5, pady=10)

    def create_intensity_tab(self):
        tab = self.tab_view.add("Intensity")
        
        ctk.CTkLabel(tab, text="Signal Filtering").pack(pady=5, anchor="w", fill="x")
        sf_frame = ctk.CTkFrame(tab)
        sf_frame.pack(fill="x", padx=5, pady=2)
        self.signal_filter_window = self.create_param_entry(sf_frame, "Signal Filter Window", "7", 0, 0, 'signal_filter_window')

        ctk.CTkLabel(tab, text="Background ROI (px)").pack(pady=(10,5), anchor="w", fill="x")
        bg_frame = ctk.CTkFrame(tab)
        self.bg_rec_x = self.create_param_entry(bg_frame, "rec_x", "330", 0, 0, 'rec_x')
        self.bg_rec_y = self.create_param_entry(bg_frame, "rec_y", "50", 0, 1, 'rec_y')
        self.bg_rec_w = self.create_param_entry(bg_frame, "rec_w", "200", 1, 0, 'rec_w')
        self.bg_rec_h = self.create_param_entry(bg_frame, "rec_h", "5", 1, 1, 'rec_h')
        bg_frame.pack(fill="x", padx=5, pady=2)
        
        ctk.CTkButton(tab, text="Select BG ROI on Plot", command=self.select_roi_on_plot).pack(fill="x", padx=5, pady=5)

        ctk.CTkLabel(tab, text="Background Analysis").pack(pady=5, anchor="w")
        bg_analysis_frame = ctk.CTkFrame(tab)
        self.bg_filter_window = self.create_param_entry(bg_analysis_frame, "BG Filter Window", "25", 0, 0, 'bg_filter_window')
        self.threshold_sigma = self.create_param_entry(bg_analysis_frame, "Threshold (σ)", "3", 0, 1, 'threshold_sigma')
        bg_analysis_frame.pack(fill="x", padx=5, pady=2)

        self.analyze_intensity_btn = ctk.CTkButton(tab, text="Analyze Intensity & Background", command=self.run_intensity_analysis, state="disabled")
        self.analyze_intensity_btn.pack(fill="x", padx=5, pady=10)

    def create_detection_tab(self):
        tab = self.tab_view.add("Detection")

        ctk.CTkLabel(tab, text="Intensity Step-Finding (AutoStepFinder)").pack(pady=5, anchor="w")
        asf_frame = ctk.CTkFrame(tab)
        self.s_max_thresh = self.create_param_entry(asf_frame, "S_max Thresh", "0.15", 0, 0, 's_max_thresh')
        self.error_tol = self.create_param_entry(asf_frame, "Error Tolerance", "2.0", 0, 1, 'error_tol')
        self.min_event_pts = self.create_param_entry(asf_frame, "Min Event (pts)", "3", 1, 0, 'min_event_pts')
        asf_frame.pack(fill="x", padx=5, pady=2)
        
        self.run_stepfinder_btn = ctk.CTkButton(tab, text="Run Intensity Step-Finding", command=self.run_intensity_stepfinder_gui, state="disabled")
        self.run_stepfinder_btn.pack(fill="x", padx=5, pady=10)

        ctk.CTkLabel(tab, text="Basepair Change-Point Detection").pack(pady=5, anchor="w")
        cpd_frame = ctk.CTkFrame(tab)
        self.cpd_win_size = self.create_param_entry(cpd_frame, "Window Size", "9", 0, 0, 'cpd_win_size')
        self.cpd_sigma = self.create_param_entry(cpd_frame, "Sigma Noise", "0.08", 0, 1, 'cpd_sigma')
        cpd_frame.pack(fill="x", padx=5, pady=2)

        self.run_cpd_btn = ctk.CTkButton(tab, text="Run BP Change-Point Detection", command=self.run_bp_detection_gui, state="disabled")
        self.run_cpd_btn.pack(fill="x", padx=5, pady=10)

    def create_plot_panel(self):
        """Creates the right-side panel for displaying matplotlib plots."""
        plot_frame = ctk.CTkFrame(self)
        plot_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nswe")
        plot_frame.grid_rowconfigure(0, weight=1)
        plot_frame.grid_columnconfigure(0, weight=1)

        self.fig = plt.Figure(figsize=(8, 6), dpi=100)
        self.fig.set_tight_layout(True)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        
        toolbar_frame = ctk.CTkFrame(plot_frame)
        toolbar_frame.pack(side="bottom", fill="x")
        
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side="left", padx=10)

        save_plot_button = ctk.CTkButton(toolbar_frame, text="Save Current Plot", command=self.save_current_plot)
        save_plot_button.pack(side="right", padx=10, pady=5)

    def create_log_panel(self):
        """Creates the bottom panel for logging output."""
        log_frame = ctk.CTkFrame(self, height=150)
        log_frame.grid(row=1, column=1, padx=10, pady=10, sticky="nswe")
        
        self.log_textbox = ctk.CTkTextbox(log_frame, state="disabled", wrap="word")
        self.log_textbox.pack(fill="both", expand=True, padx=5, pady=5)

    # --- Helper & Utility Methods ---

    def log(self, message):
        """Appends a message to the log text box."""
        self.log_textbox.configure(state="normal")
        self.log_textbox.insert("end", str(message) + "\n")
        self.log_textbox.configure(state="disabled")
        self.log_textbox.see("end")
        self.update_idletasks() # Refresh GUI

    def select_file(self, var, filetypes=[("TDMS files", "*.tdms"), ("All files", "*.*")]):
        filename = filedialog.askopenfilename(filetypes=filetypes)
        if filename:
            var.set(filename)
    
    def select_directory(self, var):
        dirname = filedialog.askdirectory()
        if dirname:
            var.set(dirname)

    def create_param_entry(self, parent, label_text, default_value, row, col, param_key):
        frame = ctk.CTkFrame(parent)
        frame.grid(row=row, column=col, padx=5, pady=2, sticky="ew")
        ctk.CTkLabel(frame, text=label_text).pack(side="left", padx=(5,0))
        entry = ctk.CTkEntry(frame, width=70)
        entry.insert(0, default_value)
        entry.pack(side="right", padx=(0,5))
        self.param_widgets[param_key] = entry
        return entry

    def get_params(self):
        """Retrieves all parameters from the GUI, converting to correct types."""
        try:
            params = {
                'kymo_cycle': self.kymo_cycle_entry.get().replace("#", ""),
                'kymo_xlim_left': int(self.kymo_xlim_left.get()),
                'kymo_xlim_right': int(self.kymo_xlim_right.get()),
                'kymo_ylim_top': int(self.kymo_ylim_top.get()),
                'kymo_ylim_bottom': int(self.kymo_ylim_bottom.get()),
                'x_offset': int(self.x_offset_entry.get()),
                'y_offset': int(self.y_offset_entry.get()),
                'x_search': int(self.x_search_range.get()),
                'y_search': int(self.y_search_range.get()),
                'signal_filter_window': int(self.signal_filter_window.get()),
                'rec_x': int(self.bg_rec_x.get()),
                'rec_y': int(self.bg_rec_y.get()),
                'rec_w': int(self.bg_rec_w.get()),
                'rec_h': int(self.bg_rec_h.get()),
                'bg_filter_window': int(self.bg_filter_window.get()),
                'threshold_sigma': float(self.threshold_sigma.get()),
                's_max_thresh': float(self.s_max_thresh.get()),
                'error_tol': float(self.error_tol.get()),
                'min_event_pts': int(self.min_event_pts.get()),
                'cpd_win_size': int(self.cpd_win_size.get()),
                'cpd_sigma': float(self.cpd_sigma.get())
            }
            # Basic validation
            if params['signal_filter_window'] % 2 == 0 or params['signal_filter_window'] <= 3:
                raise ValueError("Signal filter window must be an odd number > 3.")
            if params['bg_filter_window'] % 2 == 0 or params['bg_filter_window'] <= 3:
                raise ValueError("Background filter window must be an odd number > 3.")
            return params
        except ValueError as e:
            messagebox.showerror("Parameter Error", f"Invalid parameter value: {e}. Please enter valid numbers.")
            return None

    def save_current_plot(self):
        """Saves the currently displayed plot to a file."""
        out_dir = self.output_dir.get()
        if not os.path.isdir(out_dir):
             messagebox.showwarning("Output Directory", "Please select a valid output directory first.")
             return
        filename = filedialog.asksaveasfilename(
            initialdir=out_dir,
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("PDF files", "*.pdf"), ("All files", "*.*")],
            title="Save Plot As..."
        )
        if not filename:
            return
        try:
            self.fig.savefig(filename, dpi=300, bbox_inches='tight')
            self.log(f"Plot saved to {filename}")
        except Exception as e:
            self.log(f"Error saving plot: {e}")
            messagebox.showerror("Save Error", f"Could not save the plot:\n{e}")
            
    def on_select_roi(self, eclick, erelease):
        """Callback for RectangleSelector."""
        x1, y1 = int(eclick.xdata), int(eclick.ydata)
        x2, y2 = int(erelease.xdata), int(erelease.ydata)
        
        # Update GUI entries
        self.param_widgets['rec_x'].delete(0, 'end'); self.param_widgets['rec_x'].insert(0, str(min(x1, x2)))
        self.param_widgets['rec_y'].delete(0, 'end'); self.param_widgets['rec_y'].insert(0, str(min(y1, y2)))
        self.param_widgets['rec_w'].delete(0, 'end'); self.param_widgets['rec_w'].insert(0, str(abs(x1 - x2)))
        self.param_widgets['rec_h'].delete(0, 'end'); self.param_widgets['rec_h'].insert(0, str(abs(y1 - y2)))
        
        self.rect_selector.set_active(False)
        self.canvas.draw_idle()
        self.log("Background ROI updated from plot selection.")

        # Re-plot to show the new rectangle
        params = self.get_params()
        if not params: return
        self.plot_kymograph_with_trace(rect={'x': params['rec_x'], 'y': params['rec_y'], 'w': params['rec_w'], 'h': params['rec_h']})

    def select_roi_on_plot(self):
        if 'g' not in self.pipeline.data:
            messagebox.showwarning("Plot Not Ready", "Please load data to display the kymograph first.")
            return
        
        ax = self.fig.get_axes()[0]
        self.log("Select a rectangle on the plot for the background ROI.")
        
        # The selector needs to be stored as an instance variable
        self.rect_selector = RectangleSelector(ax, self.on_select_roi, useblit=True, button=[1], minspanx=5, minspany=5, spancoords='data', interactive=True)
            
    # --- Core Analysis Logic (Refactored from original script) ---
    
    def load_and_visualize(self):
        """Loads TDMS and trace data, then shows the initial kymograph."""
        kymo_file = self.kymo_filename.get()
        trace_file = self.trace_filename.get()
        if not os.path.exists(kymo_file):
            messagebox.showerror("File Error", "Kymograph file not found.")
            return
        
        self.log(f"Loading kymograph: {kymo_file}")
        try:
            tif_path = self.pipeline.load_data(kymo_file, trace_file, self.output_dir.get())
            self.log(f"Kymograph saved as TIF: {tif_path}")

            # Enable next steps if data is loaded
            self.analyze_intensity_btn.configure(state="normal")
            self.run_full_analysis_btn.configure(state="normal")
            if self.pipeline.data.get('trace') is not None:
                self.run_cpd_btn.configure(state="normal")

        except Exception as e:
            self.log(f"Error loading data: {e}")
            messagebox.showerror("Data Loading Error", f"Failed to load data:\n{e}")
            return

        self.plot_kymograph_with_trace()

    def plot_kymograph_with_trace(self, trace_coords=None, trace_label="Force Trace", trace_color='yellow', rect=None):
        """Plots the base kymograph and optionally overlays a trace and rectangle."""
        params = self.get_params()
        if not params or 'g' not in self.pipeline.data:
            return

        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.imshow(self.pipeline.data['g'].astype('uint16'), cmap='gray', vmax=12, aspect="auto")
        
        if trace_coords and 'x' in trace_coords and 'y' in trace_coords:
            ax.plot(trace_coords['x'], trace_coords['y'], color=trace_color, linewidth=2, label=trace_label)
            ax.legend()
        
        if rect:
            patch = patches.Rectangle((rect['x'], rect['y']), rect['w'], rect['h'], linewidth=1, edgecolor='r', facecolor='none')
            ax.add_patch(patch)

        ax.set_xlabel("Time/px")
        ax.set_ylabel("Position/px")
        ax.set_title("DNA Polymerase on ssDNA")
        ax.set_xlim(params['kymo_xlim_left'], params['kymo_xlim_right'])
        ax.set_ylim(params['kymo_ylim_bottom'], params['kymo_ylim_top'])
        self.canvas.draw()
        
    def run_auto_alignment(self):
        """Performs the automatic alignment procedure."""
        params = self.get_params()
        if not all(k in self.pipeline.data for k in ['trace', 'g', 'time_per_line', 'px_size']):
            messagebox.showwarning("Data Missing", "Please load both kymo and trace data first.")
            return
        if self.pipeline.data['trace'] is None:
            messagebox.showwarning("Trace Data Missing", "Auto-alignment requires a trace file.")
            return

        try:
            x_offset_optimal, y_offset_optimal = self.pipeline.run_auto_alignment(params)
            
            self.x_offset_entry.delete(0, 'end'); self.x_offset_entry.insert(0, str(x_offset_optimal))
            self.y_offset_entry.delete(0, 'end'); self.y_offset_entry.insert(0, str(y_offset_optimal))

            # Re-plot with new alignment
            x_cali = 1000 / self.pipeline.data['time_per_line']
            y_cali = 1000 / self.pipeline.data['px_size']
            trace = self.pipeline.data['trace']
            trace_time_opt = trace['time'] / 1000 * x_cali + x_offset_optimal
            pos_opt = pd.to_numeric(trace['junction_position_all'], errors='coerce') * y_cali + y_offset_optimal
            self.plot_kymograph_with_trace(
                trace_coords={'x': trace_time_opt, 'y': pos_opt},
                trace_label="Optimal Alignment",
                trace_color='cyan'
            )
        except RuntimeError as e:
            self.log(str(e))
            messagebox.showerror("Alignment Error", str(e))
        except Exception as e:
            self.log(f"An unexpected error occurred during alignment: {e}")
            messagebox.showerror("Alignment Error", f"An unexpected error occurred: {e}")

    def run_full_analysis(self):
        """Executes the entire analysis pipeline from start to finish."""
        self.log("\n" + "="*40)
        self.log("--- STARTING FULL ANALYSIS PIPELINE ---")
        
        # 1. Check inputs, get parameters, set up paths
        params = self.get_params()
        if not params: return
        
        out_dir = self.output_dir.get()
        os.makedirs(out_dir, exist_ok=True)
        self.log(f"Results will be saved in: '{out_dir}'")
        file_stem = os.path.splitext(os.path.basename(self.kymo_filename.get()))[0]
        
        if 'g' not in self.pipeline.data:
            self.load_and_visualize()
            if 'g' not in self.pipeline.data:
                self.log("Analysis cancelled. Please load kymograph data.")
                return

        try:
            # --- STAGE 1: INTENSITY ANALYSIS ---
            if not self.run_intensity_analysis(silent=True): return
            
            # --- STAGE 2: STEP FINDING ---
            if not self.run_intensity_stepfinder_gui(silent=True): return

            # --- STAGE 3: BP CHANGE POINT ---
            if self.pipeline.data.get('trace') is not None:
                if not self.run_bp_detection_gui(silent=True): return

            # --- STAGE 4: FINAL PLOTS & DATA SAVING ---
            self.log("\n[Final Stage] Generating final plots and saving all results...")
            self._generate_final_plots(params, file_stem, out_dir)
            self._save_consolidated_excel(params, file_stem, out_dir)
            self._save_analysis_log(params, file_stem, out_dir)

            self.log("\n--- ANALYSIS COMPLETE ---")
            messagebox.showinfo("Success", "Full analysis pipeline completed successfully!")

        except Exception as e:
            self.log(f"\n--- ANALYSIS FAILED ---")
            self.log(f"An error occurred: {e}")
            import traceback
            self.log(traceback.format_exc())
            messagebox.showerror("Analysis Error", f"The analysis failed with the following error:\n\n{e}")

    def run_intensity_analysis(self, silent=False):
        """Runs Stage 1: Intensity trace extraction and background analysis."""
        if not silent: self.log("\n--- Running: Intensity & Background Analysis ---")
        params = self.get_params()
        if not params: return False
        
        if self.pipeline.data.get('trace') is None:
            if not silent: messagebox.showwarning("Trace Data Missing", "Intensity analysis requires a trace file.")
            return False
        
        try:
            self.pipeline.calculate_intensity_trace(params)
            self.pipeline.analyze_background(params)
            if not silent:
                self._plot_intensity_trace(params)
                self.log("--- Intensity & Background Analysis Complete ---")
            
            self.run_stepfinder_btn.configure(state="normal") # Enable next step
            return True
        except Exception as e:
            self.log(f"Error during intensity analysis: {e}")
            if not silent: messagebox.showerror("Error", f"Intensity analysis failed:\n{e}")
            return False

    def run_intensity_stepfinder_gui(self, silent=False):
        """Runs intensity step-finding and plots the result."""
        if 'filtered_intensity' not in self.pipeline.results:
            messagebox.showwarning("Prerequisite Missing", "Please run 'Analyze Intensity & Background' first.")
            return False
            
        if not silent: self.log("\n--- Running: Intensity Step-Finding ---")
        params = self.get_params()
        if not params: return False

        try:
            file_stem = os.path.splitext(os.path.basename(self.kymo_filename.get()))[0]
            self.pipeline.run_intensity_stepfinder(params, self.output_dir.get(), file_stem)
            if not silent:
                self._plot_detection_results(params)
                self.log("--- Intensity Step-Finding Complete ---")
            return True
        except Exception as e:
            self.log(f"Error during step-finding: {e}")
            if not silent: messagebox.showerror("Error", f"Step-finding failed:\n{e}")
            return False

    def run_bp_detection_gui(self, silent=False):
        """Runs basepair change-point detection and plots the result."""
        if self.pipeline.data.get('trace') is None:
            messagebox.showwarning("Prerequisite Missing", "Please load a trace file for this analysis.")
            return False
        
        if not silent: self.log("\n--- Running: Basepair Change-Point Detection ---")
        params = self.get_params()
        if not params: return False
        
        try:
            self.pipeline.run_bp_changepoint(params, self.output_dir.get(), self.trace_filename.get())
            if not silent:
                self._plot_bp_changepoints()
                self.log("--- Basepair Change-Point Detection Complete ---")
            return True
        except Exception as e:
            self.log(f"Error during BP detection: {e}")
            if not silent: messagebox.showerror("Error", f"BP detection failed:\n{e}")
            return False

    def _plot_intensity_trace(self, params):
        """Generates a plot for the intensity analysis step."""
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        
        time_s = self.pipeline.results['time_s']
        raw_intensity = self.pipeline.results['raw_intensity']
        filtered_intensity = self.pipeline.results['filtered_intensity']
        bg_filter = self.pipeline.results['bagrnd_signal_filter']
        threshold = np.average(bg_filter) + params['threshold_sigma'] * np.std(bg_filter)
        
        ax.plot(time_s, raw_intensity, color='lightgray', label='Raw Intensity')
        ax.plot(time_s, filtered_intensity, color='cornflowerblue', label='Filtered Intensity')
        ax.axhline(threshold, color='red', linestyle='--', label=f'{params["threshold_sigma"]}σ Threshold')
        
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Intensity (a.u.)")
        ax.set_title("Intensity Along DNAp Trajectory")
        ax.legend()
        ax.set_xlim(time_s[0], time_s[-1])
        self.canvas.draw()
        self.log("Plotted intensity trace results.")

    def _plot_detection_results(self, params):
        """Generates a plot for the detection step."""
        self.fig.clear()
        ax = self.fig.add_subplot(111)

        time_s = self.pipeline.results['time_s']
        filtered_intensity = self.pipeline.results['filtered_intensity']
        final_fit = self.pipeline.results['final_fit']
        
        ax.plot(time_s, filtered_intensity, color='lightgray', label='Filtered Intensity')
        ax.plot(time_s, final_fit, color='black', linewidth=1.5, label='Step Fit')
        
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Intensity (a.u.)")
        ax.set_title("Intensity Step-Finding Results")
        ax.legend()
        ax.set_xlim(time_s[0], time_s[-1])
        self.canvas.draw()
        self.log("Plotted step-detection results.")

    def _plot_bp_changepoints(self):
        self.fig.clear()
        ax = self.fig.add_subplot(111)

        if 'trace' not in self.pipeline.data or self.pipeline.data['trace'] is None:
            self.log("No trace data to plot.")
            return

        df_trace = self.pipeline.data['trace']
        bp_time_s = df_trace['time'].values / 1000 # Convert to seconds
        bp_signal = df_trace['basepairs'].values
        
        ax.plot(bp_time_s, bp_signal, label='Basepairs Data', color='green', alpha=0.7)
        
        if self.pipeline.results.get('bp_segments') is not None:
            segments = self.pipeline.results['bp_segments']
            start_indices = (segments[:, 0] - 1).astype(int)
            end_indices = (segments[:, 2] - 1).astype(int)
            start_indices = np.clip(start_indices, 0, len(bp_time_s) - 1)
            end_indices = np.clip(end_indices, 0, len(bp_time_s) - 1)

            for i in range(len(start_indices)):
                segment_time = bp_time_s[start_indices[i]:end_indices[i]+1]
                segment_bp = bp_signal[start_indices[i]:end_indices[i]+1]
                
                if len(segment_time) > 1:
                    p = np.polyfit(segment_time, segment_bp, 1)
                    fit_line = np.polyval(p, segment_time)
                    ax.plot(segment_time, fit_line, color='red', linewidth=2.5, label='Segment Fit' if i==0 else "")

        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Basepairs")
        ax.set_title("Basepair Change-Point Detection")
        ax.legend()
        self.canvas.draw()
        self.log("Plotted basepair change-point results.")

    def _generate_final_plots(self, params, file_stem, output_dir):
        cycle_str = f"cycle#{params['kymo_cycle']}"
        self.log("Generating and saving final plots...")
        results_dir = os.path.join(output_dir, "results")
        os.makedirs(results_dir, exist_ok=True)

        # Plot 1: Kymograph with final trace overlay
        if 'trace' in self.pipeline.data and self.pipeline.data['trace'] is not None and 'trace_time_px' in self.pipeline.results:
            self.plot_kymograph_with_trace(
                trace_coords={'x': self.pipeline.results['trace_time_px'], 'y': self.pipeline.results['position_px']},
                trace_label="Optimal Alignment", trace_color='lime',
                rect={'x': params['rec_x'], 'y': params['rec_y'], 'w': params['rec_w'], 'h': params['rec_h']}
            )
            plot_path = os.path.join(results_dir, f"{file_stem}-{cycle_str}-overlap-optimizing.png")
            self.fig.savefig(plot_path, dpi=300, bbox_inches='tight')
            self.log(f"Saved optimized overlap plot: {plot_path}")
        
        # Plot: Background ROI on Kymo
        self.plot_kymograph_with_trace(
            rect={'x': params['rec_x'], 'y': params['rec_y'], 'w': params['rec_w'], 'h': params['rec_h']}
        )
        plot_path = os.path.join(results_dir, f"{file_stem}-{cycle_str}-background Intensity.png")
        self.fig.savefig(plot_path, dpi=300, bbox_inches='tight')
        self.log(f"Saved background ROI plot: {plot_path}")
        
        # Plot 2: Raw Bp-intensity overlay
        if 'trace' in self.pipeline.data and self.pipeline.data['trace'] is not None:
            fig_bp_int = plt.figure(figsize=(8,3))
            ax1 = fig_bp_int.add_subplot(111)
            ax2 = ax1.twinx()
            bp_time = self.pipeline.data['trace']['time']/1000
            bp = self.pipeline.data['trace']['basepairs']
            x1 = self.pipeline.results['time_s']
            all_intensity = self.pipeline.results['raw_intensity']
            ax1.plot(bp_time,bp,color='black',linewidth=1, label='basepairs')
            ax1.set_xlabel('Time/s')
            ax1.set_ylabel('Basepairs', color='black')
            ax2.plot(x1, all_intensity, color='green', linewidth=0.2, label='fluorescence intensity')
            ax2.set_ylabel('Fluorescence intensity(a.u.)', color='green')
            ax1.set_xlim((params['kymo_xlim_left']-params['x_offset'])/self.pipeline.results['x_cali'], (params['kymo_xlim_right']-params['x_offset'])/self.pipeline.results['x_cali'])
            fig_bp_int.tight_layout()
            plot_path = os.path.join(results_dir, f"{file_stem}-{cycle_str}-raw Bp-intensity along DNAp Trajectory.png")
            fig_bp_int.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close(fig_bp_int)
            self.log(f"Saved raw Bp-Intensity plot: {plot_path}")
        
        # Plot: Background Intensity Trace
        if 'bagrnd_raw_intensity' in self.pipeline.results:
            fig_bg, ax_bg = plt.subplots(figsize=(8,3))
            ax_bg.plot(self.pipeline.results['bagrnd_time_px'], self.pipeline.results['bagrnd_raw_intensity'], linestyle='solid', linewidth=0.3, label='Background Signal')
            ax_bg.plot(self.pipeline.results['bagrnd_time_px'], self.pipeline.results['bagrnd_signal_filter'], linestyle='solid', linewidth=1, label='Background Signal Filter')
            threshold_line = np.average(self.pipeline.results['bagrnd_signal_filter']) + params['threshold_sigma'] * np.std(self.pipeline.results['bagrnd_signal_filter'])
            ax_bg.axhline(threshold_line, color='k', linestyle='--',label= f"{params['threshold_sigma']} Sigma&filterSize {params['bg_filter_window']}")
            ax_bg.set_xlabel("Time/px")
            ax_bg.set_ylabel("Intensity/photon")
            ax_bg.set_title("Background Intensity")
            ax_bg.legend()
            fig_bg.tight_layout()
            plot_path = os.path.join(results_dir, f"{file_stem}-{cycle_str}-background Intensity along DNA polymerase Trajectory.png")
            fig_bg.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close(fig_bg)
            self.log(f"Saved background intensity trace plot: {plot_path}")

        # Plot: Filtered Intensity Trace with Threshold
        if 'filtered_intensity' in self.pipeline.results:
            fig_int, ax_int = plt.subplots(figsize=(8, 3))
            time_s = self.pipeline.results['time_s']
            raw_intensity = self.pipeline.results['raw_intensity']
            filtered_intensity = self.pipeline.results['filtered_intensity']
            bg_filter = self.pipeline.results['bagrnd_signal_filter']
            threshold = np.average(bg_filter) + params['threshold_sigma'] * np.std(bg_filter)
            ax_int.plot(time_s, raw_intensity, linestyle='solid', linewidth=0.3, label='Protein Trace Signal')
            ax_int.plot(time_s, filtered_intensity, linestyle='solid', linewidth=1, label='Protein Trace Signal Filter')
            ax_int.axhline(threshold, color='grey', linestyle='--', label=f"{params['threshold_sigma']} Sigma threshold")
            ax_int.set_xlabel("Time/s")
            ax_int.set_ylabel("Intensity/photon")
            ax_int.set_title("Intensity along DNA polymerase Trajectory")
            ax_int.legend()
            ax_int.set_xlim((params['kymo_xlim_left']-params['x_offset'])/self.pipeline.results['x_cali'], (params['kymo_xlim_right']-params['x_offset'])/self.pipeline.results['x_cali'])
            fig_int.tight_layout()
            plot_path = os.path.join(results_dir, f"{file_stem}-{cycle_str}-Intensity along DNA polymerase Trajectory-1.png")
            fig_int.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close(fig_int)
            self.log(f"Saved filtered intensity trace plot: {plot_path}")

        # Plot: Triple-axis final plot
        if 'trace' not in self.pipeline.data or self.pipeline.data['trace'] is None or 'final_fit' not in self.pipeline.results:
            self.log("Skipping final correlated plot as required data (trace or intensity fit) is not available.")
        else:
            fig_triple = plt.figure(figsize=(10, 5))
            host = host_subplot(111, axes_class=AA.Axes)
            plt.subplots_adjust(right=0.7)
            par1 = host.twinx()
            par2 = host.twinx()
            offset = 40
            new_fixed_axis = par2.get_grid_helper().new_fixed_axis
            par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))
            par2.axis["right"].toggle(all=True)

            bp_time = self.pipeline.data['trace']['time'] / 1000
            bp = self.pipeline.data['trace']['basepairs']
            x1 = self.pipeline.results['time_s']
            final_fit = self.pipeline.results['final_fit']
            binarized = self.pipeline.results['signal_binarized_filtered']
            threshold = np.mean(self.pipeline.results['bagrnd_signal_filter']) + (params['threshold_sigma'] * np.std(self.pipeline.results['bagrnd_signal_filter']))
            
            host.set_xlabel("Time (s)")
            host.set_ylabel("Basepairs (bp)")
            par1.set_ylabel("Intensity (a.u.)")
            par2.set_ylabel("Binarized State")
            
            p0, = host.plot(bp_time, bp, "green", label='Basepairs', linewidth=1.5)
            p1, = par1.plot(x1, final_fit, color='black', linewidth=1.5, label='Step-Fit Intensity')
            par1.plot(x1, self.pipeline.results['raw_intensity'], color='lightgray', linewidth=0.8, label='Raw Intensity')
            par1.axhline(threshold, color='grey', linestyle='--', label=f'{params["threshold_sigma"]}σ Threshold')
            p2 = par2.fill_between(x1, binarized, 0, color='orange', alpha=0.3, label='DNAp Bound State')
            
            host.set_xlim((params['kymo_xlim_left']-params['x_offset'])/self.pipeline.results['x_cali'], (params['kymo_xlim_right']-params['x_offset'])/self.pipeline.results['x_cali'])
            host.axis["left"].label.set_color(p0.get_color())
            par1.axis["right"].label.set_color(p1.get_color())
            par2.axis["right"].label.set_color('orange')
            host.legend(handles=[p0, p1, p2], loc='best')

            correlated_plot_path = os.path.join(output_dir, 'stepfinding_results', f"{file_stem}-{cycle_str}-Intensity_along_DNAp-filtered_all_correlated.png")
            fig_triple.savefig(correlated_plot_path, dpi=300, bbox_inches='tight')
            self.log(f"Generated and saved final correlated plot to: {correlated_plot_path}")
            plt.close(fig_triple) # Close to free memory

    def _save_consolidated_excel(self, params, file_stem, output_dir):
        if self.pipeline.results.get('bp_segments') is None:
            self.log("Skipping Excel report generation as no BP segments were found.")
            return

        base_filename = os.path.splitext(os.path.basename(self.trace_filename.get()))[0]
        folder_save = os.path.join(output_dir, 'ChangePoints_Results/')
        excel_filename = os.path.join(folder_save, f"{base_filename}-correlated_analysis.xlsx")

        bp_time = self.pipeline.data['trace']['time'].values
        bp_signal = self.pipeline.data['trace']['basepairs'].values
        segments_array = self.pipeline.results['bp_segments']

        start_indices = (segments_array[:, 0] - 1).astype(int)
        end_indices = (segments_array[:, 2] - 1).astype(int)
        cp_startTime = bp_time[start_indices]
        cp_endTime = bp_time[end_indices]
        cp_startBasepair = bp_signal[start_indices]
        cp_endBasepair = bp_signal[end_indices]

        with pd.ExcelWriter(excel_filename) as writer:
            df1 = pd.DataFrame({'time_s': bp_time, 'raw_basepair': bp_signal})
            df1.to_excel(writer, sheet_name='raw_data', index=False)
            
            df2 = pd.DataFrame({
                'cp_startTime_s': cp_startTime, 'cp_startBasepair': cp_startBasepair,
                'cp_endTime_s': cp_endTime, 'cp_endBasepair': cp_endBasepair,
            })
            df2.to_excel(writer, sheet_name='change_points', index=False)
            
            cp_Time = np.append(cp_startTime, cp_endTime[-1])
            cp_Basepair = np.append(cp_startBasepair, cp_endBasepair[-1])
            interp_df = pd.DataFrame({"Time": cp_Time, "Basepair": cp_Basepair}).drop_duplicates()
            set_interp = interp1d(interp_df['Time'], interp_df['Basepair'], kind='linear', fill_value="extrapolate")
            cp_basepair_interp = set_interp(bp_time)
            df3 = pd.DataFrame({'time_s': bp_time, 'interpolated_basepair': cp_basepair_interp})
            df3.to_excel(writer, sheet_name='interpolated_data', index=False)

            df4 = pd.DataFrame({'time_s': self.pipeline.results['time_s'], 'step_intensity': self.pipeline.results['final_fit']})
            df4.to_excel(writer, sheet_name='step_intensity', index=False)
    
            df5 = pd.DataFrame({'time_s': self.pipeline.results['time_s'], 'binarized_intensity': self.pipeline.results['signal_binarized_filtered']})
            df5.to_excel(writer, sheet_name='binarized_intensity', index=False)
        
        self.log(f"Combined analysis data saved to {excel_filename}")

    def _save_analysis_log(self, params, file_stem, output_dir):
        base_filename = os.path.splitext(os.path.basename(self.trace_filename.get()))[0]
        folder_save = os.path.join(output_dir, 'ChangePoints_Results/')
        excel_filename = os.path.join(folder_save, f"{base_filename}-correlated_analysis.xlsx")
        log_path = excel_filename.replace('.xlsx', '-analysis_log_data.txt')

        log_data = {
            "kymo_cycle": params['kymo_cycle'],
            "kymo_xlim_left": params['kymo_xlim_left'],
            "kymo_xlim_right": params['kymo_xlim_right'],
            "x_offset": params['x_offset'],
            "y_offset": params['y_offset'],
            "rec_x": params['rec_x'],
            "rec_y": params['rec_y'],
            "rec_w": params['rec_w'],
            "rec_h": params['rec_h'],
            "bagrnd_filter_window": params['bg_filter_window'],
            "signal_filter_window": params['signal_filter_window'],
            "threshold_sigma": params['threshold_sigma'],
            "Intensity_changepoint#": len(self.pipeline.results.get('final_steps', [])),
            "windowsize_cp_basepair": params['cpd_win_size'],
            "sigmavalue_basepair": params['cpd_sigma']
        }
        
        with open(log_path, 'w') as f:
            for key, value in log_data.items():
                f.write(f"{key}: {value}\n")
        self.log(f"Analysis log saved to: {log_path}")

    def save_parameters(self):
        params_to_save = {key: widget.get() for key, widget in self.param_widgets.items()}
        filename = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            title="Save Parameters As..."
        )
        if not filename: return
        try:
            with open(filename, 'w') as f:
                json.dump(params_to_save, f, indent=4)
            self.log(f"Parameters saved to {filename}")
        except Exception as e:
            messagebox.showerror("Save Error", f"Could not save parameters:\n{e}")

    def load_parameters(self):
        filename = filedialog.askopenfilename(
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            title="Load Parameters"
        )
        if not filename: return
        try:
            with open(filename, 'r') as f:
                loaded_params = json.load(f)
            for key, value in loaded_params.items():
                if key in self.param_widgets:
                    self.param_widgets[key].delete(0, 'end')
                    self.param_widgets[key].insert(0, str(value))
            self.log(f"Parameters loaded from {filename}")
        except Exception as e:
            messagebox.showerror("Load Error", f"Could not load parameters:\n{e}")

if __name__ == "__main__":
    app = AnalysisGUI()
    app.mainloop()

