#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cv2
from matplotlib.colors import Normalize, LinearSegmentedColormap
import seaborn as sns
from scipy.interpolate import interp1d
from scipy import interpolate
from more_itertools import chunked
import argparse
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import json


def load_data(cor_file, cor_file2, cor_file3):
    """Loads all necessary data files for the analysis."""
    raw_basepair = pd.read_excel(cor_file, sheet_name='raw_data', engine='openpyxl')
    cp_basepair = pd.read_excel(cor_file, sheet_name='change_points', engine='openpyxl')
    raw_intensity = np.loadtxt(cor_file2)
    filtered_intensity = np.genfromtxt(cor_file3, delimiter=',', skip_header=1)
    step_intensity = pd.read_excel(cor_file, sheet_name='step_intensity', engine='openpyxl')
    binarized_intensity = pd.read_excel(cor_file, sheet_name='binarized_intensity', engine='openpyxl')

    data = {
        'time_intens': binarized_intensity['time_s'].to_numpy(),
        'intensity_raw': raw_intensity.T[1],
        'intensity_filtered': filtered_intensity.T[1],
        'intensity': binarized_intensity['binarized_intensity'].to_numpy(),
        'intensity_step': step_intensity['step_intensity'].to_numpy(),
        'bp_time': raw_basepair['time_s'].to_numpy()/1000,
        'bp': raw_basepair['raw_basepair'].to_numpy(),
        'cp_basepair': cp_basepair,
        'cp_basepair_interp': pd.read_excel(cor_file, sheet_name='interpolated_data', engine='openpyxl')['interpolated_basepair'].to_numpy()
    }
    return data


def plot_correlation_overview(data, output_dir, name_prefix, show_plots=True):
    """Plots and saves an overview of the correlated force and image data."""
    font = {'family': 'Arial', 'weight': 'normal', 'size': 16}

    # Plot 1: Raw basepairs
    fig1 = plt.figure(figsize=(8, 3))
    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)
    par1 = host.twinx()
    par2 = host.twinx()
    offset = 25
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))
    par2.axis["right"].toggle(all=True)

    host.set_xlabel("Time/s", fontdict=font)
    host.set_ylabel("Basepair/bp", fontdict=font)
    par1.set_ylabel("Step-like Intensity", fontdict=font)
    par2.set_ylabel("Binarized Intensity", fontdict=font)

    host.plot(data['bp_time'], data['bp'], "lime", label='Basepairs', linewidth=1)
    p1, = par1.plot(data['time_intens'], data['intensity_step'], color='black', linewidth=1, label='Step-like intensity')
    par2.plot(data['time_intens'], data['intensity'], color='lightgray', linewidth=0.1, label='Binarized Intensity')
    par2.fill_between(data['time_intens'], data['intensity'], 0, color='lightgreen', alpha=0.1)
    par1.plot(data['time_intens'], data['intensity_raw'], color='lightgray', linewidth=0.6, linestyle='solid', markersize=0.5, label='Raw Intensity')

    host.autoscale()
    host.margins(0.1)
    par1.set_ylim(10, 77)
    par2.set_ylim(0, 1)
    par2.axis["right"].label.set_color(p1.get_color())
    par2.axis('off')
    par1.yaxis.set_ticklabels([])
    par1.yaxis.set_ticks([])
    plt.tight_layout()
    figure_path = os.path.join(output_dir, f"{name_prefix}-all_correlated_data-replot.eps")
    plt.savefig(figure_path, format='eps', dpi=300, bbox_inches='tight')
    print(f"Replot figure saved to: {figure_path}")
    if show_plots: plt.show()
    plt.close(fig1)

    # Plot 2: Segmented basepairs by rate
    fig2 = plt.figure(figsize=(8, 3))
    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)
    par1 = host.twinx()
    par2 = host.twinx()
    offset = 25
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))
    par2.axis["right"].toggle(all=True)

    host.set_xlabel("Time/s", fontdict=font)
    host.set_ylabel("Basepair/bp", fontdict=font)
    par1.set_ylabel("Step-like Intensity", fontdict=font)
    par2.set_ylabel("Binarized Intensity", fontdict=font)
    host.plot(data['bp_time'], data['bp'], color="lightgray", label='_nolegend_', linewidth=1)
    for index, row in data['cp_basepair'].iterrows():
        host.plot([row['cp_startTime_s']/1000, row['cp_endTime_s']/1000],
                  [row['cp_startBasepair'], row['cp_endBasepair']],
                  linewidth=1.5)

    p1, = par1.plot(data['time_intens'], data['intensity_step'], color='black', linewidth=1, label='Step-like intensity')
    par2.plot(data['time_intens'], data['intensity'], color='lightgray', linewidth=0.1, label='Binarized Intensity')
    par2.fill_between(data['time_intens'], data['intensity'], 0, color='lightgreen', alpha=0.1)
    par1.plot(data['time_intens'], data['intensity_raw'], color='lightgray', linewidth=0.6, linestyle='solid', markersize=0.5, label='Raw Intensity')

    host.autoscale()
    host.margins(0.1)
    par1.set_ylim(10, 77)
    par2.set_ylim(0, 1)
    par2.axis["right"].label.set_color(p1.get_color())
    par2.axis('off')
    par1.yaxis.set_ticklabels([])
    par1.yaxis.set_ticks([])
    plt.tight_layout()
    figure_path = os.path.join(output_dir, f"{name_prefix}-all_correlated_data-replot_ColoredRates.eps")
    plt.savefig(figure_path, format='eps', dpi=300, bbox_inches='tight')
    print(f"Replot figure saved to: {figure_path}")
    if show_plots: plt.show()
    plt.close(fig2)


def plot_activity_burst(data, output_dir, name_prefix, show_plots=True):
    """Plots the DNAp activity burst rate against binarized intensity."""
    fig = plt.figure(figsize=(8,3))
    plt.rcParams.update({'font.size': 14, 'axes.labelsize': 14, 'xtick.labelsize': 12, 'ytick.labelsize': 12})

    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)
    par1 = host.twinx()
    par2 = host.twinx()
    offset = 25
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))
    par2.axis["right"].toggle(all=True)

    host.set_xlabel("Time/s")
    host.set_ylabel("DNAp Activity Burst(bp/s)")
    par2.set_ylabel("Binarized Intensity")

    p0, = host.plot(data['bp_time'][:-1], np.diff(data['cp_basepair_interp'])/np.diff(data['bp_time']), "red", label="DNAp Activity Burst", linewidth=1)
    p2, = par2.plot(data['time_intens'], data['intensity'], color='lightgray', linewidth=0.1, label='Binarized Intensity')
    par2.fill_between(data['time_intens'], data['intensity'], 0, color='lightgreen', alpha=0.1)

    host.autoscale()
    host.margins(0.1)
    par2.set_ylim(0, 1)
    par2.axis["right"].label.set_color(p0.get_color())
    par2.axis('off')
    par1.yaxis.set_ticklabels([])
    par1.yaxis.set_ticks([])

    plt.tight_layout()
    figure_path = os.path.join(output_dir, f"{name_prefix}-all_correlated_data-replot_Rate.eps")
    plt.savefig(figure_path, format='eps', dpi=300, bbox_inches='tight')
    print(f"Replot_rate figure saved to: {figure_path}")
    if show_plots: plt.show()
    plt.close(fig)


def find_binding_events(intensity, time_intens, filter_dots=3):
    """Finds binding events from binarized intensity data."""
    index = np.where(intensity == 0)[0]
    index = np.insert(index, 0, -1)
    index = np.append(index, len(intensity))
    
    time_list = [time_intens[index[i]+1:index[i+1]] for i in range(len(index)-1)]
    
    time_lst = [x for x in time_list if x.size >= filter_dots]
    
    start_time = [arr[0] for arr in time_lst]
    end_time = [arr[-1] for arr in time_lst]

    return start_time, end_time


def distinguish_exo_pol_events(start_time, end_time, time_intens, intensity, force_transi_time):
    """Distinguishes between exonuclease and polymerase events based on force transition time."""
    func = interp1d(time_intens, intensity, kind='nearest', fill_value="extrapolate")
    force_transi_intensity = int(func(force_transi_time))
    
    start_time_exo, end_time_exo = [], []
    start_time_pol, end_time_pol = [], []

    if force_transi_intensity == 0:
        for i, t_end in enumerate(end_time):
            if t_end < force_transi_time:
                start_time_exo.append(start_time[i])
                end_time_exo.append(t_end)
        for i, t_start in enumerate(start_time):
            if t_start > force_transi_time:
                start_time_pol.append(t_start)
                end_time_pol.append(end_time[i])
    else: # force_transi_intensity == 1
        crossing_event_idx = -1
        for i, (t_start, t_end) in enumerate(zip(start_time, end_time)):
            if t_start <= force_transi_time < t_end:
                crossing_event_idx = i
                break

        for i, t_end in enumerate(end_time):
            if i < crossing_event_idx:
                start_time_exo.append(start_time[i])
                end_time_exo.append(t_end)

        if crossing_event_idx != -1:
            start_time_exo.append(start_time[crossing_event_idx])
            end_time_exo.append(force_transi_time)
            start_time_pol.append(force_transi_time)
            end_time_pol.append(end_time[crossing_event_idx])
        
        for i, t_start in enumerate(start_time):
             if i > crossing_event_idx:
                start_time_pol.append(t_start)
                end_time_pol.append(end_time[i])
                
    return start_time_exo, end_time_exo, start_time_pol, end_time_pol


def find_dark_events(start_times, end_times, boundary_start, boundary_end, first_interval_start_time):
    """Identifies dark (unbound) periods between bound events."""
    dark_start, dark_end = [], []
    if not start_times:
        if boundary_start < boundary_end:
            dark_start.append(boundary_start)
            dark_end.append(boundary_end)
        return dark_start, dark_end

    if start_times[0] > boundary_start:
        dark_start.append(boundary_start)
        dark_end.append(start_times[0])

    for i in range(len(end_times) - 1):
        dark_start.append(end_times[i])
        dark_end.append(start_times[i+1])
    
    if end_times[-1] < boundary_end:
        dark_start.append(end_times[-1])
        dark_end.append(boundary_end)

    return dark_start, dark_end


def calculate_segment_rates(time, basepairs, start_times, end_times):
    """Calculates rates (bp/s) for given time segments."""
    rates = []
    for start, end in zip(start_times, end_times):
        indtemp = np.where((time <= end) & (time >= start))
        time_ROI = time[indtemp]
        basepairs_ROI = basepairs[indtemp]
        if len(time_ROI) > 1:
            duration_time = time_ROI[-1] - time_ROI[0]
            processivity = basepairs_ROI[-1] - basepairs_ROI[0]
            if duration_time > 0:
                rates.append(processivity / duration_time)
            else:
                rates.append(0.0) # or np.nan
        else:
            rates.append(0.0) # or np.nan
    return np.array(rates)


def save_segmented_data(output_dir, name_prefix, analysis_results):
    """Saves all segmented track data to an Excel file."""
    excel_path = os.path.join(output_dir, f"{name_prefix}-SegmentedTrack.xlsx")

    def create_df_safely(data_dict, sheet_name, writer):
        try:
            it = iter(data_dict.values())
            the_len = len(next(it))
            if not all(len(l) == the_len for l in it):
                print(f"Warning: Mismatched lengths in data for sheet '{sheet_name}'. Skipping.")
                return
            if the_len > 0:
                df = pd.DataFrame(data_dict)
                df.to_excel(writer, sheet_name=sheet_name, index=False)
            else:
                print(f"Info: No data to write for sheet '{sheet_name}'.")
        except (StopIteration, TypeError):
             print(f"Info: No data to write for sheet '{sheet_name}'.")

    with pd.ExcelWriter(excel_path) as writer:
        create_df_safely(analysis_results['bound_exo'], 'bound_exo_track', writer)
        create_df_safely(analysis_results['bound_pol'], 'bound_pol_track', writer)
        create_df_safely(analysis_results['dark_exo'], 'dark_exo_track', writer)
        create_df_safely(analysis_results['dark_pol'], 'dark_pol_track', writer)

    print(f"Segmented track data saved to: {excel_path}")


def extract_event_segments(data, start_times, end_times, dark_segment_time):
    """Extracts data segments for each event, padded with dark time."""
    segments = []
    start_times_segment = np.array(start_times) - dark_segment_time
    end_times_segment = np.array(end_times) + dark_segment_time
    
    time = np.round(data['bp_time'], 3)
    basepairs = np.round(data['bp'], 3)
    basepairs_fitted = np.round(data['cp_basepair_interp'], 3)
    
    for i, start_time in enumerate(start_times_segment):
        end_time = end_times_segment[i]
        ind_force = np.where((time <= end_time) & (time >= start_time))
        ind_img = np.where((data['time_intens'] <= end_time) & (data['time_intens'] >= start_time))
        
        force_data = np.vstack((
            time[ind_force],
            basepairs[ind_force],
            basepairs_fitted[ind_force]
        )).T
        
        image_data = np.vstack((
            data['time_intens'][ind_img],
            data['intensity'][ind_img],
            data['intensity_step'][ind_img],
            data['intensity_raw'][ind_img]
        )).T
        
        if force_data.shape[0] > 1 and image_data.shape[0] > 1:
            segments.append(np.array([force_data, image_data], dtype=object))
            
    return segments


def segment_plots(segment, ax=None):
    if ax is None: ax = plt.gca()
        
    time_ROI, basepairs_ROI, basepairs_fitted_ROI = segment[0].T
    time_intens_ROI, intensity_ROI, intensity_step_ROI, intensity_raw_ROI = segment[1].T
    
    twin1 = ax.twinx()
    twin2 = ax.twinx()
    
    ax.plot(time_ROI, basepairs_fitted_ROI, "black", label="step fitted data", linewidth=1)
    ax.plot(time_ROI, basepairs_ROI, "lightgrey", label="raw data", linewidth=0.5)
    ax.set_xlabel('Time/s')
    ax.set_ylabel('Basepairs')
    ax.set_xlim(segment[0][0, 0], segment[0][-1, 0])
    
    twin1.plot(time_intens_ROI, intensity_raw_ROI, color='lightgreen', linewidth=0.5)
    twin2.plot(time_intens_ROI, intensity_step_ROI, color='green', linewidth=1)
    twin1.yaxis.set_ticklabels([])
    twin1.yaxis.set_ticks([])
    twin2.set_ylabel('Fluorescence Intensity (a.u.)')
    
    twin3 = ax.twinx()
    twin3.set_ylim(0.02, 1.05)
    twin3.fill_between(time_intens_ROI, intensity_ROI, 0, step="pre", alpha=0.2, color='gray')
    twin3.axis('off')

    return ax


def segment_rate_plots(segment, ax=None):
    if ax is None: ax = plt.gca()

    time_ROI, basepairs_ROI, basepairs_fitted_ROI = segment[0].T
    time_intens_ROI, intensity_ROI, intensity_step_ROI, intensity_raw_ROI = segment[1].T

    ax1 = ax.twinx()
    ax2 = ax.twinx()

    ax.plot(time_ROI[:-1], np.diff(basepairs_fitted_ROI) / np.diff(time_ROI), "red", label="step fitted data", linewidth=1)
    ax.set_xlabel('Time/s')
    ax.set_ylabel('Activity Burst(bp/s)')
    ax.set_xlim(segment[0][0, 0], segment[0][-1, 0])

    ax1.plot(time_intens_ROI, intensity_raw_ROI, color='lightgreen', linewidth=0.5)
    ax1.plot(time_intens_ROI, intensity_step_ROI, color='green', linewidth=0.5)
    ax1.set_ylabel('Fluorescence Intensity (a.u.)')
    
    ax2.set_ylim(0.02, 1.05)
    ax2.fill_between(time_intens_ROI, intensity_ROI, 0, step="pre", alpha=0.1, color='gray')
    ax2.axis('off')
    
    return ax


def plot_all_segments(segments, segment_type, output_dir, name_prefix, show_plots, plot_function, suptitle, figsize, ncols):
    """Generic function to create and save subplot grids for segments."""
    if not segments:
        print(f"No {segment_type} segments to plot.")
        return

    nrows = (len(segments) + ncols - 1) // ncols
    fig, axs = plt.subplots(nrows, ncols, figsize=figsize, facecolor='w', edgecolor='k', squeeze=False)
    fig.subplots_adjust(hspace=0.8, wspace=0.3)
    fig.suptitle(suptitle)

    for i, seg in enumerate(segments):
        ax = axs.flat[i]
        plot_function(seg, ax)
    
    for i in range(len(segments), nrows * ncols):
        axs.flat[i].set_visible(False)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    png_path = os.path.join(output_dir, f"{name_prefix}-{segment_type}-segments.png")
    fig.savefig(png_path, dpi=300)
    print(f"{segment_type.capitalize()} segments figure saved to: {png_path}")

    if show_plots:
        plt.show()
    plt.close(fig)


def analyze_correlation_heatmap(data, split_idx, output_dir, name_prefix, show_plots):
    """Generates, analyzes, and plots the correlation heatmap."""
    df_activity_data = pd.DataFrame({'Time': data['bp_time'][:-1], 'Activity': np.diff(data['cp_basepair_interp'])/np.diff(data['bp_time'])})
    bins = [-np.inf, -15, 20, np.inf]
    labels = [1, 0, 1]
    df_activity_data['BinarizedActivity'] = pd.cut(df_activity_data['Activity'], bins=bins, labels=labels, ordered=False).astype(int)

    df_intensity_data = pd.DataFrame({'Time': data['time_intens'], 'Intensity': data['intensity']})

    f = interp1d(df_activity_data['Time'], df_activity_data['BinarizedActivity'], kind='nearest', bounds_error=False, fill_value="extrapolate")
    df_activity_data_interpolated = pd.DataFrame({'Time': df_intensity_data['Time'], 
                                                  'InterpolatedBinarizedActivity': f(df_intensity_data['Time']).astype(int)})

    df_merged = pd.merge(df_intensity_data, df_activity_data_interpolated, on='Time')

    conditions = [
        (df_merged['Intensity'] == 1) & (df_merged['InterpolatedBinarizedActivity'] == 1),
        (df_merged['Intensity'] == 0) & (df_merged['InterpolatedBinarizedActivity'] == 0),
        (df_merged['Intensity'] == 0) & (df_merged['InterpolatedBinarizedActivity'] == 1),
        (df_merged['Intensity'] == 1) & (df_merged['InterpolatedBinarizedActivity'] == 0),
    ]
    values = [1, -1, 0.5, -0.5]
    df_merged['CustomCorrelation'] = np.select(conditions, values)

    # --- Matrix Generation ---
    matrix_size = len(df_merged)
    matrix = np.zeros((matrix_size, matrix_size))
    
    # This logic from the original notebook fills square blocks on the anti-diagonal.
    start_idx = 0
    for i in range(1, matrix_size):
        if df_merged['CustomCorrelation'].iat[i] != df_merged['CustomCorrelation'].iat[i-1]:
            end_idx = i
            val = df_merged['CustomCorrelation'].iat[i-1]
            matrix[matrix_size - end_idx:matrix_size - start_idx, start_idx:end_idx] = val
            start_idx = i
    # Last block
    end_idx = matrix_size
    val = df_merged['CustomCorrelation'].iat[matrix_size-1]
    matrix[matrix_size - end_idx:matrix_size - start_idx, start_idx:end_idx] = val

    # Flip vertically to move anti-diagonal to main diagonal for easier plotting/splitting
    matrix_flipped = np.flip(matrix, 0)

    # --- Plot Full Heatmap ---
    fig_full = plt.figure(figsize=(7, 6))
    ax_full = fig_full.add_subplot(111)
    cmap = LinearSegmentedColormap.from_list('custom', [(0, 'cyan'), (0.25, 'blue'), (0.5, 'white'), (0.75, 'pink'), (1, 'red')], N=5)
    norm = Normalize(vmin=-1, vmax=1)
    
    X, Y = np.meshgrid(np.arange(matrix_flipped.shape[1] + 1), np.arange(matrix_flipped.shape[0] + 1))
    im = ax_full.pcolormesh(X, Y, matrix_flipped, cmap=cmap, norm=norm, shading='auto', rasterized=True)
    ax_full.set_xlabel("Fluorescence Measurement Time (pixel)")
    ax_full.set_ylabel("DNAp Activity Time (pixel)")
    ax_full.set_title('Correlation Heatmap')
    ax_full.set_aspect('equal')
    fig_full.colorbar(im, ax=ax_full, shrink=0.8)

    full_heatmap_base = os.path.join(output_dir, f"{name_prefix}-correlation_heatmap_full")
    plt.savefig(f"{full_heatmap_base}.png", format='png', dpi=300, bbox_inches='tight')
    plt.savefig(f"{full_heatmap_base}.eps", format='eps', dpi=300, bbox_inches='tight')
    print(f"Full heatmap saved to: {full_heatmap_base}.(png/eps)")
    if show_plots: plt.show()
    plt.close(fig_full)


    # --- Cluster Counting ---
    diag_values = df_merged['CustomCorrelation'].to_numpy()
    def count_clusters(diagonal_values, min_cluster_length=3):
        cluster_counts = {-1: 0, -0.5: 0, 0.5: 0, 1: 0}
        if len(diagonal_values) == 0: return cluster_counts

        current_val = diagonal_values[0]
        current_len = 0
        for val in diagonal_values:
            if val == current_val:
                current_len += 1
            else:
                if current_val != 0 and current_len >= min_cluster_length:
                    if current_val in cluster_counts:
                        cluster_counts[current_val] += 1
                current_val = val
                current_len = 1
        if current_val != 0 and current_len >= min_cluster_length:
            if current_val in cluster_counts:
                cluster_counts[current_val] += 1
        return cluster_counts
    
    cluster_counts_1 = count_clusters(diag_values[:split_idx])
    cluster_counts_2 = count_clusters(diag_values[split_idx:])
    
    heatmap_data_path = os.path.join(output_dir, f"{name_prefix}-cluster_counts.json")
    with open(heatmap_data_path, 'w') as f:
        json.dump({"exo_cluster_counts": cluster_counts_1, "pol_cluster_counts": cluster_counts_2}, f)
    print(f"Cluster counts data saved to: {heatmap_data_path}")

    # --- Plot Split Heatmaps ---
    fig_split, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    matrix_1 = matrix_flipped[:split_idx, :split_idx]
    matrix_2 = matrix_flipped[split_idx:, split_idx:]
    
    X1, Y1 = np.meshgrid(np.arange(matrix_1.shape[1] + 1), np.arange(matrix_1.shape[0] + 1))
    ax1.pcolormesh(X1, Y1, matrix_1, cmap=cmap, norm=norm, shading='auto', rasterized=True)
    ax1.set_title("Exo Events")
    ax1.set_xlabel("Fluorescence Measurement Time (pixel)")
    ax1.set_ylabel("DNAp Activity Time (pixel)")
    ax1.set_aspect('equal')

    x_coords2 = np.arange(split_idx, matrix_flipped.shape[1] + 1)
    y_coords2 = np.arange(split_idx, matrix_flipped.shape[0] + 1)
    X2, Y2 = np.meshgrid(x_coords2, y_coords2)

    if X2.shape[0] > matrix_2.shape[0] and X2.shape[1] > matrix_2.shape[1]:
        ax2.pcolormesh(X2[:matrix_2.shape[0]+1, :matrix_2.shape[1]+1], Y2[:matrix_2.shape[0]+1, :matrix_2.shape[1]+1], matrix_2, cmap=cmap, norm=norm, shading='auto', rasterized=True)
    else:
        ax2.pcolormesh(matrix_2, cmap=cmap, norm=norm, shading='auto', rasterized=True) # Fallback
    
    ax2.set_title("Pol Events")
    ax2.set_xlabel("Fluorescence Measurement Time (pixel)")
    ax2.set_aspect('equal')
    
    fig_split.tight_layout()
    heatmap_filename_base = os.path.join(output_dir, f"{name_prefix}-correlation_exo_pol_heatmap")
    plt.savefig(f"{heatmap_filename_base}.eps", format='eps', dpi=300, bbox_inches='tight')
    plt.savefig(f"{heatmap_filename_base}.png", format='png', dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to: {heatmap_filename_base}.(eps/png)")

    if show_plots: plt.show()
    plt.close(fig_split)


def main():
    # --- Path setup ---
    # To make the script runnable from anywhere, we define default paths
    # relative to the script's location, assuming a structure like:
    # project_root/
    #  |- code/ (this script)
    #  |- data/
    #  |- results/
    script_dir = os.path.dirname(os.path.realpath(__file__))
    project_root = os.path.dirname(script_dir) # This will be the 'code' directory's parent

    default_cor_file = os.path.join(project_root, 'results/ChangePoints_Results/OTdata_example_30nM DNAp + trx  + 625uM dNTPs-cycle#01-processedData-correlated_analysis.xlsx')
    default_cor_file2 = os.path.join(project_root, 'results/image_data_example 30nM DNAp + trx  + 625uM dNTPs-cycle#1-Intensity_along_DNAp-filtered.txt')
    default_cor_file3 = os.path.join(project_root, 'results/stepfinding_results/image_data_example 30nM DNAp + trx  + 625uM dNTPs-cycle#1-Intensity_along_DNAp-filtered_fit_trace.csv')
    default_output_dir = os.path.join(project_root, 'results/final_analysis')

    parser = argparse.ArgumentParser(description="Correlated segment analysis of DNA polymerase activity.")
    parser.add_argument('--cor_file', type=str, default=default_cor_file,
                        help='Path to the correlated data Excel file from step 2.')
    parser.add_argument('--cor_file2', type=str, default=default_cor_file2,
                        help='Path to the raw intensity data file.')
    parser.add_argument('--cor_file3', type=str, default=default_cor_file3,
                        help='Path to the filtered intensity trace file.')
    parser.add_argument('-o', '--output_dir', type=str, default=default_output_dir, help='Output directory for results and figures.')
    parser.add_argument('--force_transition_time', type=float, default=62.5, help='Time of force transition in seconds.')
    parser.add_argument('--dark_segment_time', type=float, default=1.5, help='Additional dark time to include before and after a bound event.')
    parser.add_argument('--heatmap_split_pixel', type=int, default=260, help='Pixel index to split the correlation heatmap for exo/pol analysis.')
    parser.add_argument('--no-show-plots', action='store_true', help="Do not display plots, only save them.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    base_filename = os.path.basename(args.cor_file)
    name_prefix = base_filename.split('#')[0] if '#' in base_filename else os.path.splitext(base_filename)[0]

    data = load_data(args.cor_file, args.cor_file2, args.cor_file3)

    plot_correlation_overview(data, args.output_dir, name_prefix, show_plots=not args.no_show_plots)
    plot_activity_burst(data, args.output_dir, name_prefix, show_plots=not args.no_show_plots)

    start_time, end_time = find_binding_events(data['intensity'], data['time_intens'])
    
    start_time_exo, end_time_exo, start_time_pol, end_time_pol = distinguish_exo_pol_events(
        start_time, end_time, data['time_intens'], data['intensity'], args.force_transition_time
    )

    dark_exo_start, dark_exo_end = find_dark_events(start_time_exo, end_time_exo, data['time_intens'][0], args.force_transition_time, data['time_intens'][0])
    dark_pol_start, dark_pol_end = find_dark_events(start_time_pol, end_time_pol, args.force_transition_time, data['time_intens'][-1], args.force_transition_time)

    rate_bound_exo = calculate_segment_rates(data['bp_time'], data['bp'], start_time_exo, end_time_exo)
    rate_bound_pol = calculate_segment_rates(data['bp_time'], data['bp'], start_time_pol, end_time_pol)
    rate_unbound_exo = calculate_segment_rates(data['bp_time'], data['bp'], dark_exo_start, dark_exo_end)
    rate_unbound_pol = calculate_segment_rates(data['bp_time'], data['bp'], dark_pol_start, dark_pol_end)

    analysis_results = {
        'bound_exo': {'track_Nmr': np.arange(1, len(start_time_exo) + 1), 'start_time/s': start_time_exo, 'end_time/s': end_time_exo, 'track_duration/s': (np.array(end_time_exo) - np.array(start_time_exo)), 'track_rate': rate_bound_exo},
        'bound_pol': {'track_Nmr': np.arange(1, len(start_time_pol) + 1), 'start_time/s': start_time_pol, 'end_time/s': end_time_pol, 'track_duration/s': (np.array(end_time_pol) - np.array(start_time_pol)), 'track_rate': rate_bound_pol},
        'dark_exo': {'track_Nmr': np.arange(1, len(dark_exo_start) + 1), 'dark_exo_start/s': dark_exo_start, 'dark_exo_end/s': dark_exo_end, 'track_duration/s': (np.array(dark_exo_end) - np.array(dark_exo_start)), 'track_rate': rate_unbound_exo},
        'dark_pol': {'track_Nmr': np.arange(1, len(dark_pol_start) + 1), 'dark_pol_start/s': dark_pol_start, 'dark_pol_end/s': dark_pol_end, 'dark_track_duration/s': (np.array(dark_pol_end) - np.array(dark_pol_start)), 'track_rate': rate_unbound_pol}
    }
    save_segmented_data(args.output_dir, name_prefix, analysis_results)
    
    segment_exo = extract_event_segments(data, start_time_exo, end_time_exo, args.dark_segment_time)
    segment_pol = extract_event_segments(data, start_time_pol, end_time_pol, args.dark_segment_time)
    
    exo_npy_path = os.path.join(args.output_dir, f"{name_prefix}-exo-segments.npy")
    np.save(exo_npy_path, segment_exo, allow_pickle=True)
    print(f"Exo segments saved to: {exo_npy_path}")

    pol_npy_path = os.path.join(args.output_dir, f"{name_prefix}-pol-segments.npy")
    np.save(pol_npy_path, segment_pol, allow_pickle=True)
    print(f"Pol segments saved to: {pol_npy_path}")

    # Plotting segments
    plot_all_segments(segment_exo, 'exo', args.output_dir, name_prefix, not args.no_show_plots, segment_plots, 'Exo Segments', (12, 8), 3)
    plot_all_segments(segment_pol, 'pol', args.output_dir, name_prefix, not args.no_show_plots, segment_plots, 'Pol Segments', (12, 8), 4)
    plot_all_segments(segment_exo, 'exo-rate', args.output_dir, name_prefix, not args.no_show_plots, segment_rate_plots, 'Exo Segments Activity Burst Rate', (12, 8), 3)
    plot_all_segments(segment_pol, 'pol-rate', args.output_dir, name_prefix, not args.no_show_plots, segment_rate_plots, 'Pol Segments Activity Burst Rate', (12, 8), 4)

    analyze_correlation_heatmap(data, args.heatmap_split_pixel, args.output_dir, name_prefix, show_plots=not args.no_show_plots)


if __name__ == '__main__':
    main()

