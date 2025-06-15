#!/usr/bin/env python
# coding: utf-8
"""
Correlation Image Force Analysis Script
Converted from Jupyter notebook: 2_Correlation_image_force.ipynb

This script performs correlation analysis between DNA polymerase trace from images
and force measurements, including change-point detection and visualization.
"""
# import all the libraries 
# python==3.8; jupyterlab==3.0.12; lumicks.pylake==0.8.1; matplotlib==3.3.4; more-itertools==8.7.0;
# npTDMS==1.1.0; numpy==1.20.1; opencv-python==4.5.1.48; pandas==1.2.3; scipy==1.6.1; tifffile==2021.3.5
import os
import sys
import numpy as np
import cv2
import pylab as pl
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import collections  as mc
from matplotlib import colors as mcolors
import pandas as pd
from nptdms import TdmsFile
from scipy.signal import savgol_filter
import tifffile as tif
from scipy import interpolate
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

# Matplotlib configuration for script execution
plt.ioff()  # Turn off interactive mode

def main():
    """Main function to run the correlation analysis"""
    
    # Determine project root to make file paths relative to it
    # Assumes the script is in a subfolder of the project root (e.g., 'python_scripts')
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    print("Starting Correlation Image Force Analysis...")

    # read raw image data of .tdms file with TdmsFile function
    # kymo_filename = input('please type in the file name:' )
    kymo_filename_rel = 'example_dataset/image data/image_data_example 30nM DNAp + trx  + 625uM dNTPs.tdms'
    kymo_filename = os.path.join(project_root, kymo_filename_rel)

    # kymo_cycle = str('#') + input('please type in the processing cycle(e.g: 1):' )
    kymo_cycle = str('#1') 

    # access .tdms file and convert .tdms to image data
    metadata = TdmsFile.read_metadata(kymo_filename)
    # print(metadata.properties)
    width = metadata.properties['Pixels per line']
    px_size = float(metadata.properties['Scan Command.Scan Command.scanning_axes.0.pix_size_nm'])
    px_dwell_time = float(metadata.properties['Scan Command.PI Fast Scan Command.pixel_dwell_time_ms'])
    inter_frame_time = float(metadata.properties['Scan Command.PI Fast Scan Command.inter_frame_wait_time_ms'])

    tdms_file = TdmsFile(kymo_filename)
    kymo_time = tdms_file['Data']['Time (ms)'][:]
    kymo_time = np.array([int(i) for i in kymo_time])
    kymo_position = tdms_file['Data']['Actual position X (um)'][:]
    kymo_position = np.array([int(i) for i in kymo_position])
    height = len(kymo_time)/width
    time_per_line = kymo_time[-1]/height # this is calculated by (time_per_line = acquisition time/all lines)

    chn_r = tdms_file['Data']['Pixel ch 1'][:]
    chn_r = np.array([int(i) for i in chn_r])
    chn_g = tdms_file['Data']['Pixel ch 2'][:]
    chn_g = np.array([int(i) for i in chn_g])
    chn_b = tdms_file['Data']['Pixel ch 3'][:]
    chn_b = np.array([int(i) for i in chn_b])

    chn_rgb = np.vstack((chn_r, chn_g, chn_b)).T
    img = chn_rgb.reshape((int(height), int(width), 3))
    img = img.transpose((1, 0, 2))
    img = img.astype(np.uint16)

    # show tif image
    plt.figure(kymo_filename)
    plt.imshow(img.astype('uint16'),vmax = 10)
    plt.xlabel('Time/px')
    plt.ylabel('Position/px')
    plt.tight_layout()
    # plt.show()
    plt.close()

    # save .tdms file to tif image to a specific folder
    output_dir = os.path.join(os.path.dirname(kymo_filename), 'results')
    os.makedirs(output_dir, exist_ok=True)
    file_stem = os.path.splitext(os.path.basename(kymo_filename))[0]

    # Save the tif image to a specific folder
    plot_filename = f"{file_stem}.tiff"
    plot_path = os.path.join(output_dir, plot_filename)

    tif.imwrite(plot_path, img)

    # Plot DNA polymerase from image data
    b,g,r = cv2.split(img)

    # the kymograph can be plotted in a normal image or in a reverse black-white image
    vmax =12
    fig, ax=plt.subplots(figsize=(10,4))
    ax.imshow(g.astype('uint16'),cmap='gray',vmax = vmax,aspect ="auto")
    # plt.imshow(100-g.astype('uint16'),cmap='gray',vmin = 95,aspect ="auto")

    plt.xlabel("Time/px")
    plt.ylabel("Position/px")
    # plt.title("DNA Polymerase on ssDNA")

    # this number defines the ROI (region of interest) of image in px; should be tuned to properly display the image
    kymo_xlim_left = 5
    kymo_xlim_right = 550
    kymo_ylim_top = 25
    kymo_ylim_bottom = 83

    ax.set_xlim(kymo_xlim_left,kymo_xlim_right)
    ax.set_ylim(kymo_ylim_bottom,kymo_ylim_top)
    
    # read and plot the DNA trace calculated from force measurement (1_CalculatingDNApTrace_OT.py)
    # trace_file = input('please type in the file name:' )
    trace_file_rel = 'example_dataset/force data/results/OTdata_example_30nM DNAp + trx  + 625uM dNTPs-cycle#01-processedData.xlsx'
    trace_file = os.path.join(project_root, trace_file_rel)
    trace = pd.read_excel(trace_file, engine='openpyxl')
    print(trace.head())

    # the following step is intended to plot the DNA junction on top of DNAp trajectory
    # The time of image is used for reference and displayed
    # Tune x_offset and y_offset to get the proper alignment
    plt.xlabel("Time/px")
    plt.ylabel("Position/px")
    plt.title("DNA Polymerase Trace")
    x_offset_searching = -28 # this is done by align starting time in pixel （larger absolute value moves to left）
    y_offset_searching = 25 # this is done by align starting position in pixel 
    x_cali = 1000/time_per_line  # this is calculated by (1s/time per line， e.g.：163.8ms)
    y_cali = 1000/px_size # this is calculated by (1um/pixel size, in this case pixlsize = 75nm)
    trace_time = trace['time']/1000 * x_cali + x_offset_searching
    trace_time = trace_time.dropna()
    position = pd.to_numeric(trace['junction_position_all'], errors='coerce') 
    position = position * y_cali + y_offset_searching
    position = position.dropna()

    ax.plot(trace_time,position,'yellow',linewidth = 2, label = 'first trial correlation')
    # plt.ylim(25,85)
    # plt.legend()
    # plt.savefig(kymo_filename.replace('.tdms', '-cycle')+ kymo_cycle+'-overlap.eps', format='eps', dpi=300,bbox_inches='tight')

    # Save the plot
    plot_filename = f"{file_stem}-cycle{kymo_cycle}-overlap.png"
    plot_output_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_output_path, dpi=300,bbox_inches='tight')
    # plt.show()
    plt.close(fig)

    # In case after running the above code, you see the DNAp moving in an opposite direction, go back to step1 (1_CalculatingDNApTrace_OT.py) and recalculate the junction position

    # OPTIONAL STEP
    # this step is intended to automatically find the optimal x_offset and y_offset by searching in a confined space, 
    # wiht a searching goal of maximizing the intensity along DNAp trace
    x_offset_ls = []
    y_offset_ls = []
    intens_sum_ls = []
    for i_x in np.arange(-5,6):
        for i_y in np.arange (-2,3):
            x_offset = x_offset_searching + i_x
            y_offset = y_offset_searching + i_y
            # print(x_offset,y_offset)
            trace_time = trace['time']/1000 * x_cali + x_offset
            trace_time = trace_time.dropna()
            position = pd.to_numeric(trace['junction_position_all'], errors='coerce') 
            position = position * y_cali + y_offset
            position = position.dropna()
            

            # In order to calculate the the intensity of protein trace as a function of time (px), we first need the locate the protein trace based on DNA trajectory
            # this is done by reading the coordinates (time, position) of protein trace in kymograph based on the overlapping of force data and image data.
            # Since the sampling rate is different between force data and image data, therefore interpolate is performed
            func = interpolate.interp1d(trace_time,position,kind='slinear',fill_value="extrapolate")

            # x defines the pixel in time
            # x = np.arange(round(trace_time.iloc[0]),round(trace_time.iloc[-1]))
            # to exclude part of the trajectory that is too close to the bright bead
            x = np.arange(round(trace_time.iloc[100]),round(trace_time.iloc[-1]))
            # x1 here traces back to time in unit of s
            x1 = (x-x_offset)/x_cali
            y = func(x)
            y = np.rint(y).astype(int)
            # the sum of the intensity along DNAp
            intens_sum = np.sum(g[y,x])
            
            x_offset_ls.append(x_offset)
            y_offset_ls.append(y_offset)
            intens_sum_ls.append(intens_sum)
            
    x_offset_ls = np.array(x_offset_ls)       
    y_offset_ls = np.array(y_offset_ls)
    intens_sum_ls = np.array(intens_sum_ls)

    # this step is to find out the index of the maximum value of intensity
    intens_ind = np.where(intens_sum_ls[:] == np.amax(intens_sum_ls[:]))
    x_offset_optimal = x_offset_ls[intens_ind]
    y_offset_optimal = y_offset_ls[intens_ind]
    print(x_offset_optimal)
    print(y_offset_optimal)
    trace_time = trace['time']/1000 * x_cali + x_offset_optimal
    trace_time = trace_time.dropna()
    position = pd.to_numeric(trace['junction_position_all'], errors='coerce') 
    position = position * y_cali + y_offset_optimal
    position = position.dropna()
    ax.plot(trace_time,position,'b',linewidth = 2, label = 'optimal correlation')
    ax.legend()   

    # Save the plot
    plot_filename = f"{file_stem}-cycle{kymo_cycle}-overlap-optimizing.png"
    plot_output_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_output_path, dpi=300,bbox_inches='tight')
    # plt.show()
    plt.close()

    # In case the optimazation fails, you can manually set the optimal offset, otherwise don't run this code
    x_offset_optimal = x_offset_searching
    y_offset_optimal = y_offset_searching 

    # OPTIONAL STEP
    # the fllowing code is intended to plot the intensity as a function of x-offset and y-offset in 3d
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #surf = ax.plot_trisurf(Xs, Ys, Zs, cmap=cm.jet, linewidth=0)
    surf = ax.plot_trisurf(x_offset_ls, y_offset_ls, intens_sum_ls, cmap=cm.jet, linewidth=0)
    fig.colorbar(surf)
    ax.set_xlabel('X-offset')
    ax.set_ylabel('Y-offset')
    ax.set_zlabel('Intensity of DNAp Trace')

    # plt.show()
    plt.close(fig)

    # In order to calculate the the intensity of protein trace as a function of time (px), we first need the locate the protein trace based on DNA trajectory
    # this is done by reading the coordinates (time, position) of protein trace in kymograph based on the overlapping of force data and image data.
    # Since the sampling rate is different between force data and image data, therefore interpolate is performed
    func = interpolate.interp1d(trace_time,position,kind='slinear',fill_value="extrapolate")
    # x defines the pixel in time
    x = np.arange(round(trace_time.iloc[0]),round(trace_time.iloc[-1]))
    # x1 here traces back to time in unit of s
    x1 = (x-x_offset_optimal)/x_cali
    y = func(x)
    y = np.rint(y).astype(int)
    plt.figure(figsize=(10,3))
    plt.plot(x,y)
    plt.xlabel("Time/px")
    plt.ylabel("Position/px")
    plt.tight_layout()
    # plt.show()
    plt.close()

    # this step is to plot the intensity on the axis of DNAp trace
    # img has been split into 3 channels, b,g,r; so here green channel is used
    ins = g[y,x]
    plt.figure(figsize=(8,3))
    plt.plot(x1,ins)
    plt.xlabel("Time/s")
    plt.ylabel("Intensity/photon")
    plt.title("Intensity along DNA polymerase Trajectory")
    plt.tight_layout()
    # plt.ylim(0,30)
    # plt.show()
    plt.close()

    # sum up 5 pixels around DNAp trace
    list = [y-2, y-1, y, y+1, y+2]
    intensity = np.zeros((len(list),len(x)), dtype=int)
    all_intensity = np.zeros((len(x),), dtype=int)
    for i, value in enumerate(list):
        intensity[i] = g[value,x]
    #     print(intensity)
        all_intensity += g[value,x]
        
    # print(intensity)
    # print(all_intensity)
    plt.figure(figsize=(8,3))
    plt.plot(x1,all_intensity)
    plt.xlabel("Time/s")
    plt.ylabel("Intensity/photon")
    plt.title("Intensity along DNA polymerase Trajectory")
    plt.tight_layout()
    plt.xlim((kymo_xlim_left-x_offset_optimal)/x_cali,(kymo_xlim_right-x_offset_optimal)/x_cali)
    # plt.ylim(0,70)
    # plt.show()
    plt.close()

    # this code block is intended to overlap the basepair-time and raw intensity -time plots
    basepair = pd.read_excel(trace_file)
    bp_time = basepair['time']/1000
    bp_time = bp_time.dropna()
    bp = basepair['basepairs']
    bp = bp.dropna()

    fig, ax1 = plt.subplots(figsize=(8,3))

    ax2 = ax1.twinx()
    ax1.plot(bp_time,bp,color='black',linewidth=1, label='basepairs')
    ax1.set_xlabel('Time/s')
    ax1.set_ylabel('Basepairs', color='black')

    ax2.plot(x1,all_intensity,color='green',linewidth=0.2, label='fluorescence intensity')

    ax2.set_ylabel('Fluorescence intensity(a.u.)', color='green')
    plt.xlim((kymo_xlim_left-x_offset_optimal)/x_cali,(kymo_xlim_right-x_offset_optimal)/x_cali)
    # plt.legend()
    # plt.xlim(35,70)
    plt.tight_layout()

    # Save the plot
    plot_filename = f"{file_stem}-cycle{kymo_cycle}-rawBp-intensity_along_DNAp_Trajectory.png"
    plot_output_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_output_path, dpi=300,bbox_inches='tight')
    # plt.show()
    plt.close(fig)

    # plt.close()

    # ## To extract the DNAp-bound events

    # Apply a Savitzky-Golay filter to an array.
    signal = all_intensity
    signal_filter_window = 7
    signal_filter = savgol_filter(signal[:], signal_filter_window, 3) # filter window size of 47 gives us a sampling rate of 0.1HZ
    # dt = x[1]-x[0]
    # signal_filter_grad = np.gradient(signal_filter, dt)
    # plt.figure()

    plt.figure(figsize=(8,3))
    plt.plot(x1,all_intensity)
    plt.plot(x1,signal_filter)
    plt.xlabel("Time/s")
    plt.ylabel("Intensity/photon")
    plt.title("Intensity along DNA polymerase Trajectory")
    plt.tight_layout()
    plt.xlim((kymo_xlim_left-x_offset_optimal)/x_cali,(kymo_xlim_right-x_offset_optimal)/x_cali)
    # plt.ylim(0,30)
    # plt.show()
    plt.close()

    # plt.savefig(kymo_filename.replace('.tdms', '-cycle') + kymo_cycle + 'Intensity along DNA polymerase Trajectory'+'.eps', format='eps', dpi=300, bbox_inches='tight')

    # Save the plot
    plot_filename = f"{file_stem}-cycle{kymo_cycle}-Intensity along DNA polymerase Trajectory.png"
    plot_output_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_output_path, dpi=300,bbox_inches='tight')
    # plt.show()
    plt.close()

    # the following code blocks are to find out the intensity of the background, by selecting a region in the kymopgraph and calculating the intensity
    # step 1: selecting the ROI of background
    plt.figure(figsize=(10,3))
    # the kymograph can be plotted in a normal way or in a reverse way
    # plt.imshow(g.astype('uint16'),cmap='gray',vmax = 10,aspect ="auto")
    plt.imshow(g.astype('uint16'),cmap='gray',vmax = vmax,aspect ="auto")
    currentAxis=plt.gca()
    # for consistentence, background window is selected with a fixed size (5px * 200px)
    rec_x,rec_y,rec_w,rec_h = [330,50,len(list),200]
    rect=patches.Rectangle((rec_x,rec_y),rec_h,rec_w,linewidth=1,edgecolor='r',facecolor='none',fill=False)
    currentAxis.add_patch(rect)

    plt.xlabel("Time/px")
    plt.ylabel("Position/px")
    plt.title("DNA Polymerase Trace")
    plt.ylim(82,22)
    plt.xlim(kymo_xlim_left,kymo_xlim_right)
    # plt.xlim(0,560)
    plt.tight_layout()

    # plt.savefig(kymo_filename.replace('.tdms', '-cycle') + kymo_cycle + '-background Intensity'+'.eps', format='eps', dpi=300, bbox_inches='tight')
    # plt.close()

    # Save the plot
    plot_filename = f"{file_stem}-cycle{kymo_cycle}-background Intensity.png"
    plot_output_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_output_path, dpi=300,bbox_inches='tight')
    # plt.show()
    plt.close()

    # step 2: calculating the intensity of the background
    list = [rec_y,rec_y+1,rec_y+2,rec_y+3,rec_y+4]
    rec_x_list = np.arange(rec_x, rec_x+rec_h)
    bagrnd_intensity = np.zeros((rec_w,len(rec_x_list)), dtype=int)
    bagrnd_all_intensity = np.zeros((len(rec_x_list),), dtype=int)
    for i, value in enumerate(list):
        bagrnd_intensity[i] = g[value,rec_x_list]
    #     print(intensity)
        bagrnd_all_intensity += g[value,rec_x_list]


    fig, ax=plt.subplots(figsize=(8,3))
    ax.plot(rec_x_list,bagrnd_all_intensity,linestyle='solid',linewidth=0.3, markersize=0.5,label='Background Signal')

    ax.set_xlim(rec_x,rec_x+rec_h)
    ax.set_xlabel("Time/px")
    ax.set_ylabel("Intensity/photon")
    ax.set_title("Background Intensity")
    plt.tight_layout()
    # plt.ylim(0,30)
    # plt.show()

    # Apply a Savitzky-Golay filter to an array.
    bagrnd_filter_window = 25
    threshold_sigma  = 3

    bagrnd_signal = bagrnd_all_intensity 
    bagrnd_signal_filter = savgol_filter(bagrnd_signal[:], bagrnd_filter_window, 3) # filter window size of 25 gives us a sampling rate of 0.2HZ
    # dt = x[1]-x[0]
    # signal_filter_grad = np.gradient(signal_filter, dt)
    # plt.figure()
    ax.plot(rec_x_list,bagrnd_signal_filter,linestyle='solid',linewidth=1, markersize=1,label='Background Signal Filter')
    ax.axhline(np.average(bagrnd_signal_filter)+threshold_sigma*np.std(bagrnd_signal_filter), color='k', linestyle='--',label= str(threshold_sigma) + ' Sigma&filterSize ' + str(bagrnd_filter_window))
    ax.axhline(np.average(bagrnd_signal_filter)-threshold_sigma*np.std(bagrnd_signal_filter), color='k', linestyle='--')

    # plt.axhline(np.average(bagrnd_signal_filter)+np.std(bagrnd_signal_filter), color='r', linestyle='--',label='1 sigma')
    # plt.axhline(np.average(bagrnd_signal_filter)-np.std(bagrnd_signal_filter), color='r', linestyle='--')

    # plt.ylim(0,80)
    ax.set_xlim(rec_x,rec_x+rec_h)
    ax.legend(loc='upper right')
    fig
    # plt.tight_layout()

    # plt.savefig(kymo_filename.replace('.tdms', '-cycle') + kymo_cycle + 'background Intensity along DNA polymerase Trajectory'+'.eps', format='eps', dpi=300, bbox_inches='tight')

    # Save the plot
    plot_filename = f"{file_stem}-cycle{kymo_cycle}-background Intensity along DNA polymerase Trajectory.png"
    plot_output_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_output_path, dpi=300,bbox_inches='tight')
    # plt.show()
    plt.close()

    plt.figure(figsize=(10,1))
    # print(image)
    # plt.imshow(intensity.astype('uint8'),aspect ="auto")
    plt.imshow(bagrnd_intensity.astype('uint8'),cmap='gray',vmax = vmax,interpolation='nearest', aspect='auto')

    plt.axis('off')
    # plt.xlim(0, )
    plt.tight_layout()
    # plt.ylim(0,7)

    # Save the plot
    plot_filename = f"{file_stem}-cycle{kymo_cycle}-Background Intensity.png"
    plot_output_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_output_path, dpi=300,bbox_inches='tight')
    # plt.show()
    plt.close()

    plt.figure(figsize=(10,1))
    # print(image)
    # plt.imshow(intensity.astype('uint8'),aspect ="auto")
    plt.imshow(intensity.astype('uint8'),cmap='gray',vmax = vmax,interpolation='nearest', aspect='auto')

    plt.axis('off')
    # plt.xlim(0, )
    plt.tight_layout()
    # plt.ylim(0,7)

    # Save the plot
    plot_filename = f"{file_stem}-cycle{kymo_cycle}-cropped DNAp Trace.png"
    plot_output_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_output_path, dpi=300,bbox_inches='tight')
    # plt.show()
    plt.close()

    plt.figure(figsize=(8,3))
    plt.plot(x1,all_intensity)
    plt.plot(x1,signal_filter)
    # plt.show()
    plt.close()

    # plot the DNAP intensity in a raw data, filtered data
    plt.figure(figsize=(8,3))
    # plt.plot(rec_x_list,bagrnd_signal_filter,linestyle='solid',linewidth=1, markersize=1,label='Background Signal Filter')
    plt.plot(x1,all_intensity,linestyle='solid',linewidth=0.3, markersize=0.5,label='Protein Trace Signal')
    plt.plot(x1,signal_filter,linestyle='solid',linewidth=1, markersize=0.5,label='Protein Trace Signal Filter')
    # plt.plot(x1,data_stepfinder['FinalFit'],color='black',linewidth=1, label='Step Fit')
    # plt.axhline(np.average(bagrnd_signal_filter),color = 'b',label='Background Signal Filter Mean')
    plt.axhline(np.average(bagrnd_signal_filter)+threshold_sigma*np.std(bagrnd_signal_filter), color='grey', linestyle='--',label=str(threshold_sigma) + ' Sigma threshold')
    plt.axhline(np.average(bagrnd_signal_filter)-threshold_sigma*np.std(bagrnd_signal_filter), color='grey', linestyle='--')

    plt.xlabel("Time/s")
    plt.ylabel("Intensity/photon")
    plt.title("Intensity along DNA polymerase Trajectory")
    plt.legend()
    plt.tight_layout()
    plt.xlim((kymo_xlim_left-x_offset_optimal)/x_cali,(kymo_xlim_right-x_offset_optimal)/x_cali)
    # plt.xlim(32,40)
    # plt.ylim(3500,6000)
    # plt.ylim((0,60))
    # plt.show()
    plt.close()

    # plt.savefig(kymo_filename.replace('.tdms', '-cycle') + kymo_cycle + 'Intensity along DNA polymerase Trajectory-1'+'.eps', format='eps', dpi=300, bbox_inches='tight')

    # Save the plot
    plot_filename = f"{file_stem}-cycle{kymo_cycle}-Intensity along DNA polymerase Trajectory-1.png"
    plot_output_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_output_path, dpi=300,bbox_inches='tight')
    # plt.show()
    plt.close()

    print(np.average(bagrnd_signal_filter)+threshold_sigma*np.std(bagrnd_signal_filter))

    # plot the DNAP intensity in a raw data, filtered data and overlap the basepair-time
    basepair = pd.read_excel(trace_file)
    bp_time = basepair['time']/1000
    bp_time = bp_time.dropna()
    bp = basepair['basepairs']
    bp = bp.dropna()

    fig, ax1 = plt.subplots(figsize=(8,3))
    font = {'family': 'Arial', 'weight': 'normal', 'size': 16}

    ax2 = ax1.twinx()
    ax1.plot(bp_time,bp,color='black',linewidth=1, label='basepairs')
    ax1.set_xlabel('Time/s',fontdict=font)
    ax1.set_ylabel('Basepairs', color='green',fontdict=font)

    ax2.plot(x1,all_intensity,color='gray',linestyle='solid',linewidth=0.2, markersize=0.5,label='Protein Trace Signal')
    ax2.plot(x1,signal_filter,color='orange',linestyle='solid',linewidth=1, markersize=0.5,label='Protein Trace Signal Filter')
    ax2.axhline(np.average(bagrnd_signal_filter)+threshold_sigma*np.std(bagrnd_signal_filter), color='grey', linestyle='--',label=str(threshold_sigma) + ' Sigma threshold')
    ax2.axhline(np.average(bagrnd_signal_filter)-threshold_sigma*np.std(bagrnd_signal_filter), color='grey', linestyle='--')

    ax2.set_xlabel("Time/s",fontdict=font)
    ax2.set_ylabel("Intensity/photon", color='black',fontdict=font)


    plt.xlim((kymo_xlim_left-x_offset_optimal)/x_cali,(kymo_xlim_right-x_offset_optimal)/x_cali)

    # plt.legend()
    # plt.xlim(20,26)

    # plt.title("Intensity along DNA polymerase Trajectory")
    plt.tight_layout()

    # plt.savefig(kymo_filename.replace('.tdms', '-cycle') + kymo_cycle + 'Intensity along DNA polymerase Trajectory'+'.eps', format='eps', dpi=300, bbox_inches='tight')

    # Save the plot
    plot_filename = f"{file_stem}-cycle{kymo_cycle}-Intensity along DNA polymerase Trajectory.png"
    plot_output_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_output_path, dpi=300,bbox_inches='tight')
    # plt.show()
    plt.close(fig)

    print(np.average(bagrnd_signal_filter)+threshold_sigma*np.std(bagrnd_signal_filter))
    # plt.close()

    # this data block is used to save the intensity data for step-finding analysis
    df2 = pd.DataFrame({'time serials': x1, 'filtered intensity': signal_filter})
    intensity_filtered_path = os.path.join(output_dir, f"{file_stem}-cycle{kymo_cycle}-Intensity_along_DNAp-filtered.txt")
    np.savetxt(intensity_filtered_path, df2.values, fmt='%1.3f')

    # This step intends to detect the steps of the fluroescence intensity along DNAp trace
    # Using AutoStepfinder: https://doi.org/10.1016/j.patter.2021.100256.
    print("\n--- Starting Intensity Step-Finding ---")

    # Make sure the autostepfinder module is in the Python path
    sys.path.append(os.path.join(project_root, 'autostepfinder_GUI'))
    from autostepfinder import AutoStepFinder, StepFinderParameters

    # Data is already in memory
    time_data = x1
    signal_data = signal_filter

    # Create Output Directory and Filename Stem for step-finding results
    output_dir_steps_rel = os.path.join('example_dataset/image data', 'stepfinding_results/')
    output_dir_steps = os.path.join(project_root, output_dir_steps_rel)
    os.makedirs(output_dir_steps, exist_ok=True)
    file_stem_step = os.path.splitext(os.path.basename(intensity_filtered_path))[0]
    print(f"Step-finding results will be saved in: '{output_dir_steps}/'")
    print(f"Output file prefix: '{file_stem_step}'")

    # Set Analysis Parameters for step finding
    time_resolution = np.mean(np.diff(time_data))
    print(f"Detected time resolution: {time_resolution:.4f} s/point")

    params = StepFinderParameters(
        s_max_threshold=0.15,
        resolution=time_resolution,
        fit_range=10000,
        fit_mode='mean',
        local_step_merge=True,
        error_tolerance=2.0,
        overshoot=1
    )

    # Run the Analysis
    print("\nRunning AutoStepFinder analysis...")
    finder = AutoStepFinder(params=params)
    final_fit, final_steps, s_curves, n_found_steps = finder.run(signal_data)
    print(f"Analysis complete. Found {len(final_steps)} steps.")

    # Save Results to Files
    print("\nSaving step-finding results...")
    trace_df = pd.DataFrame({'time': time_data, 'data': signal_data, 'fit': final_fit})
    trace_filename = os.path.join(output_dir_steps, f"{file_stem_step}_fit_trace.csv")
    trace_df.to_csv(trace_filename, index=False, float_format='%.4f')
    print(f"  - Fit trace saved to: {trace_filename}")

    if not final_steps.empty:
        properties_filename = os.path.join(output_dir_steps, f"{file_stem_step}_step_properties.csv")
        final_steps.to_csv(properties_filename, index=False, float_format='%.4f')
        print(f"  - Step properties saved to: {properties_filename}")
        print("\nStep Properties Table:")
        print(final_steps)

    # Visualize and Save Plots
    # Note: Saving plots as images requires the 'kaleido' package (pip install kaleido)
    print("\nGenerating and saving step-finding plots...")

    # Fit Plot
    fig_fit = go.Figure()
    fig_fit.add_trace(go.Scatter(x=time_data, y=signal_data, mode='lines', name='Original Data', line=dict(color='lightgray', width=1.5)))
    fig_fit.add_trace(go.Scatter(x=time_data, y=final_fit, mode='lines', name='Step Fit', line=dict(color='black', width=2)))
    fig_fit.update_layout(title="AutoStepFinder Analysis", xaxis_title=f"Time (s)", yaxis_title="Intensity", hovermode="x unified", height=400, font=dict(family="serif", size=14))
    fit_plot_filename = os.path.join(output_dir_steps, f"{file_stem_step}_fit_plot.png")
    fig_fit.write_image(fit_plot_filename, scale=3, width=1000, height=400)
    print(f"  - Fit plot saved to: {fit_plot_filename}")

    # Dwell Time Distribution
    if not final_steps.empty:
        dwell_times_sec = final_steps['dwell_time_after'] * params.resolution
        fig_dwell = px.histogram(x=dwell_times_sec, nbins=20, title="Dwell Time Distribution (After Step)", labels={'x': 'Dwell Time (s)'})
        fig_dwell.update_layout(yaxis_title="Count", height=400, font=dict(family="serif", size=14))
        dwell_plot_filename = os.path.join(output_dir_steps, f"{file_stem_step}_dwell_time_dist.png")
        fig_dwell.write_image(dwell_plot_filename, scale=3, width=800, height=400)
        print(f"  - Dwell time plot saved to: {dwell_plot_filename}")
    else:
        print("  - No steps found, skipping dwell time plot.")

    # Plotting Raw Data, Step-Fit, and Threshold
    fig_step_threshold = plt.figure(figsize=(8, 3))
    font_plot = {'family': 'serif', 'weight': 'normal', 'size': 14}
    plt.plot(time_data, signal_data, linestyle='solid', linewidth=0.7, label='Raw Signal', color='lightgray')
    plt.plot(time_data, final_fit, color='black', linewidth=1.5, label='Step Fit')
    plt.axhline(np.average(bagrnd_signal_filter)+threshold_sigma*np.std(bagrnd_signal_filter), color='grey', linestyle='--',label=str(threshold_sigma) + ' Sigma threshold')
    plt.axhline(np.average(bagrnd_signal_filter)-threshold_sigma*np.std(bagrnd_signal_filter), color='grey', linestyle='--')
    plt.xlabel("Time (s)", fontdict=font_plot)
    plt.ylabel("Intensity", fontdict=font_plot)
    plt.title("Intensity Trace with Step Fit", fontdict=font_plot)
    plt.legend()
    plt.tight_layout()
    plot_filename = os.path.join(output_dir_steps, f"{file_stem_step}_intensity_vs_fit.png")
    plt.savefig(plot_filename, dpi=300,bbox_inches='tight')
    print(f"Plot saved to: {plot_filename}")
    plt.close(fig_step_threshold)

    # Binarize Intensity Based on Threshold
    threshold = np.mean(bagrnd_signal_filter) + (threshold_sigma + 1.1) * np.std(bagrnd_signal_filter)
    signal_binarized = np.where(final_fit > threshold, 1, 0)
    fig_binarized = plt.figure(figsize=(8, 3))
    plt.plot(time_data, signal_binarized, color='black', linewidth=1.5)
    plt.xlabel("Time (s)", fontdict=font_plot)
    plt.ylabel("Binarized State", fontdict=font_plot)
    plt.title("Binarized Intensity Trace", fontdict=font_plot)
    plt.ylim(-0.1, 1.1)
    plt.tight_layout()
    plot_filename = os.path.join(output_dir_steps, f"{file_stem_step}_binarized_trace.png")
    plt.savefig(plot_filename, dpi=300)
    print(f"Plot saved to: {plot_filename}")
    plt.close(fig_binarized)

    # Filter Short-Lived Events from Binarized Trace
    def filter_short_events(binary_trace, min_length_points):
        """Removes consecutive runs of 0s or 1s shorter than min_length_points."""
        filtered_trace = np.copy(binary_trace)
        for state_to_filter in [1, 0]:
            padded_trace = np.pad(filtered_trace, 1, mode='constant', constant_values=1 - state_to_filter)
            diff = np.diff(padded_trace)
            starts = np.where(diff == (1 if state_to_filter == 1 else -1))[0]
            ends = np.where(diff == (-1 if state_to_filter == 1 else 1))[0]
            for start, end in zip(starts, ends):
                if end - start < min_length_points:
                    filtered_trace[start:end] = 1 - state_to_filter
        return filtered_trace

    min_event_duration_points = 3
    signal_binarized_filtered = filter_short_events(signal_binarized, min_event_duration_points)
    fig_bin_filtered = plt.figure(figsize=(8, 3))
    plt.plot(time_data, signal_binarized_filtered, label="Filtered Binarized Trace", color='black', linewidth=1.5)
    plt.xlabel("Time (s)", fontdict=font_plot)
    plt.ylabel("Binarized State", fontdict=font_plot)
    plt.title(f"Binarized Trace (Events < {min_event_duration_points} points removed)", fontdict=font_plot)
    plt.ylim(-0.1, 1.1)
    plt.tight_layout()
    plot_filename = os.path.join(output_dir_steps, f"{file_stem_step}_binarized_filtered.png")
    plt.savefig(plot_filename, dpi=300)
    print(f"Plot saved to: {plot_filename}")
    plt.close(fig_bin_filtered)

    # Kymograph with dual time axis (pixels and seconds)
    fig_kymo_dual_axis = plt.figure(figsize=(8, 3))
    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(bottom=0.2)
    par = host.twiny()
    offset = -30
    new_fixed_axis = host.get_grid_helper().new_fixed_axis
    par.axis["bottom"] = new_fixed_axis(loc="bottom", axes=par, offset=(0, offset))
    par.axis["top"].set_visible(False)
    host.axis["right"].set_visible(False)
    host.imshow(g.astype('uint16'), cmap='gray', vmax=vmax, aspect="auto")
    host.set_xlabel("Time(pixel)", fontdict=font_plot)
    host.set_ylabel("Position(pixel)", fontdict=font_plot)
    host.set_xlim((kymo_xlim_left, kymo_xlim_right))
    host.set_ylim((kymo_ylim_bottom, kymo_ylim_top))
    par.set_xlabel("Time(s)", fontdict=font_plot)
    par.set_xlim((kymo_xlim_left - x_offset_optimal) / x_cali, (kymo_xlim_right - x_offset_optimal) / x_cali)
    host.plot(trace_time, position, 'lime', linewidth=2)
    plt.tight_layout()
    plot_filename = os.path.join(output_dir_steps, f"{file_stem_step}_DNA_Polymerase_Trace_2.png")
    plt.savefig(plot_filename, dpi=300)
    print(f"Kymograph plot saved to: {plot_filename}")
    plt.close(fig_kymo_dual_axis)

    # Triple-axis plot: Basepairs, Step-like Intensity, Binarized Intensity
    fig_triple = plt.figure(figsize=(8, 3))
    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)
    par1 = host.twinx()
    par2 = host.twinx()
    par3 = host.twinx()
    offset = 25
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))
    par2.axis["right"].toggle(all=True)
    host.set_xlim((kymo_xlim_left - x_offset_optimal) / x_cali, (kymo_xlim_right - x_offset_optimal) / x_cali)
    host.set_xlabel("Time/s", fontdict=font_plot)
    host.set_ylabel("Basepairs/bp", fontdict=font_plot)
    par1.set_ylabel("Step-like Intensity", fontdict=font_plot)
    par2.set_ylabel("Binarized Intensity", fontdict=font_plot)
    p0, = host.plot(bp_time, bp, "lime", label='Basepairs', linewidth=1.5)
    p1, = par1.plot(x1, final_fit, color='black', linewidth=1, label='Step-like intensity')
    p2, = par2.plot(x1, signal_binarized_filtered, color='lightgray', linewidth=0.1, label='Binarized Intensity')
    par2.fill_between(x1, signal_binarized_filtered, 0, color='lightgreen', alpha=0.1)
    p3, = par3.plot(x1, all_intensity, color='lightgray', linewidth=0.6, linestyle='solid', markersize=0.5, label='Raw Intensity')
    par3.set_yticklabels([])
    par3.set_yticks([])
    par1.axhline(np.average(bagrnd_signal_filter) + threshold_sigma * np.std(bagrnd_signal_filter), color='grey', linestyle='--')
    par1.axhline(np.average(bagrnd_signal_filter) - threshold_sigma * np.std(bagrnd_signal_filter), color='grey', linestyle='--')
    par1.set_ylim(7, 87)
    par2.set_ylim(0, 1)
    host.axis["left"].label.set_color('green')
    par1.axis["right"].label.set_color(p1.get_color())
    par2.axis["right"].label.set_color(p1.get_color())
    plot_filename = os.path.join(output_dir_steps, f"{file_stem_step}_all_correlated.png")
    plt.savefig(plot_filename, dpi=300)
    print(f"All-correlated plot saved to: {plot_filename}")
    plt.close(fig_triple)

    # --- Basepair Change-Point Detection ---
    # This section analyzes the basepair trace to find periods of activity and pausing.
    # Ensure the ChangePointDetection module is in the Python path.
    sys.path.append(os.path.join(project_root, 'ChangePointDetection_slope_GUI'))
    from bp_detection.bp_batch_segments import bp_batch_segments_for

    # Load file
    file_path_rel = 'example_dataset/force data/results/OTdata_example_30nM DNAp + trx  + 625uM dNTPs-cycle#01-processedData.xlsx'
    file_path = os.path.join(project_root, file_path_rel)

    # Create Output Directory and Filename Stem
    folder_save_rel = os.path.join('example_dataset/image data', 'ChangePoints_Results/')
    folder_save = os.path.join(project_root, folder_save_rel)
    os.makedirs(folder_save, exist_ok=True)
    base_filename = os.path.splitext(os.path.basename(file_path))[0]

    # Analysis Parameters (Tunable)
    window_size = 9      # Smaller window_size gives more segments (default = 6)
    sigma_noise = 0.08   # Smaller sigma_noise gives more segments (default = 0.06)

    # Column Names 
    time_col = 'time'
    signal_col = 'basepairs'

    df = pd.read_excel(file_path)
    data = np.zeros((len(df), 5))
    data[:, 0] = np.arange(1, len(df) + 1)
    data[:, 1] = df[time_col].values
    data[:, 4] = df[signal_col].values

    # --- 1. Calculate First Derivative with a Moving-Window Filter ---
    print("--- Starting Change-Point Detection ---")
    velocity_par = []
    time_s = []
    nbpts = window_size

    # The loop for calculating the derivative requires the data length to be > 2 * nbpts.
    if len(data) > 2 * nbpts:
        for it in range(nbpts, len(data[:, 0]) - nbpts):
            window_indices = np.arange(it - nbpts, it + nbpts + 1).astype(int)
            window_x = data[window_indices, 1]  # Time values from window
            window_y = data[window_indices, 4]  # Signal values from window
            
            # Linear regression on the window to find the slope
            p = np.polyfit(window_x, window_y, 1)
            velocity_par.append(p[0])
            time_s.append(data[it, 1])
    else:
        print(f"Warning: Data length ({len(data)}) is too short for the window size ({nbpts}).")
        print("The derivative cannot be calculated. Consider reducing 'window_size'.")

    # --- 2. Detect Change-Points ---
    velocity_par_all = np.zeros((len(data[:, 0]), 2))
    segments = []  # Initialize segments as an empty list

    if velocity_par:  # Only proceed if the derivative calculation was successful
        velocity_par_all[:nbpts, 1] = velocity_par[0]
        velocity_par_all[nbpts:nbpts+len(time_s), 1] = velocity_par
        velocity_par_all[nbpts+len(time_s):, 1] = velocity_par[-1]
        velocity_par_all[:, 0] = data[:, 0]
        
        sigma_here = sigma_noise
        segments, CPs = bp_batch_segments_for(
            [data[:, 0]], 
            [velocity_par_all[:, 1]], 
            sigma=sigma_here,
            linear=False
        )

    # --- 3. Process and Plot Results ---
    # Define file paths for outputs
    excel_filename = os.path.join(folder_save, f"{base_filename}-correlated_analysis.xlsx")
    bp_time = df[time_col].values
    bp_signal = df[signal_col].values

    if segments and segments[0].size > 0:
        segments_array = segments[0]
        
        output_df = pd.DataFrame({
            'cp_startIndex': segments_array[:, 0],
            'cp_endIndex': segments_array[:, 2]
        })
        print(f"Change-point detection complete. Found {len(output_df)} segments.")

        fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

        # Plot 1: Original trace
        axs[0].plot(data[:, 0], data[:, 4], '.k', markersize=2)
        axs[0].set_title('Original Data')
        axs[0].set_ylabel(signal_col)

        # Plot 2: First derivative with detected segments
        axs[1].plot(data[:, 0], velocity_par_all[:, 1], '-.k', linewidth=0.5)
        for i in range(segments_array.shape[0]):
            axs[1].plot(segments_array[i, [0, 2]], segments_array[i, [1, 3]], linewidth=3)
        axs[1].set_title('First Derivative with Segments')
        axs[1].set_ylabel('Slope (1st Derivative)')

        # Plot 3: Original trace with fitted segments
        segments_h = segments_array[:, [0, 2]].T
        axs[2].plot(data[:, 0], data[:, 4], '.k', markersize=2)
        for iseg in range(segments_h.shape[1]):
            # FIX: Cast indices to integer before using them to index the numpy array.
            start_idx = int(segments_h[0, iseg] - 1)
            end_idx = int(segments_h[1, iseg] - 1)
            
            if 0 <= start_idx < len(data) and 0 <= end_idx < len(data):
                axs[2].plot(data[[start_idx, end_idx], 0], data[[start_idx, end_idx], 4], linewidth=3)
                
        axs[2].set_title('Original Data with Fitted Segments')
        axs[2].set_ylabel(signal_col)
        axs[2].set_xlabel('Time Index')

        plt.tight_layout()
        fig_filename = os.path.join(folder_save, f"{base_filename}_WinSize_{nbpts}_Sigma_{sigma_here}_summary.png")
        plt.savefig(fig_filename, dpi=150)
        print(f"Summary plot saved to {fig_filename}")
        # plt.show()
        plt.close(fig)

        # --- Data Preparation for subsequent plots/analysis ---
        start_indices = (output_df['cp_startIndex'] - 1).astype(int)
        end_indices = (output_df['cp_endIndex'] - 1).astype(int)
        cp_startTime = bp_time[start_indices]
        cp_endTime = bp_time[end_indices]
        cp_startBasepair = bp_signal[start_indices]
        cp_endBasepair = bp_signal[end_indices]

        # --- Interactive Plotly Plot ---
        print("\n--- Generating Interactive Step-Fit Plot with Plotly ---")
        fig_interactive = go.Figure()
        fig_interactive.add_trace(go.Scatter(
            x=bp_time, y=bp_signal, mode='lines',
            line=dict(color='lightgray', width=1.5), name='Raw Data'
        ))
        for i in range(len(cp_startTime)):
            fig_interactive.add_trace(go.Scatter(
                x=[cp_startTime[i], cp_endTime[i]], y=[cp_startBasepair[i], cp_endBasepair[i]],
                mode='lines', line=dict(width=3), name=f'Segment {i + 1}'
            ))
        fig_interactive.update_layout(
            title="Interactive ChangePoints for Basepairs vs. Time",
            xaxis_title="Time (s)", yaxis_title=signal_col,
            legend_title="Trace", template="plotly_white"
        )
        html_fig_filename = os.path.join(folder_save, f"{base_filename}-InteractiveStepPlot.html")
        fig_interactive.write_html(html_fig_filename)
        print(f"Interactive plot saved to {html_fig_filename}")
        # fig_interactive.show()

        # --- Interpolated Data Plot ---
        print("\n--- Generating Interpolated Data Plot ---")
        cp_Time = np.append(cp_startTime, cp_endTime[-1])
        cp_Basepair = np.append(cp_startBasepair, cp_endBasepair[-1])
        interp_df = pd.DataFrame({"Time": cp_Time, "Basepair": cp_Basepair}).drop_duplicates()
        set_interp = interp1d(interp_df['Time'], interp_df['Basepair'], kind='linear', fill_value="extrapolate")
        cp_basepair_interp = set_interp(bp_time)

        plt.figure(figsize=(8, 3))
        plt.plot(bp_time, cp_basepair_interp, "red", label="Step-Fitted Data", linewidth=2)
        plt.plot(bp_time, bp_signal, "lightgrey", label="Raw Data", linewidth=1)
        plt.legend()
        plt.xlabel("Time (s)")
        plt.ylabel("Interpolated Basepairs")
        plt.title("Activity Bursts as a Function of Time")
        plt.tight_layout()
        interp_fig_filename = os.path.join(folder_save, f"{base_filename}-InterpolatedPlot.png")
        plt.savefig(interp_fig_filename, dpi=150)
        print(f"Interpolation plot saved to {interp_fig_filename}")
        # plt.show()
        plt.close()

        # --- Saving All Analyzed Data to Excel ---
        print("\n--- Saving All Analyzed Data to Excel ---")
        with pd.ExcelWriter(excel_filename) as writer:
            # Sheet 1: Raw data used in analysis
            df1 = pd.DataFrame({'time_s': bp_time, 'raw_basepair': bp_signal})
            df1.to_excel(writer, sheet_name='raw_data', index=False)
            
            # Sheet 2: Detected change-points with time and basepair values
            df2 = pd.DataFrame({
                'cp_startTime_s': cp_startTime,
                'cp_startBasepair': cp_startBasepair,
                'cp_endTime_s': cp_endTime,
                'cp_endBasepair': cp_endBasepair,
            })
            df2.to_excel(writer, sheet_name='change_points', index=False)
            
            # Sheet 3: Interpolated data based on change-points
            df3 = pd.DataFrame({'time_s': bp_time, 'interpolated_basepair': cp_basepair_interp})
            df3.to_excel(writer, sheet_name='interpolated_data', index=False)

            df4 = pd.DataFrame({'time_s':x1, 'step_intensity':final_fit})
            df4.to_excel(writer, sheet_name ='step_intensity', index=False)
    
            df5 = pd.DataFrame({'time_s':x1, 'binarized_intensity':signal_binarized_filtered})
            df5.to_excel(writer, sheet_name ='binarized_intensity', index=False)
        print(f"Combined analysis data saved to {excel_filename}")

    else:
        print("\nNo basepair change-points were detected.")
        print("Consider adjusting 'sigma_noise' or 'window_size'.")
        print("Saving raw data to Excel.")
        with pd.ExcelWriter(excel_filename) as writer:
            df1 = pd.DataFrame({'time_s': bp_time, 'raw_basepair': bp_signal})
            df1.to_excel(writer, sheet_name='raw_data', index=False)
        print(f"Raw data saved to {excel_filename}")
        
    print("\nAnalysis complete.")

    # Save all the analysis parameters in the analysis log in .txt format
    analysis_log_data = {"kymo_cycle": kymo_cycle,
                         "kymo_xlim_left": kymo_xlim_left,
                         "kymo_xlim_right": kymo_xlim_right,
                         "x_offset": x_offset_optimal,
                         "y_offset": y_offset_optimal,
                         "rec_x": rec_x,
                         "rec_y": rec_y,
                         "rec_w": rec_w,
                         "rec_h": rec_h,
                         "bagrnd_filter_window": bagrnd_filter_window,
                         "signal_filter_window": signal_filter_window,
                         "threshold_sigma": threshold_sigma,
                         "Intensity_changepoint#": len(final_steps),
                         "windowsize_cp_basepair": window_size,
                         "sigmavalue_basepair": sigma_noise}

    # The original notebook code contained comments about MATLAB, which are not relevant here.

    log_filename = excel_filename.replace('.xlsx', '-analysis_log_data.txt')
    with open(log_filename, 'w') as text_file:
        # Saving as a more readable format (like JSON) could be a future improvement.
        for key, value in analysis_log_data.items():
            text_file.write(f'{key}: {value}\n')
        
    print(f"Analysis log saved to {log_filename}")


if __name__ == "__main__":
    main()

