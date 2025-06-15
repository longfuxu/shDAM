import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import Cursor
import pandas as pd
from nptdms import TdmsFile
from sympy import coth
from fast_pwl_fit import FastPWLFit
from scipy.signal import savgol_filter
import matplotlib as mpl

# Define numerical coth function for backup or direct calculation
def numerical_coth(x):
    """Numerical implementation of hyperbolic cotangent for arrays"""
    return 1.0 / np.tanh(np.array(x, dtype=float))

class OTDataAnalyzer:
    def __init__(self, root):
        self.root = root
        self.root.title("Optical Tweezers Data Analyzer")
        self.root.geometry("1400x900")
        
        # Configure root grid weights for scaling
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_columnconfigure(1, weight=3)
        
        # Variables for parameters
        self.filename_var = tk.StringVar()
        self.cycle_var = tk.StringVar(value="01")
        self.time_from_var = tk.StringVar()
        self.time_to_var = tk.StringVar()
        self.time_from_exo_var = tk.StringVar()
        self.time_to_exo_var = tk.StringVar()
        self.time_from_pol_var = tk.StringVar()
        self.time_to_pol_var = tk.StringVar()
        self.bead_size_var = tk.StringVar(value="1.76")
        self.segment_number_var = tk.StringVar(value="20")  # New variable for segment number
        self.ssb_factor_var = tk.StringVar(value="0.03")  # SSB factor for 10pN
        self.total_basepairs_var = tk.StringVar(value="8393")  # Total basepairs
        
        # Model parameters - updated to match the notebook values
        self.C_var = tk.StringVar(value="440")
        self.g0_var = tk.StringVar(value="-637")
        self.g1_var = tk.StringVar(value="17")
        self.Lc_var = tk.StringVar(value="2.85056")
        self.Lp_var = tk.StringVar(value="56")
        self.S_var = tk.StringVar(value="1500")
        self.Lss_var = tk.StringVar(value="4.69504")
        self.b_var = tk.StringVar(value="1.5")
        self.Sss_var = tk.StringVar(value="800")
        
        # Force parameters for exo and pol phases
        self.exo_force_var = tk.StringVar(value="50")
        self.pol_force_var = tk.StringVar(value="20")
        
        # Data storage
        self.time_data = None
        self.force_data = None
        self.distance_data = None
        self.active_time_entry = None
        
        # Data storage for exo and pol phases
        self.time_range_exo = None
        self.force_range_exo = None
        self.distance_range_exo = None
        self.time_range_pol = None
        self.force_range_pol = None
        self.distance_range_pol = None
        
        self.create_gui()
        
    def create_gui(self):
        # Create main frames with proper weights
        main_frame = ttk.Frame(self.root)
        main_frame.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=5, pady=5)
        main_frame.grid_columnconfigure(1, weight=3)
        main_frame.grid_rowconfigure(0, weight=1)
        
        input_frame = ttk.LabelFrame(main_frame, text="Input Parameters", padding="10")
        input_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
        
        plot_frame = ttk.LabelFrame(main_frame, text="Plots", padding="10")
        plot_frame.grid(row=0, column=1, padx=5, pady=5, sticky="nsew")
        plot_frame.grid_rowconfigure(0, weight=1)
        plot_frame.grid_columnconfigure(0, weight=1)
        
        # File selection
        file_frame = ttk.Frame(input_frame)
        file_frame.grid(row=0, column=0, columnspan=3, sticky="ew", pady=5)
        ttk.Label(file_frame, text="TDMS File:").pack(side=tk.LEFT)
        ttk.Entry(file_frame, textvariable=self.filename_var, width=40).pack(side=tk.LEFT, padx=5)
        ttk.Button(file_frame, text="Browse", command=self.browse_file).pack(side=tk.LEFT)
        
        # Basic parameters
        params_frame = ttk.Frame(input_frame)
        params_frame.grid(row=1, column=0, columnspan=3, sticky="ew", pady=5)
        
        ttk.Label(params_frame, text="Cycle Number:").grid(row=0, column=0, sticky="w", pady=2)
        ttk.Entry(params_frame, textvariable=self.cycle_var, width=10).grid(row=0, column=1, sticky="w", padx=5)
        
        # Add total basepairs input in the first row on the right
        ttk.Label(params_frame, text="Total Basepairs:").grid(row=0, column=2, sticky="w", pady=2)
        ttk.Entry(params_frame, textvariable=self.total_basepairs_var, width=10).grid(row=0, column=3, sticky="w", padx=5)
        
        # Time entries in two columns
        # Left column - Start times
        ttk.Label(params_frame, text="Exo Start Time (ms):").grid(row=1, column=0, sticky="w", pady=2)
        exo_start_entry = ttk.Entry(params_frame, textvariable=self.time_from_exo_var, width=10)
        exo_start_entry.grid(row=1, column=1, sticky="w", padx=5)
        exo_start_entry.bind("<FocusIn>", lambda e: self.set_active_time_entry("exo_start"))
        
        ttk.Label(params_frame, text="Pol Start Time (ms):").grid(row=2, column=0, sticky="w", pady=2)
        pol_start_entry = ttk.Entry(params_frame, textvariable=self.time_from_pol_var, width=10)
        pol_start_entry.grid(row=2, column=1, sticky="w", padx=5)
        pol_start_entry.bind("<FocusIn>", lambda e: self.set_active_time_entry("pol_start"))
        
        # Right column - End times
        ttk.Label(params_frame, text="Exo End Time (ms):").grid(row=1, column=2, sticky="w", pady=2)
        exo_end_entry = ttk.Entry(params_frame, textvariable=self.time_to_exo_var, width=10)
        exo_end_entry.grid(row=1, column=3, sticky="w", padx=5)
        exo_end_entry.bind("<FocusIn>", lambda e: self.set_active_time_entry("exo_end"))
        
        ttk.Label(params_frame, text="Pol End Time (ms):").grid(row=2, column=2, sticky="w", pady=2)
        pol_end_entry = ttk.Entry(params_frame, textvariable=self.time_to_pol_var, width=10)
        pol_end_entry.grid(row=2, column=3, sticky="w", padx=5)
        pol_end_entry.bind("<FocusIn>", lambda e: self.set_active_time_entry("pol_end"))
        
        # Force parameters
        ttk.Label(params_frame, text="Exo Force (pN):").grid(row=4, column=0, sticky="w", pady=2)
        ttk.Entry(params_frame, textvariable=self.exo_force_var, width=10).grid(row=4, column=1, sticky="w", padx=5)
        
        ttk.Label(params_frame, text="Pol Force (pN):").grid(row=4, column=2, sticky="w", pady=2)
        ttk.Entry(params_frame, textvariable=self.pol_force_var, width=10).grid(row=4, column=3, sticky="w", padx=5)
        
        ttk.Label(params_frame, text="Bead Size:").grid(row=5, column=0, sticky="w", pady=2)
        ttk.Entry(params_frame, textvariable=self.bead_size_var, width=10).grid(row=5, column=1, sticky="w", padx=5)
        
        # Add segment number input
        ttk.Label(params_frame, text="Segment Number:").grid(row=5, column=2, sticky="w", pady=2)
        ttk.Entry(params_frame, textvariable=self.segment_number_var, width=10).grid(row=5, column=3, sticky="w", padx=5)
        
        # Model parameters frame
        model_frame = ttk.LabelFrame(input_frame, text="Model Parameters", padding="5")
        model_frame.grid(row=2, column=0, columnspan=3, pady=10, sticky="nsew")
        
        # tWLC parameters
        ttk.Label(model_frame, text="tWLC Parameters:").grid(row=0, column=0, columnspan=2, sticky="w")
        params = [
            ("C:", self.C_var), ("g0:", self.g0_var), ("g1:", self.g1_var),
            ("Lc:", self.Lc_var), ("Lp:", self.Lp_var), ("S:", self.S_var)
        ]
        for i, (label, var) in enumerate(params):
            ttk.Label(model_frame, text=label).grid(row=i+1, column=0, sticky="w")
            ttk.Entry(model_frame, textvariable=var, width=10).grid(row=i+1, column=1, padx=5)
            
        # FJC parameters
        ttk.Label(model_frame, text="FJC Parameters:").grid(row=0, column=2, columnspan=2, sticky="w")
        fjc_params = [
            ("Lss:", self.Lss_var), ("b:", self.b_var), ("Sss:", self.Sss_var)
        ]
        for i, (label, var) in enumerate(fjc_params):
            ttk.Label(model_frame, text=label).grid(row=i+1, column=2, sticky="w")
            ttk.Entry(model_frame, textvariable=var, width=10).grid(row=i+1, column=3, padx=5)
        
        # Buttons frame with grid layout
        buttons_frame = ttk.Frame(input_frame)
        buttons_frame.grid(row=3, column=0, columnspan=3, pady=10, sticky="ew")
        buttons_frame.grid_columnconfigure(0, weight=1)
        buttons_frame.grid_columnconfigure(1, weight=1)
        
        ttk.Button(buttons_frame, text="Fit to Model", command=self.fit_to_model).grid(row=0, column=0, columnspan=2, sticky="ew", pady=2, padx=2)
        
        # Create a frame for the tracking buttons
        track_buttons_frame = ttk.Frame(buttons_frame)
        track_buttons_frame.grid(row=1, column=0, columnspan=2, sticky="ew", pady=2)
        track_buttons_frame.grid_columnconfigure(0, weight=1)
        track_buttons_frame.grid_columnconfigure(1, weight=1)
        
        ttk.Button(track_buttons_frame, text="Track ss/dsDNA Junction (forward)", 
                  command=lambda: self.track_junction(reverse=False, use_exo_pol=True)).grid(row=0, column=0, sticky="ew", padx=2)
        ttk.Button(track_buttons_frame, text="Track ss/dsDNA Junction (reverse)", 
                  command=lambda: self.track_junction(reverse=True, use_exo_pol=True)).grid(row=0, column=1, sticky="ew", padx=2)
        
        # Add a button to save DNA polymerase trace
        ttk.Button(buttons_frame, text="Save DNA Polymerase Trace", 
                  command=self.plot_dna_polymerase_trace).grid(row=2, column=0, columnspan=2, sticky="ew", pady=2, padx=2)
        
        # Add Linear Segment Fitting button
        ttk.Button(buttons_frame, text="Linear Segment Fitting", 
                  command=self.linear_segment_fitting).grid(row=3, column=0, columnspan=2, sticky="ew", pady=2, padx=2)
        
        ttk.Button(buttons_frame, text="Save Data", 
                  command=self.save_data).grid(row=4, column=0, columnspan=2, sticky="ew", pady=2, padx=2)
        
        # Configure plot frame
        self.fig = plt.Figure(figsize=(12, 8))
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        
        # Add toolbar
        toolbar_frame = ttk.Frame(plot_frame)
        toolbar_frame.grid(row=1, column=0, sticky="ew")
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()
        
    def set_active_time_entry(self, entry_type):
        self.active_time_entry = entry_type
        
    def browse_file(self):
        filename = filedialog.askopenfilename(filetypes=[("TDMS files", "*.tdms")])
        if filename:
            self.filename_var.set(filename)
            self.load_and_plot_data()
            
    def on_click(self, event):
        if event.inaxes and self.active_time_entry:
            if self.active_time_entry == "start":
                self.time_from_var.set(f"{event.xdata:.1f}")
            elif self.active_time_entry == "end":
                self.time_to_var.set(f"{event.xdata:.1f}")
            elif self.active_time_entry == "exo_start":
                self.time_from_exo_var.set(f"{event.xdata:.1f}")
            elif self.active_time_entry == "exo_end":
                self.time_to_exo_var.set(f"{event.xdata:.1f}")
            elif self.active_time_entry == "pol_start":
                self.time_from_pol_var.set(f"{event.xdata:.1f}")
            elif self.active_time_entry == "pol_end":
                self.time_to_pol_var.set(f"{event.xdata:.1f}")
            
    def load_and_plot_data(self):
        try:
            # Read TDMS file
            tdms_file = TdmsFile(self.filename_var.get())
            
            # Extract data
            self.time_data = np.array(tdms_file['FD Data']['Time (ms)'][:])
            self.force_data = np.array(tdms_file['FD Data']['Force Channel 0 (pN)'][:])
            self.distance_data = np.array(tdms_file['FD Data']['Distance 1 (um)'][:])
            
            # Clear previous plots
            self.fig.clear()
            
            # Create initial force vs time plot
            ax = self.fig.add_subplot(111)
            ax.plot(self.time_data, self.force_data)
            ax.set_xlabel('Time (ms)')
            ax.set_ylabel('Force (pN)')
            ax.set_title('Click to select time points (select input box first)')
            
            # Add cursor
            cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
            
            # Connect click event
            self.canvas.mpl_connect('button_press_event', self.on_click)
            
            self.fig.tight_layout()
            self.canvas.draw()
            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while loading the file: {str(e)}")
            
    def tWLC(self, F):
        C = float(self.C_var.get())
        g0 = float(self.g0_var.get())
        g1 = float(self.g1_var.get())
        Lc = float(self.Lc_var.get())
        Lp = float(self.Lp_var.get())
        S = float(self.S_var.get())
        return Lc * (1 - 0.5 * (4.1/(F*Lp))**0.5 + C*F/(-(g0+g1*F)**2 + S*C))
    
    def FJC(self, F):
        Lss = float(self.Lss_var.get())
        b = float(self.b_var.get())
        Sss = float(self.Sss_var.get())
        EEDss = []
        for Fext in F:
            x = Lss * (float(coth(Fext * b / 4.1).evalf()) - 4.1 / (Fext * b)) * (1 + Fext / Sss)
            EEDss.append(x)
        return np.array(EEDss)
    
    def get_analysis_data(self):
        """Get analysis data for the ROI time range"""
        time_from = float(self.time_from_var.get())
        time_to = float(self.time_to_var.get())
        indtemp = np.where((self.time_data <= time_to) & (self.time_data >= time_from))
        
        time_range = self.time_data[indtemp]
        force_range = self.force_data[indtemp]
        distance_range = self.distance_data[indtemp]
        
        bead_size = float(self.bead_size_var.get())
        dsDNA_ref = self.tWLC(45)
        ssDNA_ref = float(self.Lss_var.get()) * (float(coth(45 * float(self.b_var.get()) / 4.1).evalf()) - 4.1 / (45 * float(self.b_var.get()))) * (1 + 45 / float(self.Sss_var.get()))
        
        ssDNA_percentage = (distance_range - bead_size - dsDNA_ref)/(ssDNA_ref - dsDNA_ref)
        total_basepairs = float(self.total_basepairs_var.get())
        basepairs = (1-ssDNA_percentage) * total_basepairs
        junction_position = (ssDNA_percentage * ssDNA_ref) * (distance_range - bead_size) / ((ssDNA_percentage * ssDNA_ref) + (1 - ssDNA_percentage) * dsDNA_ref)
        
        return time_range, force_range, distance_range, ssDNA_percentage, basepairs, junction_position
    
    def get_exo_pol_analysis_data(self):
        """Get separate analysis data for exo and pol phases"""
        try:
            # Get exo data
            time_from_exo = float(self.time_from_exo_var.get())
            time_to_exo = float(self.time_to_exo_var.get())
            indtemp_exo = np.where((self.time_data <= time_to_exo) & (self.time_data >= time_from_exo))
            
            self.time_range_exo = self.time_data[indtemp_exo]
            self.force_range_exo = self.force_data[indtemp_exo]
            self.distance_range_exo = self.distance_data[indtemp_exo]
            
            # Get pol data
            time_from_pol = float(self.time_from_pol_var.get())
            time_to_pol = float(self.time_to_pol_var.get())
            indtemp_pol = np.where((self.time_data <= time_to_pol) & (self.time_data >= time_from_pol))
            
            self.time_range_pol = self.time_data[indtemp_pol]
            self.force_range_pol = self.force_data[indtemp_pol]
            self.distance_range_pol = self.distance_data[indtemp_pol]
            
            # Combine for all data
            time_range_all = np.append(self.time_range_exo, self.time_range_pol)
            force_range_all = np.append(self.force_range_exo, self.force_range_pol)
            distance_range_all = np.append(self.distance_range_exo, self.distance_range_pol)
            
            # Calculate parameters
            bead_size = float(self.bead_size_var.get())
            ssb_factor = float(self.ssb_factor_var.get())
            exo_force = float(self.exo_force_var.get())
            pol_force = float(self.pol_force_var.get())
            
            # Calculate reference values for dsDNA and ssDNA under exo and pol forces
            dsDNA_exo_ref = self.tWLC(exo_force)
            dsDNA_pol_ref = self.tWLC(pol_force)
            
            # Calculate ssDNA reference values
            Lss = float(self.Lss_var.get())
            b = float(self.b_var.get())
            Sss = float(self.Sss_var.get())
            
            ssDNA_exo_ref = Lss * (float(coth(exo_force * b / 4.1).evalf()) - 4.1 / (exo_force * b)) * (1 + exo_force / Sss)
            ssDNA_pol_ref = Lss * (float(coth(pol_force * b / 4.1).evalf()) - 4.1 / (pol_force * b)) * (1 + pol_force / Sss)
            
            # Calculate ssDNA percentages
            ssDNA_exo_percentage = (self.distance_range_exo - bead_size - dsDNA_exo_ref)/(ssDNA_exo_ref - ssb_factor - dsDNA_exo_ref)
            ssDNA_pol_percentage = (self.distance_range_pol - bead_size - dsDNA_pol_ref)/(ssDNA_pol_ref - ssb_factor - dsDNA_pol_ref)
            ssDNA_all_percentage = np.append(ssDNA_exo_percentage, ssDNA_pol_percentage)
            
            # Calculate basepairs
            total_basepairs = float(self.total_basepairs_var.get())
            basepairs = (1 - ssDNA_all_percentage) * total_basepairs
            
            # Calculate junction positions for both exo and pol phases
            # We'll use the forward direction by default, the track_junction method will handle reverse if needed
            junction_position_exo = self.get_junction_position(
                self.distance_range_exo, 
                ssDNA_exo_percentage, 
                ssDNA_exo_ref, 
                dsDNA_exo_ref, 
                bead_size, 
                False  # Forward direction
            )
            
            junction_position_pol = self.get_junction_position(
                self.distance_range_pol, 
                ssDNA_pol_percentage, 
                ssDNA_pol_ref, 
                dsDNA_pol_ref, 
                bead_size, 
                False  # Forward direction
            )
            
            # Combine junction positions
            junction_position_all = np.append(junction_position_exo, junction_position_pol)
            
            return time_range_all, force_range_all, distance_range_all, ssDNA_all_percentage, basepairs, junction_position_all
            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during exo/pol analysis: {str(e)}")
            return None, None, None, None, None, None
    
    def fit_to_model(self):
        try:
            if self.time_data is None:
                messagebox.showerror("Error", "Please load a TDMS file first")
                return
            
            # Use exo/pol data instead of ROI data
            if (self.time_from_exo_var.get() and self.time_to_exo_var.get() and 
                self.time_from_pol_var.get() and self.time_to_pol_var.get()):
                time_range, force_range, distance_range, ssDNA_percentage, basepairs, _ = self.get_exo_pol_analysis_data()
                if time_range is None:  # Error occurred in get_exo_pol_analysis_data
                    return
            else:
                messagebox.showerror("Error", "Please set Exo and Pol time ranges first")
                return
            
            # Clear previous plots
            self.fig.clear()
            
            # Create subplots in the top half
            gs = self.fig.add_gridspec(2, 2, height_ratios=[1, 0.1])
            ax1 = self.fig.add_subplot(gs[0, 0])
            ax2 = self.fig.add_subplot(gs[0, 1])
            
            # Plot 1: Models and experimental data
            Force = np.linspace(0.1, 68, 1000)
            ax1.plot(self.tWLC(Force), Force, 'r-', label='WLC Model')
            ax1.plot(self.FJC(Force), Force, 'b-', label='FJC Model')
            # ax1.plot(self.FJC_2_parallel_ssDNA(Force), Force, 'g-', label='pFJC Model')
            ax1.plot(distance_range - float(self.bead_size_var.get()), force_range, 'k.', label='Experimental Data')
            ax1.set_xlabel('Distance (um)')
            ax1.set_ylabel('Force (pN)')
            ax1.legend()
            ax1.set_title('Model Fitting')
            
            # Plot 2: Base pairs changes
            # Use Savitzky-Golay filter to smooth basepairs
            bp_filter = savgol_filter(basepairs, 31, 3)
            
            # Plot basepair changes as a function of time
            ax2.plot(time_range/1000, basepairs, color='lightgrey', linewidth=1)
            ax2.plot(time_range/1000, bp_filter, color='green', linewidth=1, markersize=1, label='Basepairs')
            ax2.set_xlabel('Time (s)')
            ax2.set_ylabel('Base pairs')
            ax2.set_title('DNA polymerase catalyzing DNA')
            
            # Adjust spacing between subplots
            self.fig.subplots_adjust(top=0.9, bottom=0.15, left=0.1, right=0.9, wspace=0.3)
            self.canvas.draw()
            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during model fitting: {str(e)}")
    
    def get_junction_position(self, distance_range, ssDNA_percentage, ssDNA_ref, dsDNA_ref, bead_size, reverse=False):
        """
        Calculate ss/dsDNA junction position based on the direction (forward or reverse)
        SSB correction factor is included in the ssDNA_percentage calculation, but we need to apply it to ssDNA_ref too
        """
        # Get SSB factor
        ssb_factor = float(self.ssb_factor_var.get())
        # Apply SSB correction to ssDNA_ref
        corrected_ssDNA_ref = ssDNA_ref - ssb_factor
        
        if reverse:
            # Reverse direction calculation
            return (distance_range - bead_size) - (ssDNA_percentage * corrected_ssDNA_ref) * (distance_range - bead_size) / ((ssDNA_percentage * corrected_ssDNA_ref) + (1 - ssDNA_percentage) * dsDNA_ref)
        else:
            # Forward direction calculation
            return (ssDNA_percentage * corrected_ssDNA_ref) * (distance_range - bead_size) / ((ssDNA_percentage * corrected_ssDNA_ref) + (1 - ssDNA_percentage) * dsDNA_ref)
    
    def track_junction(self, reverse=False, use_exo_pol=True):
        try:
            if self.time_data is None:
                messagebox.showerror("Error", "Please load a TDMS file first")
                return
            
            # Check if we have exo/pol data
            if not (self.time_from_exo_var.get() and self.time_to_exo_var.get() and 
                    self.time_from_pol_var.get() and self.time_to_pol_var.get()):
                messagebox.showerror("Error", "Please set Exo and Pol time ranges first")
                return
                
            bead_size = float(self.bead_size_var.get())
            
            # If reverse is True, we need to recalculate the junction positions with reverse=True
            if reverse:
                # Get data without junction positions first
                time_range_all, force_range_all, distance_range_all, ssDNA_all_percentage, basepairs, _ = self.get_exo_pol_analysis_data()
                if time_range_all is None:  # Error occurred in get_exo_pol_analysis_data
                    return
                    
                # Calculate reference values for dsDNA and ssDNA under exo and pol forces
                exo_force = float(self.exo_force_var.get())
                pol_force = float(self.pol_force_var.get())
                dsDNA_exo_ref = self.tWLC(exo_force)
                dsDNA_pol_ref = self.tWLC(pol_force)
                
                # Calculate ssDNA reference values
                Lss = float(self.Lss_var.get())
                b = float(self.b_var.get())
                Sss = float(self.Sss_var.get())
                
                ssDNA_exo_ref = Lss * (float(coth(exo_force * b / 4.1).evalf()) - 4.1 / (exo_force * b)) * (1 + exo_force / Sss)
                ssDNA_pol_ref = Lss * (float(coth(pol_force * b / 4.1).evalf()) - 4.1 / (pol_force * b)) * (1 + pol_force / Sss)
                
                # Get ssDNA percentages
                ssb_factor = float(self.ssb_factor_var.get())
                ssDNA_exo_percentage = (self.distance_range_exo - bead_size - dsDNA_exo_ref)/(ssDNA_exo_ref - ssb_factor - dsDNA_exo_ref)
                ssDNA_pol_percentage = (self.distance_range_pol - bead_size - dsDNA_pol_ref)/(ssDNA_pol_ref - ssb_factor - dsDNA_pol_ref)
                
                # Calculate junction positions for reverse direction
                junction_position_exo = self.get_junction_position(
                    self.distance_range_exo, 
                    ssDNA_exo_percentage, 
                    ssDNA_exo_ref, 
                    dsDNA_exo_ref, 
                    bead_size, 
                    True  # Reverse direction
                )
                
                junction_position_pol = self.get_junction_position(
                    self.distance_range_pol, 
                    ssDNA_pol_percentage, 
                    ssDNA_pol_ref, 
                    dsDNA_pol_ref, 
                    bead_size, 
                    True  # Reverse direction
                )
                
                junction_position_all = np.append(junction_position_exo, junction_position_pol)
                junction_position_all = np.array(junction_position_all, dtype=np.float64)
            else:
                # For forward direction, we can use the junction positions directly from get_exo_pol_analysis_data
                time_range_all, force_range_all, distance_range_all, ssDNA_all_percentage, basepairs, junction_position_all = self.get_exo_pol_analysis_data()
                if time_range_all is None:  # Error occurred in get_exo_pol_analysis_data
                    return
            
            # Clear previous plots
            self.fig.clear()
            
            # Create plot
            ax = self.fig.add_subplot(111)
            
            # Plot junction movement using the updated parameters
            ax.scatter(time_range_all/1000, distance_range_all - bead_size, color='black', linestyle='dashed', s=2, label='End-to-End Distance')
            ax.scatter(time_range_all/1000, junction_position_all, color='green', s=2, label='DNA Polymerase Trace')
            
            # Fill between the lines
            ax.fill_between(np.array(time_range_all/1000), distance_range_all - bead_size, junction_position_all, color='gray', alpha=0.2)
            
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Distance (um)')
            ax.set_title('ssDNA/dsDNA Junction Position' + (' (Reverse Direction)' if reverse else ''))
            
            # Invert y-axis and move x-axis to top
            ax.invert_yaxis()
            ax.xaxis.set_ticks_position('top')
            
            self.fig.tight_layout()
            self.canvas.draw()
            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while tracking junction position: {str(e)}")
    
    def plot_dna_polymerase_trace(self):
        """Plot DNA polymerase trace using exo and pol data"""
        try:
            # Reuse the track_junction method with default parameters (forward direction)
            self.track_junction(reverse=False, use_exo_pol=True)
            
            # Save the plot
            base_dir = os.path.dirname(self.filename_var.get())
            base_name = os.path.splitext(os.path.basename(self.filename_var.get()))[0]
            cycle = self.cycle_var.get()
            results_dir = os.path.join(base_dir, 'results')
            
            # Create results directory if it doesn't exist
            if not os.path.exists(results_dir):
                os.makedirs(results_dir)
                
            output_filename = os.path.join(results_dir, f"{base_name}-cycle#{cycle}-DNApTraces.eps")
            self.fig.savefig(output_filename, format='eps', dpi=300, bbox_inches='tight')
            
            # Show success message
            messagebox.showinfo("Success", f"DNA Polymerase Trace saved to {output_filename}")
            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while plotting DNA polymerase trace: {str(e)}")
    
    def apply_savgol_filter(self, data, window_length=31, polyorder=3):
        """Apply Savitzky-Golay filter to smooth data"""
        from scipy.signal import savgol_filter
        return savgol_filter(data, window_length, polyorder)
    
    def linear_segment_fitting(self):
        try:
            if self.time_data is None:
                messagebox.showerror("Error", "Please load a TDMS file first")
                return
            
            # Use exo/pol data instead of ROI data
            if (self.time_from_exo_var.get() and self.time_to_exo_var.get() and 
                self.time_from_pol_var.get() and self.time_to_pol_var.get()):
                time_range, _, _, _, basepairs, _ = self.get_exo_pol_analysis_data()
                if time_range is None:  # Error occurred in get_exo_pol_analysis_data
                    return
            else:
                messagebox.showerror("Error", "Please set Exo and Pol time ranges first")
                return
            
            # Convert time to seconds and prepare data
            time_seconds = time_range / 1000
            
            # Setup results directory
            base_dir = os.path.dirname(self.filename_var.get())
            base_name = os.path.splitext(os.path.basename(self.filename_var.get()))[0]
            cycle = self.cycle_var.get()
            results_dir = os.path.join(base_dir, 'results')
            
            # Create results directory if it doesn't exist
            if not os.path.exists(results_dir):
                os.makedirs(results_dir)
            
            # Initialize and fit the model
            segment_number = int(self.segment_number_var.get())
            pwlf = FastPWLFit(time_seconds, basepairs)
            pwlf.fit_model(segment_number)
            pwlf.plot_fit()
            
            # Save the fit results
            output_filename = os.path.join(results_dir, f"{base_name}-cycle#{cycle}_segments.csv")
            pwlf.save_results_csv(output_filename)
            
            # Save the fitted plot
            plot_filename = os.path.join(results_dir, f"{base_name}-cycle#{cycle}_segments.png")
            plt.savefig(plot_filename, format='png', dpi=300, bbox_inches='tight')
            
            # Show success message
            messagebox.showinfo("Success", f"Linear segment fitting results saved to {results_dir}")
            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during linear segment fitting: {str(e)}")
    
    def save_data(self):
        try:
            if self.time_data is None:
                messagebox.showerror("Error", "Please load a TDMS file first")
                return
            
            # Check if we have exo/pol data
            if (self.time_from_exo_var.get() and self.time_to_exo_var.get() and 
                self.time_from_pol_var.get() and self.time_to_pol_var.get()):
                # Use exo/pol data
                time_range, force_range, distance_range, ssDNA_percentage, basepairs, junction_position = self.get_exo_pol_analysis_data()
                if time_range is None:  # Error occurred in get_exo_pol_analysis_data
                    return
                
                # Apply Savitzky-Golay filter to basepairs
                basepairs_filtered = self.apply_savgol_filter(basepairs)
                
                # Extract base filename without extension and create results directory
                base_dir = os.path.dirname(self.filename_var.get())
                base_name = os.path.splitext(os.path.basename(self.filename_var.get()))[0]
                cycle = self.cycle_var.get()
                results_dir = os.path.join(base_dir, 'results')
                
                # Create results directory if it doesn't exist
                if not os.path.exists(results_dir):
                    os.makedirs(results_dir)
                
                # Define font settings for plots
                font = {'family': 'DejaVu Sans', 'weight': 'normal', 'size': 16}
                mpl.rc('font', **font)
                
                # 1. Save Excel data
                output_filename = os.path.join(results_dir, f"{base_name}-cycle#{cycle}-processedData.xlsx")
                data = {
                    'time': time_range,
                    'ssDNA_all_percentage': ssDNA_percentage,
                    'junction_position_all': junction_position,
                    'basepairs': basepairs
                }
                df = pd.DataFrame(data)
                df.to_excel(output_filename, index=False)
                
                # 2. Save Plot 1: Basepair Change (Filtered)
                plot_filename = os.path.join(results_dir, f"{base_name}-cycle#{cycle}-BasepairChange-filtered.png")
                plt.figure(figsize=(6, 4))
                plt.xlabel('Time (s)', fontdict=font)
                plt.ylabel('Basepairs', fontdict=font)
                plt.plot(time_range/1000, basepairs, color='lightgrey', linewidth=1)
                plt.plot(time_range/1000, basepairs_filtered, color='green', linewidth=1, label='Basepairs')
                plt.tight_layout()
                plt.savefig(plot_filename, format='png', dpi=300, bbox_inches='tight')
                plt.close()
                
                # 3. Save Plot 2: ssDNA Percentage (as Basepairs)
                plot_filename = os.path.join(results_dir, f"{base_name}-cycle#{cycle}-ssDNA_percentage.png")
                plt.figure(figsize=(8, 3))
                plt.ylabel('Basepairs (bp)', fontdict=font)
                plt.xlabel('Time (s)', fontdict=font)
                plt.scatter(time_range/1000, basepairs, color='black', s=0.5, label='End-to-End Distance')
                plt.tight_layout()
                plt.savefig(plot_filename, format='png', dpi=300)
                plt.close()
                
                # 4. Save Plot 3: DNA Polymerase Traces
                plot_filename = os.path.join(results_dir, f"{base_name}-cycle#{cycle}-DNApTraces.png")
                plt.figure(figsize=(6, 4))
                plt.title('Time (s)', fontdict=font)
                plt.ylabel('Distance (µm)', fontdict=font)
                plt.scatter(time_range/1000, distance_range - float(self.bead_size_var.get()), color='black', s=2, label='End-to-End Distance')
                plt.scatter(time_range/1000, junction_position, color='green', s=2, label='DNA Polymerase Trace')
                plt.fill_between(time_range/1000, distance_range - float(self.bead_size_var.get()), junction_position, color='gray', alpha=0.2)
                plt.ylim(0, 3.8)
                ax = plt.gca()
                ax.invert_yaxis()
                ax.xaxis.set_ticks_position('top')
                plt.tight_layout()
                plt.savefig(plot_filename, format='png', dpi=300)
                plt.close()
                
                # 5. Save Plot 4: Basepair Change (Raw)
                plot_filename = os.path.join(results_dir, f"{base_name}-cycle#{cycle}-BasepairChange.png")
                plt.figure(figsize=(8, 3))
                plt.xlabel('Time (s)', fontdict=font)
                plt.ylabel('Basepairs', fontdict=font)
                plt.plot(time_range/1000, basepairs, color='red', marker='o', linestyle='dashed', linewidth=2, markersize=2, label='Basepairs')
                plt.tight_layout()
                plt.savefig(plot_filename, format='png', dpi=300, bbox_inches='tight')
                plt.close()
                
                # Set matplotlib back to defaults
                mpl.rcdefaults()
                
                # 6. Save the previous plot from the GUI as an EPS file (keeping this from original code)
                output_filename_eps = os.path.join(results_dir, f"{base_name}-cycle#{cycle}-DNApTraces.eps")
                self.track_junction(reverse=False, use_exo_pol=True)
                self.fig.savefig(output_filename_eps, format='eps', dpi=300, bbox_inches='tight')
                
                messagebox.showinfo("Success", f"All data and plots saved to {results_dir}")
            else:
                messagebox.showerror("Error", "Please set Exo and Pol time ranges first")
            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while saving data: {str(e)}")

    def save_plot_if_not_exists(self, filename):
        """Check if file exists before saving to avoid overwriting"""
        if os.path.exists(filename):
            return False
        return True

if __name__ == "__main__":
    root = tk.Tk()
    app = OTDataAnalyzer(root)
    root.mainloop()
