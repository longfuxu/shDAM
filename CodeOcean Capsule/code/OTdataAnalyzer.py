# Import libraries
from __future__ import division
import os
from sympy import *
from sympy import coth
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from nptdms import TdmsFile
# from more_itertools import chunked

# Set project root for robust file paths
# Assumes this script is in a subdirectory of the project root (e.g., 'python_scripts')
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)

# Set default file path for the TDMS file relative to the project root
filename_rel = 'data/force data/OTdata_example_30nM DNAp + trx  + 625uM dNTPs.tdms'
filename = os.path.join(project_root, filename_rel)

# Create results directory if it doesn't exist
results_dir = os.path.join(project_root, 'results')
os.makedirs(results_dir, exist_ok=True)

# Read data from the TDMS file
tdms_file = TdmsFile(filename)
# Extract and convert time, force, and distance data to numpy arrays
time = np.array(tdms_file['FD Data']['Time (ms)'][:])
force = np.array(tdms_file['FD Data']['Force Channel 0 (pN)'][:])
distance = np.array(tdms_file['FD Data']['Distance 1 (um)'][:])

# Plotting force against time
plt.figure(figsize=(6,4))
plt.plot(time, force)
plt.ylabel('force/pN')
plt.tight_layout()
plt.savefig(os.path.join(results_dir, 'force_vs_time.png'), dpi=300, bbox_inches='tight')
plt.close()

# Set default values for cycle and time ranges
cycle = '01'
time_from_exo = 10500  # starting time of exo in ms
time_to_exo = 62140  # ending time of exo in ms
time_from_pol = 62713  # starting time of pol in ms
time_to_pol = 121910  # ending time of pol in ms

# here you will analyze your interesting cycle with corresponding time range; 
# define a temporary index to compute time of ROI, and subsuquent ROI of force, distance
# indtemp_all = np.where((time <= time_to_all) & (time >= time_from_all))
indtemp_exo = np.where((time <= time_to_exo) & (time >= time_from_exo))
indtemp_pol = np.where((time <= time_to_pol) & (time >= time_from_pol))

# exo time range of ROI
time_range_exo = time[indtemp_exo]
force_range_exo = force[indtemp_exo]
distance_range_exo = distance[indtemp_exo]

# pol time range of ROI
time_range_pol = time[indtemp_pol]
force_range_pol = force[indtemp_pol]
distance_range_pol = distance[indtemp_pol]

# all time range of ROI
time_range_all = np.append(time_range_exo, time_range_pol)
force_range_all = np.append(force_range_exo, force_range_pol)
distance_range_all = np.append(distance_range_exo, distance_range_pol)

# parameters for tWLC model: Peter Gross, et al. Nature Physics volume 7, pages731–736(2011)
# dsDNA contour length Lc = 2.85056um; persistent length Lp = 56nm
# the twist rigidity C=440 pN nm2;
# the stretching modulus S=1500 pN;
# the twist–stretch coupling g(F) is given by: g(F) =g0+g1F,where g0=−637 pN nm, g1=17 nm
EEDds,Lc,F,Lp,C,g0,g1,S = symbols('EEDds Lc F Lp C g0 g1 S', real=True)
C = 440
g0= -637
g1 = 17
Lc = 2.85056
Lp = 56
S = 1500
# tWLC model expression:
def tWLC(F):
    EEDds = Lc*(1-0.5*(4.1/(F*Lp))**0.5 + C*F/(-(g0+g1*F)**2 + S*C))
    return (EEDds)

# parameters for FJC model: Smith, S. B., et al. Science 271, 795–799 (1996).
# ssDNA contour length Lss = 4.69504um,
# Kuhn length b = 1.5nm (persistent length is 0.75nm),
# the stretching modulus S=800pN
EEDss,Lss,b,Sss = symbols('EEDss Lss b Sss', real=True)
Lss = 4.69504
b = 1.5
Sss = 800
# FJC model expression:
def FJC(F):
    EEDss = []
    for Fext in F:
        x  = Lss * (coth(Fext * b / 4.1) - 4.1 / (Fext * b)) * (1 + Fext / Sss)
        EEDss.append(x)
    EEDss = np.array(EEDss)
    return (EEDss)

Force = np.linspace(0.1,68,1000)

plt.figure(figsize=(6,4))
font = {'family': 'Arial', 'weight': 'normal', 'size': 16}
plt.xlabel('Distance (um)',fontdict=font)
plt.ylabel('Force(pN)',fontdict=font)
plt.plot(tWLC(Force),Force,color='red', marker='o', linestyle='dashed',linewidth=2, markersize=2,label='tWLC Model')
plt.plot(FJC(Force),Force,color='b', marker='o', linestyle='dashed',linewidth=2, markersize=2,label='FJC Model')

# Default bead size = 1.76um, by running this code to fit the right bead size
bead_size = 1.78

# plot experimental data together with 
plt.plot(distance_range_all - bead_size, force_range_all, marker='o', linestyle='solid',linewidth=2, markersize=2,label='Experimental Data')
plt.legend()
plt.ylim(-5,72)
plt.tight_layout()
plt.savefig(os.path.join(results_dir, 'force_distance_models.png'), dpi=300, bbox_inches='tight')
plt.close()

# Default exo_force is 50pN and pol_force 20pN 
exo_force = 50
pol_force = 20

# calculating the length of dsDNA under exo and pol force
dsDNA_exo_ref = tWLC(exo_force)
dsDNA_pol_ref = tWLC(pol_force)

# calculating the length of ssDNA under exo and pol force
ssDNA_exo_ref = 4.69504 * (coth(exo_force * 1.5 / 4.1) - 4.1 / (exo_force * b)) * (1 + exo_force / 800)
ssDNA_pol_ref =4.69504 * (coth(pol_force * 1.5 / 4.1) - 4.1 / (pol_force * b)) * (1 + pol_force / 800)

# calculating ssDNA under exo and pol force
ssDNA_exo_percentage = (distance_range_exo - bead_size - dsDNA_exo_ref)/(ssDNA_exo_ref - dsDNA_exo_ref)
ssDNA_pol_percentage = (distance_range_pol - bead_size - dsDNA_pol_ref)/(ssDNA_pol_ref - dsDNA_pol_ref)
ssDNA_all_percentage = np.append(ssDNA_exo_percentage,ssDNA_pol_percentage)

# calculating basepairs
# Construct of pkyb1 DNA has a length of 8393bp
basepairs = (1-ssDNA_all_percentage) * 8393

# plot ssDNA% as a function of time
plt.figure(figsize=(8,3))
font = {'family': 'Arial', 'weight': 'normal', 'size': 16}
plt.ylabel('Basepairs',fontdict=font)
plt.xlabel('Time(s)',fontdict=font)
plt.scatter(time_range_all/1000,ssDNA_all_percentage  * 8393,color='black', linestyle='dashed',s=0.5)
plt.tight_layout()
plt.savefig(os.path.join(results_dir, 'basepairs_vs_time.png'), dpi=300, bbox_inches='tight')
plt.close()

# calculating ssDNA/dsDNA junction position under exo and pol force
junction_position_exo = (ssDNA_exo_percentage * ssDNA_exo_ref) * (distance_range_exo - bead_size) / ((ssDNA_exo_percentage * ssDNA_exo_ref) + (1 - ssDNA_exo_percentage) * dsDNA_exo_ref)
junction_position_pol = (ssDNA_pol_percentage * ssDNA_pol_ref) * (distance_range_pol - bead_size) / ((ssDNA_pol_percentage * ssDNA_pol_ref) + (1 - ssDNA_pol_percentage) * dsDNA_pol_ref)
junction_position_all = np.append(junction_position_exo, junction_position_pol)

# ## Under some occasions, DNA polymerase starts from the other side, then we need to ajust the above math expression 
# # You need to decide based on the image data
# junction_position_exo = (distance_range_exo - bead_size) - (ssDNA_exo_percentage * ssDNA_exo_ref) * (distance_range_exo - bead_size) / ((ssDNA_exo_percentage * ssDNA_exo_ref) + (1 - ssDNA_exo_percentage) * dsDNA_exo_ref)
# junction_position_pol = (distance_range_pol - bead_size) - (ssDNA_pol_percentage * ssDNA_pol_ref) * (distance_range_pol - bead_size) / ((ssDNA_pol_percentage * ssDNA_pol_ref) + (1 - ssDNA_pol_percentage) * dsDNA_pol_ref)
# junction_position_all = np.append(junction_position_exo, junction_position_pol)
# plot DNA polymerase trace as a function of time
plt.figure(figsize=(6,4))
font = {'family': 'Arial', 'weight': 'normal', 'size': 16}

plt.title('Time (s)',fontdict=font)
plt.ylabel('Distance(um)',fontdict=font)
plt.scatter(time_range_all/1000,distance_range_all - bead_size,color='black', linestyle='dashed',s=2,label='End-to-End Distance')
plt.scatter(time_range_all/1000,junction_position_all,color='green',s=2,label='DNA Polymerase Trace')
junction_position_all = np.array(junction_position_all, dtype=np.float64)
# plt.fill_between(np.array(time_range_all/1000), np.array(junction_position_all),alpha=0.5)
# Fill the area between the two data lines
plt.fill_between(np.array(time_range_all/1000), distance_range_all - bead_size, junction_position_all,color='gray',alpha=0.2)

# plt.legend()
# plt.ylim(0,3.8)
# plt.xlim(0,159)

ax = plt.gca()
ax.invert_yaxis()
ax.xaxis.set_ticks_position('top')


plt.tight_layout()
plt.savefig(os.path.join(results_dir, 'DNAp_traces.png'), dpi=300, bbox_inches='tight')
plt.close()

# plt.savefig(filename.replace('.tdms', '-cycle#') + cycle + '-DNApTraces'+'.eps', dpi=300)
# plt.savefig(filename.replace('.tdms', '-cycle#') + cycle + '-DNApTraces'+'.png', dpi=300)

# plot basepair changes as a function of time
plt.figure(figsize=(6,4))
font = {'family': 'Arial', 'weight': 'normal', 'size': 16}

plt.xlabel('Time (s)',fontdict=font)
plt.ylabel('Basepairs',fontdict=font)
plt.plot(time_range_all/1000,basepairs,color='red', marker='o', linestyle='dashed',linewidth=2, markersize=2,label='Basepairs')

# plt.legend()
# plt.title('BasePair Changes as a Function of Time',fontdict=font)
# plt.xlim(0,158)
# plt.ylim(6500,8500)
plt.tight_layout()
plt.savefig(os.path.join(results_dir, 'basepair_changes.png'), dpi=300, bbox_inches='tight')
plt.close()

# plt.savefig(filename.replace('.tdms', '-cycle#') + cycle + '-BasepairChange'+'.eps', format='eps', dpi=300, bbox_inches='tight')
# plt.savefig(filename.replace('.tdms', '-cycle#') + cycle + '-BasepairChange'+'.png', format='png', dpi=300, bbox_inches='tight')

# plt.close()

# Save all the analyzed data in an excel file
base_name = os.path.splitext(os.path.basename(filename))[0]
excel_filename = os.path.join(results_dir, f"{base_name}-cycle#{cycle}-processedData.xlsx")

data = {'time': time_range_all,
        'ssDNA_all_percentage': ssDNA_all_percentage,
        'junction_position_all': junction_position_all,
        'basepairs': basepairs}
df = pd.DataFrame(data)

with pd.ExcelWriter(excel_filename) as writer:
    df.to_excel(writer, index=False)

print(f"Processed data saved to: {excel_filename}")

