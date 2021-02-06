"""Script for plotting Joint Probability Density Function (JPDF)
Source: Quadrant Analysis in Turbulence Research:History and Evolution By James M. Wallace
JPDF of u_z and Temp_fluctuation.
u_z is the vertical velocity
T_prime = T- T_mean
"""

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import math

#%% USER DEFINED FUNCTIONS
def cart2pol(x, y, z):
    """Conversion from Cartesian to Cylinderical
    """
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    z = z
    return(rho, phi, z)

#%% USER INPUT

total_points = 800    #total probes from the history file recorded
probe_per_line = 40   #probe per line
probe_1 = 0           #python starts from 0 so prob start from 0 tp 799, give the first probe number to +40
num_probe = 8

probes_point = np.array([0, 80])
probes_point = probes_point.astype(int)
filename = "/usr/scratch4/DATA-TEST/cyl_Ra07_2.his"
 
#%% MAIN PROGRAM
data_cordi = pd.read_table(filename, delim_whitespace=True,  skiprows=1, 
                           nrows=total_points, header=None)
data_cordi = np.array(data_cordi)
first_point_probe = 0

x_n = data_cordi[first_point_probe:first_point_probe+probe_per_line,0]
y_n = data_cordi[first_point_probe:first_point_probe+probe_per_line,1]
z_n = data_cordi[first_point_probe:first_point_probe+probe_per_line,2]

# Converting all 800 points from Cartesian to Polar
# i.e. x,y,z -> r, theta, z
data_cordi_cyl = np.array(cart2pol(data_cordi[:,0], data_cordi[:,1], 
                                   data_cordi[:,2]) ).T

# Following is just to convert our array to pandas DF
data_cordi_cyl_pan = pd.DataFrame(data_cordi_cyl)
data_cordi_cyl_pan = data_cordi_cyl_pan.rename(columns={0: "r", 1: "theta", 2: "z"})


z_H = data_cordi_cyl_pan["z"].unique()
#So 5 plane' z = 0.01049561, 0.25113645, 0.5       , 0.74886358, 0.9895044

# Check the points on circumeference
# Here, I need to use range of "r" due to some machine lever numerical difference
theta_H = data_cordi_cyl_pan["z"].unique()
theta_cir = data_cordi_cyl_pan.loc[(data_cordi_cyl_pan["z"] == z_H[0]) &
                       (data_cordi_cyl_pan['r'] >= .24) & 
                       (data_cordi_cyl_pan['r'] <= .26)]

#%% Read from history point
# Read the data from file by skipping first 801 points

"""
This history point contains 6 line of data
 0 , 1, 2, 3, 4, 5 
time, U, V, W, Pres,Temp, 
first rows has the sensor coordinates
"""
data = pd.DataFrame(np.loadtxt(filename, skiprows = total_points+1, 
                               usecols = (0,1,2,3,4,5)))
time_steps = int(data.shape[0]/total_points)

data = data.rename(columns={0: "time", 1: "u", 2: "v", 3: "w", 4: "p", 5: "T"})
# I am adding the cylinderical coordinate to the data

data = pd.concat([data, pd.DataFrame(np.tile(data_cordi_cyl, (time_steps, 1)))], 
                   axis=1)

data = data.rename(columns={0: "r", 1: "theta", 2: "z"})

# Extract time for plotting latter
t = data.time.unique()    



# Here one can change the plane by changing z_H[0] to desired plane
H_0 = np.array(data.loc[(data["z"] == z_H[0]) ]) #<-----

#H_0 = H_0[H_0[:,7].argsort()] 
#H_0 = H_0[H_0[:,0].argsort(kind='mergesort')]



# Filter out tempR

which_var = 5

T_plane_pan0 = pd.DataFrame(H_0[:,which_var])  #<-----
#T_plane_pan0 = pd.DataFrame(T_plane_pan0.values.reshape([time_steps,num_probe])).copy().transpose()
## Sorting as per tempR because, higher tempR would be on one side
## Doing this because theta is not as per tempR value
#T_plane_pan0 = T_plane_pan0.sort_values(1, ascending= False)
#T_plane_pan = T_plane_pan0.copy()



# Filter out velocity-z
which_var = 3 

Uz_plane_pan0 = pd.DataFrame(H_0[:,which_var])  #<-----
#Uz_plane_pan0 = pd.DataFrame(Uz_plane_pan0.values.reshape([time_steps,num_probe])).copy().transpose()
#
## HERE I AM MATCHING THE INDEX SAME AS TEMPR, instead of sorting
#Uz_plane_pan0 = Uz_plane_pan0.reindex(T_plane_pan0.index)
#
## Sorting as per tempR because, higher tempR would be on one side
## Doing this because theta is not as per tempR value
#Uz_plane_pan = Uz_plane_pan0.copy()


#%% Joint probability density function

# People used Uz and T_prime
uz4JPDF = Uz_plane_pan0.values.flatten()
T4JPDF = T_plane_pan0.values.flatten() - np.mean(T_plane_pan0.values.flatten())

#out_mid_plane = np.array([uz4JPDF, T4JPDF]).T
#np.savetxt("out_mid_plane.csv", out_mid_plane, delimiter = ",")


# Basically, we create 2 columns and perform this analysis with 2D Histogram
# So one can replace with any variable after the any kind of filtering
data4JPDF= [uz4JPDF, T4JPDF]
data4JPDF= np.transpose(data4JPDF)


numBins = 50  # number of bins in each dimension

sandR=np.max(np.abs(data4JPDF[:,0])) 
sandQ=np.max(np.abs(data4JPDF[:,1]))

xi = np.linspace(-sandR, sandR, 101)
yi = np.linspace(-sandQ, sandQ, 101)

Rx, Qy = np.meshgrid(xi, yi)

hst, edges = np.histogramdd(data4JPDF, bins = [xi, yi])

dx = xi[1]-xi[0]
dy = yi[1]-yi[0]
area = dx*dy;

pdfData = hst/sum(sum(hst))/area;

xvv, yvv = np.meshgrid(np.linspace(-sandR, sandR, 100), 
                       np.linspace(-sandQ, sandQ, 100))

plt.contourf(xvv, yvv, np.log10(pdfData))
plt.xlabel("Uz")
plt.ylabel("TPrime")
plt.xlim([-.4, .4])
plt.ylim([-.4, .4])
plt.colorbar()
plt.show()
