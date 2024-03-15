'''
Loads the raw data files and generates a simple version of the map in figure 3b.
The file ending with pci.dat contains the raw data while the file ending with log.log contains
measurement information and metadata.

In this measurement source/drain voltage (Vsd) and gate voltage (Vg) were applied to the device 
while current flow through the junction was recorded.

This produces a 2d data map with Vsd and Vg axis
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
from scipy import ndimage
from os.path import join

path = '' #location of data files
path = 'C:/Users/jedki/Dropbox/Gabor-Lui Collaborations/Revised manuscript/Revisions for Editorial Requests/Data and Code/Dev2_dark_transport_fig3b'
datafile_name = 'CPS_2021_03_18_1_pci.dat'
logfile_name = 'CPS_2021_03_18_1_log.log'

#load information from log file into dictionary
with (open(join(path, logfile_name), 'r')) as f:

    info = {}
    for line in f:
        print(line)
        s = line.split(':')
        info[s[0]] = s[1]

#Generate x and y arrays defining the axis of the data map
#In this measurement the fast axis is Vsd and the slow axis is Vg (can see this in the logfile)
xval = np.linspace(float(info["Fast Axis Start"]), float(info["Fast Axis End"]), int(info["nx"]))
yval = np.linspace(float(info["Slow Axis Start"]), float(info["Slow Axis End"]), int(info["ny"]))

#load data map
data = np.loadtxt(join(path, datafile_name))

#Multiply by the pre-amplifier gain to calibrate to nA
data = data * float(info['Pre-Amp Gain']) * 1e9

#In this case the data is negative, but we are not concerned with the sign.
#I'm multiplying by -1 to make it easier to think about
data = data * -1

#Generating a simple version of fig3b
colormap = 'bone'
smooth = 1

#bounds on colormap scaling
zlimit = [0, 8]
zlimit = [] #if no bounds are provided sets it to the maximum range of the dataset

#apply some gaussian smoothing to deal with high frequency noise
if smooth > 1:
    data = ndimage.gaussian_filter(data, smooth)


#Get colormap and color normalization
cmap = pl.get_cmap(colormap)
if len(zlimit) == 0:
    cnorm = mpl.colors.Normalize(np.min(data), np.max(data))
else:
    cnorm = mpl.colors.Normalize(zlimit[0], zlimit[1])
#define the edges of the map
xt = [xval[0], xval[-1], yval[0], yval[-1]]

#generate figure and plot
pl.figure()
pl.imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
pl.xlabel(r'$V_{sd}$ (V)')
pl.ylabel(r'$V_{g}$ (V)')
pl.show()