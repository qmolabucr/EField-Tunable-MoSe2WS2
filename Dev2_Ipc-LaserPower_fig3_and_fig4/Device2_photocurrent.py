'''
Loads the raw data files and generates maps like the one in figure 3c.
The complete measurement is held in three files:
pci.npy contains the raw photocurrent data
pow.npy contains the laser power data
log.log contains measurement information and metadata.

In this measurement the heterojunction overlap was stimulated with a laser while
source/drain voltage (Vsd) and gate voltage (Vg) were applied to the device. 
Photocurrent flow through the junction was recorded. The measurement was repeated
with increasing laser power. This data is fit as described in the main text
to produce the map in fig3e and fig4.

This produces a 3d data cube with Vsd, Vg, and laser power axis
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
from scipy import ndimage
from os.path import join

path = '' #location of data files
path = 'C:/Users/jedki/Dropbox/Gabor-Lui Collaborations/Revised manuscript/Revisions for Editorial Requests/Data and Code/Dev2_Ipc-LaserPower_fig3_and_fig4'
datafile_name = 'CPS_2021_03_25_51_pci.npy'
powerfile_name = 'CPS_2021_03_25_51_pow.npy'
logfile_name = 'CPS_2021_03_25_51_log.log'

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

#The z axis, laser power, is generated from the pow.npy file
powdata = np.load(join(path, powerfile_name))
#We can simply take the average power for each Vsd-Vg map
pval = []
for i in range(powdata[0,0,:].size):
    pval.append(np.mean(powdata[:,:,i]))
pval = np.asanyarray(pval)

#load data map
data = np.load(join(path, datafile_name))
#The numpy array holding the data has shape (Vg, Vsd, Power) = (60, 115, 30)
#That is to say the measurement iterated across a 60x115 map of Vg and Vsd values
#then repeated the measurement 29 more times with increasing power.
#In this example we will just use the 16th scan (as that is the one shown in the main text)
#The lock in amplifier also introduces an offset in the data, we can correct this by universally subtracting
#the average of the first, essentially zero power, data map
data = data[:,:,16] - np.mean(data[:,:,0])

#Multiply the data by the pre-amplifier gain to calibrate to nA
data = data * float(info['Pre-Amp Gain']) * 1e9

#In this case the data is negative, but we are not concerned with the sign.
#I'm multiplying by -1 to make it easier to think about
# data = data * -1

#Generating a simple version of fig3b
colormap = 'plasma'
smooth = 1.5

#bounds on colormap scaling
zlimit = []
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