'''
Produces figure 2

The script is organized in a way to work through the figure one panel at a time
'''

import matplotlib as mpl
import matplotlib.pyplot as pl
import numpy as np
from scipy import ndimage
from scipy import interpolate as interp
from os.path import join, dirname
from os import listdir
import figure_generator


print(' - Generating figure axes')
ax, fig = figure_generator.fig_2()

# Figure display variables
fontsize = 7
ticksize = 7
insetfontsize = 5
mapcolor_name = 'terrain_spec' #Colorscale for map


'''
----------------------------------------------------        Fig 2b -- Charge dependant PL     -----------------------
'''

path = '' #location of data files
path = join(dirname(__file__), 'dev1_charge-density-dependent_PL', 'data')

#Empty lists to hold axis information and data
gate = []   #Vbg=Vtg voltages
wval = []   #wavelength
data = []

#pull every filename in the data directory
f = listdir(path)
print(' - Loading Charge Dependent PL Data')
for file in f:
    #if no wavelength data exists, generate it from the first file loaded. Every file has the same wavelength axis.
    if data == []:
        d = np.loadtxt(join(path, file), skiprows=1)
        wval = d[:,0]

    #get the gate voltage from the file name and append it to the votlage axis
    name = file.split("_")
    v1 = name[1].split("=")[1]
    v1 = float(v1[:len(v1)-5])

    gate.append(v1)

    #load the data and append it to the data map
    d = np.loadtxt(join(path, file), skiprows=1)
    data.append(d[:,1])


#the datamap and gate voltage axis and now filled, we just need to sort them
#sorts data
xy = sorted(zip(gate, data))
data = [x for y, x in xy]
data = np.asarray(data)
#sorts gate voltage, labelled as yval as it will be the y-axis
yval = np.asarray(gate)
yval = np.sort(yval)

#converts wavelength axis to energy, labelled as xval as it will be the x-axis
xval = 1240/wval
#Interpolates datamap to compensate for axis inversion
f = interp.interp2d(xval, yval, data, kind = 'linear')
ev = np.linspace(xval[-1], xval[0], wval.size)
data = f(ev, yval)

#Axis labels
xlabel = "Photon Energy (eV)"
zlabel = "PL Intensity (arb. units)"
ylabel = r'$V_{\mathrm{bg}} = V_{\mathrm{tg}}$ (V)'

#Guassian smoothing of data
smooth = 0

#trims data if nonzero
xlimits = []
ylimits = []
#tick labels
xticks = [1.55, 1.60, 1.65, 1.70]
yticks = [-6, -4, -2, 0, 2, 4, 6]

#colorscale and colorbar labels
zlimit = [0,2]
zticks = [0, 1,  2]

#take log of data?
zlog = False

#                                     ~~ Data Manipulations ~~

xlen = xval.size
ylen = yval.size

# PL scale is arbitrary so normalize to 2 and offset slightly to improve colormap
data = data - np.min(data)
normval = np.max(data)
data = 2 * (data / np.max(data)) -.01

if smooth > 1:
    data = ndimage.gaussian_filter(data, smooth)

if zlog:
    data = np.sign(data) * np.log(np.abs(data)+1)


#                                      ~~ Plot Population ~~


#Get colormap and color normalization
cmap = pl.get_cmap(mapcolor_name)
if len(zlimit) == 0:
    cnorm = mpl.colors.Normalize(np.min(data), np.max(data))
else:
    cnorm = mpl.colors.Normalize(zlimit[0], zlimit[1])

#define the edges of the map
xt = [xval[-1], xval[0], yval[0], yval[-1]]

#Data map
ax[0].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')

#Scalebar
colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[-2].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = 'lower')
ax[-2].set_xticks([])
ax[-2].set_yticks(zticks, zticks, color = 'w')


#Sets map edges and labelling
if xlimits == []:
    None
else:
    ax[0].set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[0].set_ylim(ylimits[0], ylimits[1])

if len(xticks) > 1:
    ax[0].set_xticks(xticks)
if len(yticks) > 1:
    ax[0].set_yticks(yticks)

ax[0].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[0].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax[-2].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = 'w', labelcolor = 'w')
ax[-2].spines['top'].set_color('white')
ax[-2].spines['bottom'].set_color('white')
ax[-2].spines['left'].set_color('white')
ax[-2].spines['right'].set_color('white')
ax[0].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[0].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax[-2].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 8, color = 'w', rotation = 270)


'''
----------------------------------------------------        Fig 2c -- Field dependant PL     -----------------------
'''


path = join(dirname(__file__), 'dev1_electric-field-dependent_PL', 'data')

#Empty lists to hold axis information and data
gate = []   #Vbg=Vtg voltages
wval = []   #wavelength
data = []

#pull every filename in the data directory
f = listdir(path)
print(' - Loading Charge Dependent PL Data')
for file in f:
    #if no wavelength data exists, generate it from the first file loaded. Every file has the same wavelength axis.
    if data == []:
        d = np.loadtxt(join(path, file), skiprows=1)
        wval = d[:,0]

    #get the gate voltage from the file name and append it to the votlage axis
    name = file.split("_")
    v1 = name[1].split("=")[1]
    v1 = float(v1[:len(v1)-5])

    gate.append(v1)

    #load the data and append it to the data map
    d = np.loadtxt(join(path, file), skiprows=1)
    data.append(d[:,1])

#the datamap and gate voltage axis and now filled, we just need to sort them
#sorts data
xy = sorted(zip(gate, data))
data = [x for y, x in xy]
data = np.asarray(data)
#sorts gate voltage, labelled as yval as it will be the y-axis
yval = np.asarray(gate)
yval = np.sort(yval)

#converts wavelength axis to energy, labelled as xval as it will be the x-axis
xval = 1240/wval
#Interpolates datamap to compensate for axis inversion
f = interp.interp2d(xval, yval, data, kind = 'linear')
ev = np.linspace(xval[-1], xval[0], wval.size)
data = f(ev, yval)


data = data - np.min(data)
normval = np.max(data)
data = 6 * (data / np.max(data))
data = data + np.min(data) - 0.01

xlabel = "Photon Energy (eV)"
zlabel = "PL Intensity (arb. units)"
ylabel = r'$V_{\mathrm{bg}} = -V_{\mathrm{tg}}$ (V)'
ylabel1 = "Electric Field (V/nm)"

insetfontsize = 7

mapcolor_name = 'terrain_spec' #Color of map

smooth = 0

#trims data if nonzero
xlimits = []
ylimits = []
xticks = [1.55, 1.60, 1.65, 1.70]
yticks = [-8, -6, -4, -2, 0, 2, 4, 6, 8]

zlimit = [0,6]
zticks = [0, 2, 4, 6]
zlimit1 = [yval[0]/20.01, yval[-1]/20.01]
zticks_1 = [-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4]

zlog = False

#Field axis labels
ax0 = ax[1].twinx()
ax0.set_ylim(zlimit1)
ax0.set_yticks(zticks_1)

ax1 = ax[1].twiny()
ax1.set_xlim(np.min(xval), np.max(xval))
ax1.set_xticks(xticks, ['' for ll in xticks])


#                                     ~~ Data Manipulations ~~

if smooth > 1:
    data = ndimage.gaussian_filter(data, smooth)

if zlog:
    data = np.sign(data) * np.log(np.abs(data)+1)


#                                      ~~ Plot Population ~~


#Get colormap and color normalization
cmap = pl.get_cmap(mapcolor_name)
if len(zlimit) == 0:
    cnorm = mpl.colors.Normalize(np.min(data), np.max(data))
else:
    cnorm = mpl.colors.Normalize(zlimit[0], zlimit[1])
#define the edges of the map
xt = [xval[-1], xval[0], yval[0], yval[-1]]

ax[1].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[-1].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = 'lower')
ax[-1].set_xticks([])
ax[-1].set_yticks(zticks, zticks, color = 'w')

#Sets edges of image
if xlimits == []:
    None
else:
    ax[1].set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[1].set_ylim(ylimits[0], ylimits[1])

if len(xticks) > 1:
    ax[1].set_xticks(xticks)
if len(yticks) > 1:
    ax[1].set_yticks(yticks)

ax0.tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax1.tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[1].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[1].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax[-1].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = 'w', labelcolor = 'w')
ax[-1].spines['top'].set_color('white')
ax[-1].spines['bottom'].set_color('white')
ax[-1].spines['left'].set_color('white')
ax[-1].spines['right'].set_color('white')
ax[1].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[1].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax0.set_ylabel(ylabel1, fontsize = fontsize, labelpad = 3, rotation = 90)
ax[-1].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 10, color = 'w', rotation = 270)


#                                      ~~ Show Figure ~~


pl.show()