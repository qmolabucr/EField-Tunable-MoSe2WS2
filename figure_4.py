'''
Produces figure 3

The script is organized in a way to work through the figure one panel at a time
'''

import matplotlib as mpl
import matplotlib.pyplot as pl
import numpy as np
from scipy import ndimage
import scipy.optimize as op
from os.path import join, exists, dirname
from os import listdir
import figure_generator

print(' - Generating figure axes')
ax, fig = figure_generator.fig_4()


'''
----------------------------------------------------        Fig 4a -- Ipv-Vsd with Power     -----------------------
'''

#                                     ~~ Load Data ~~

path = join(dirname(__file__), 'dev2_photocurrent')

datafile_name = 'CPS_2021_03_25_51_pci.npy'
powerfile_name = 'CPS_2021_03_25_51_pow.npy'
logfile_name = 'CPS_2021_03_25_51_log.log'

#load information from log file into dictionary
with (open(join(path, logfile_name), 'r')) as f:

    info = {}
    for line in f:
        # print(line)
        s = line.split(':')
        info[s[0]] = s[1]

#Generate x and y arrays defining the axis of the data map
#In this measurement the fast axis is Vsd and the slow axis is Vg (can see this in the logfile)
xval = np.linspace(float(info["Fast Axis Start"]), float(info["Fast Axis End"]), int(info["nx"]))
yval = np.linspace(float(info["Slow Axis Start"]), float(info["Slow Axis End"]), int(info["ny"]))

#The z axis, laser power, is generated from the pow.npy file
powdata = np.load(join(path, powerfile_name))
#We can simply take the average power for each Vsd-Vg map
cval = []
for i in range(powdata[0,0,:].size):
    cval.append(np.mean(powdata[:,:,i]))
cval = np.asanyarray(cval)

#load data map
data = np.load(join(path, datafile_name))

data = data * float(info['Pre-Amp Gain'])
data = data * 1e12       # calibrates to real values
dataoffset = .0103+9.84
fieldMod = (.1205/.65) #.114 to .127
# fieldMod = (.1205/.6) #.114 to .127
xval = xval * fieldMod

#                                     ~~ Data and Panel Parameters ~~

xlabel = r'$V_{\mathrm{sd}}$ (V)'
xlabel = ''
ylabel = r'$I_{\mathrm{pc}}$  (pA)'
zlabel = r'Power ($\mu$W)'
fontsize = 7
ticksize = 7
linewidth = 3.5
markersize = 8

mapcolor_name = 'plasma' #Color of map
colormap_name = 'viridis_r' #Color of slices

#The Vg axis index used in the figure
gateindex = 14

gate = -6.5
powerslice = [0, .01, .02, .03] # Go check cval after power arranging
speccolor = []

hlines = []
vlines = [1.047]
vlcolor = 'r'
vlines = [.15, .29]
va = .85

fill = []
fill = []
gradfill = [.5, 2.4]
fillclr = ['deepskyblue', 'darkorange']
fillclr = ['red']
filla = .1

smooth = 2
savgolsmooth = 0
savgolpoly = 3

poweraxis = "z"
zeropoweroffset = True

#trims data if nonzero
xlimits = [.25, 2.25]
xlimits = [.25* fieldMod, 2.25* fieldMod]
xticks = [.05]

ylimits = [-1,160]
yticks = [20,40,60,80,100,120,140,160]


#                                     ~~ Data Manipulations ~~

xlen = xval.size
ylen = yval.size

if poweraxis == "z":
    if zeropoweroffset:
        offset_ = np.mean(data[:,:,0])
        data = data - offset_ + dataoffset
        cval = cval - np.min(cval)

if smooth > 0:
    for i in range(cval.size):
        data[:,:,i] = ndimage.gaussian_filter(data[:,:,i], smooth)

powerslice = cval[:7:1]
zticks = [0, 4, 8]

#                                      ~~ Plot Population ~~

if speccolor == []:
    cmap = pl.get_cmap(colormap_name, len(powerslice))
    clr = cmap(range(len(powerslice)))
else:
    clr = speccolor

colordata = np.asarray([np.linspace(powerslice[0]*1000, powerslice[-1]*1000, 1000), np.linspace(powerslice[0]*1000, powerslice[-1]*1000, 1000)])
colorxt=[powerslice[0]*1000, powerslice[-1]*1000, -1, 1]
ax[-1].imshow(colordata, cmap = pl.get_cmap('viridis_r'), norm = mpl.colors.Normalize(vmin = powerslice[0]*1000, vmax = powerslice[-1]*1000), extent = colorxt, aspect = 'auto', origin = 'upper')
ax[-1].set_yticks([])
ax[-1].set_xticks(zticks, zticks, color = 'k')

for i in hlines:
    ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = linewidth)
for i in range(len(vlines)):
    if i == 0:
        ax[0].axvline(vlines[i], ymin = 0, ymax = .43, linestyle = ":", c = vlcolor, linewidth = linewidth, alpha = va)
    else:
        ax[0].axvline(vlines[i], ymin = 0, ymax = .8, linestyle = ":", c = vlcolor, linewidth = linewidth, alpha = va)

x = xval
for i in range(len(powerslice)):
    powerindex = np.searchsorted(cval, powerslice[i])
    # gateindex = np.searchsorted(yval, gate)
    y = data[gateindex, :, powerindex]

    ax[0].plot(x, y, linewidth = linewidth, c = clr[i])


#Sets x and y axis parameters
if xlimits != []:
    ax[0].set_xlim(xlimits[0], xlimits[1])
if ylimits != []:
    ax[0].set_ylim(ylimits[0], ylimits[1])

ax[0].set_xticks([])
if len(yticks) > 0:
    ax[0].set_yticks(yticks)

ax[0].tick_params(axis='x', labelsize = ticksize, direction = 'in')
ax[0].tick_params(axis='y', labelsize = ticksize, direction = 'in')
ax[-1].tick_params(axis='x', labelsize = ticksize, direction = 'out', color = 'k')
ax[0].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[0].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax[-1].set_xlabel(zlabel, fontsize = fontsize, labelpad = 0, rotation = 0)


'''
----------------------------------------------------        Fig 4a inset -- Ipc - Power Scatter     -----------------------
'''

#                                            ~~ Load Data ~~

path = join(dirname(__file__), 'dev2_photocurrent')

datafile_name = 'CPS_2021_03_25_51_pci.npy'
powerfile_name = 'CPS_2021_03_25_51_pow.npy'
logfile_name = 'CPS_2021_03_25_51_log.log'

#load information from log file into dictionary
with (open(join(path, logfile_name), 'r')) as f:

    info = {}
    for line in f:
        # print(line)
        s = line.split(':')
        info[s[0]] = s[1]

#Generate x and y arrays defining the axis of the data map
#In this measurement the fast axis is Vsd and the slow axis is Vg (can see this in the logfile)
xval = np.linspace(float(info["Fast Axis Start"]), float(info["Fast Axis End"]), int(info["nx"]))
yval = np.linspace(float(info["Slow Axis Start"]), float(info["Slow Axis End"]), int(info["ny"]))

#The z axis, laser power, is generated from the pow.npy file
powdata = np.load(join(path, powerfile_name))
#We can simply take the average power for each Vsd-Vg map
cval = []
for i in range(powdata[0,0,:].size):
    cval.append(np.mean(powdata[:,:,i]))
cval = np.asanyarray(cval)

#load data map
data = np.load(join(path, datafile_name))
data = data * float(info['Pre-Amp Gain'])


#                                     ~~ Data and Panel Parameters ~~

fit = True
trimval = 150

# calibrates to pA

data = data * 1e12       
cval = cval * 1000

xlabel = r'Laser Power ($\mu$W)'
ylabel = r'$I_{\mathrm{pc}}$ (pA)'
linewidth = 2.5
markersize = 10

hlines1 = []
vlines1 = []

smooth = 1

#Vsd and Vg values to take scatter data from
gate = [-6.5]
source = [1.15]
clr = ['k']
lineclr = ['r']

poweraxis = "z"
zeropoweroffset = True
calcpower = False
zerocurve = True

xlimits1 = [-.015, .2]
xlimits1 = [0, trimval]
ylimits = [0, 120]
xticks = [0, 50, 100]
yticks = [0, 50, 100]


#                                        ~~ Data Manipulations ~~

xlen = xval.size
ylen = yval.size

if poweraxis == "z":
    if zeropoweroffset:
        offset_ = np.mean(data[:,:,0])
        data = data - offset_
        cval = cval - np.min(cval)


#                                         ~~ Plot Population ~~

def twobody(x, a, b):
    i = b * np.log(1+x*a)
    return i

for i in hlines1:
    ax[1].axhline(i, linestyle = ":", c = 'k', linewidth = linewidth)
for i in vlines1:
    ax[1].axvline(i, linestyle = ":", c = 'k', linewidth = linewidth)

x = cval
for i in range(len(gate)):
    # gateindex = np.searchsorted(yval, gate[i])
    sourceindex = np.searchsorted(xval, source[i])

    y = data[gateindex, sourceindex, :]

    if zerocurve:
        y = y - y[0]

    ax[1].scatter(x, y, s = markersize, color = clr[i], alpha = 1, zorder = 2)

    #Fits scatter
    if fit:
        try:
            print("Trying fit")
            ptrim = np.searchsorted(cval, trimval)
            par, pcov = op.curve_fit(twobody, x[:ptrim], y[:ptrim], maxfev = 3200)
            print("A = %f" % (par[0]))
            print("B = %f" % (par[1]))
            x_ = np.linspace(x[0], x[-1], 1000)
            ax[1].plot(x_, twobody(x_, *par), linewidth = linewidth, c = lineclr[i], zorder = 1)
        except RuntimeError:
            print("Fit failed")


#Sets x and y axis parameters
if xlimits1 != []:
    ax[1].set_xlim(xlimits1[0], xlimits1[1])
if ylimits != []:
    ax[1].set_ylim(ylimits[0], ylimits[1])

ax[1].set_xticks(xticks)
ax[1].set_yticks(yticks)
ax[1].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'k')
ax[1].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'k')
ax[1].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[1].set_ylabel(ylabel, fontsize = fontsize, labelpad = 12, rotation = 270)
ax[1].yaxis.tick_right()
ax[1].yaxis.set_label_position('right')


'''
----------------------------------------------------        Fig 4b -- gamma/(alpha + beta) vs Vsd     -----------------------
'''

#                                            ~~ Load Data ~~

path = join(dirname(__file__), 'dev2_photocurrent')
name = '2021_03_25_51_powercube_twobody3_'

savefile = join(path, name + ".npz")
if exists(savefile):
    f = np.load(savefile, allow_pickle=True)

    alpha = f['alpha']
    beta = f['beta']
    alphaerror = f['alphaerror']
    betaerror = f['betaerror']
    rsq = f['rsq']
    xval1 = f['xval']
    yval1 = f['yval']

#Invalidly large number filter
xval = xval1 * fieldMod
cutoff = 1e7
if cutoff > 0:
    for i in range(xval1.size):
        for p in range(yval1.size):
            if alpha[p,i] > cutoff or alpha[p,i] < 0:
                alpha[p,i] = 0

#Upper error filter
rsq_trim = 0.7
if rsq_trim > 0:
    for i in range(xval1.size):
        for p in range(yval1.size):
            if np.abs(rsq[p,i]) < rsq_trim:
                alpha[p,i] = 0

#hotpixel filter
hotpixels = True
hp_strength = 2
if hotpixels:
    for p in range(yval1.size):
        for i in range(xval1.size - 2):
            localaverage = (alpha[p,i] + alpha[p,i+2])/2
            if alpha[p,i+1] > localaverage * hp_strength:
                alpha[p,i+1] = localaverage


#                                     ~~ Data and Panel Parameters ~~

xlabel = r"$\Delta$E (V/nm)"
xlabel = "Electric Field (V/nm)"
ylabel = r'$\gamma$ / $(\alpha + \beta)$ (arb. units)'
zlabel = "Laser Power (mW)"
linewidth = 3.5
markersize = 8

colormap_name = 'viridis_r' #Color of slices

powerslice = [0, .01, .02, .03] # Go check cval after power arranging
speccolor = []

hlines = []

smooth = .5

#Axis parameters
ylog = True
ylimits = [-5.5, 12]
ylimits = [2, 2e6]
xticks = [.5, 1, 1.5, 2]
xticks = [.1, .2, .3, .4]
yticks = [ 1, 10, 100, 1000, 10000]


#                                        ~~ Data Manipulations ~~

xlen = xval.size
ylen = yval.size

if smooth > 0:
    alpha = ndimage.gaussian_filter(alpha, smooth)

powerslice = cval[:8:1]

#                                         ~~ Plot Population ~~

for i in hlines:
    ax[2].axhline(i, linestyle = "dotted", c = 'k', linewidth = linewidth)
for i in vlines:
    # ax[2].axvline(i, ymin = .26, linestyle = "dotted", c = vlcolor, linewidth = linewidth, alpha = va)
    ax[2].axvline(i, ymin = 0, linestyle = "dotted", c = vlcolor, linewidth = linewidth, alpha = va)

# gateindex = np.searchsorted(yval, gate)
y = alpha[gateindex,:]
ax[2].plot(xval, y, linewidth = linewidth, c = 'k')

#Sets x and y range and axis labels
if xlimits != []:
    ax[2].set_xlim(xlimits[0], xlimits[1])
if ylimits != []:
    ax[2].set_ylim(ylimits[0], ylimits[1])

ax[2].set_xticks(xticks)
ax[0].xaxis.tick_top()
ax[0].xaxis.set_label_position('top')
xticktop = [np.around(xtick/fieldMod, 1) for xtick in xticks]
ax[0].set_xticks(xticks, xticktop)
ax[0].set_xlabel(r'$V_{\mathrm{sd}}$  (V)', fontsize = fontsize, labelpad = 5)
ax[2].set_yticks(yticks)

if ylog:
    ax[2].set_yscale('log')

ax[2].tick_params(axis='x', labelsize = ticksize, direction = 'in')
ax[2].tick_params(axis='y', labelsize = ticksize, direction = 'in')
ax[2].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[2].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)


pl.show()