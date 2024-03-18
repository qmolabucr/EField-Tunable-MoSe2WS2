'''
Produces figure 3

The script is organized in a way to work through the figure one panel at a time
'''

import matplotlib as mpl
import matplotlib.pyplot as pl
import numpy as np
from scipy import ndimage
from os.path import join, exists, dirname
from os import listdir
import figure_generator


print(' - Generating figure axes')
ax, fig = figure_generator.fig_3()


'''
----------------------------------------------------        Fig 3b -- Dark Transport     -----------------------
'''

plotno = 1
barno = 2

#                                     ~~ Load Data ~~

path = '' #location of data files
path = join(dirname(__file__), 'dev2_dark_transport')

datafile_name = 'CPS_2021_03_18_1_pci.dat'
logfile_name = 'CPS_2021_03_18_1_log.log'

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

#load data map
data = np.loadtxt(join(path, datafile_name))

#Multiply by the pre-amplifier gain to calibrate to nA
data = data * float(info['Pre-Amp Gain']) * 1e9

#In this case the data is negative, but we are not concerned with the sign.
#I'm multiplying by -1 to make it easier to think about
data = data * -1

#                                     ~~ Data and Panel Parameters ~~

#Axis labels
xlabel = r'$V_{\mathrm{sd}}$ (V)'
zlabel = "${I}$ (nA)"
ylabel = r'$V_{\mathrm{g}}$ (V)'

#Font size
fontsize = 7
ticksize = 7
insetfontsize = 5
#Line curve
linewidth = 3
#Axis color of inset
insetcolor = 'w'
#Colormap
colormap = "bone"
#Gaussian smoothing parameter, 0 is off
smooth = 1

#map range
xlimits = [-2, 2.5]
ylimits = []
#colorscale and colorbar ticks
zlimit = [0, 8]
zticks = [0,4,8]


#                                     ~~ Data Manipulations ~~

xlen = xval.size
ylen = yval.size

if smooth > 1:
    data = ndimage.gaussian_filter(data, smooth)


#                                      ~~ Plot Population ~~

#Get colormap and color normalization
cmap = pl.get_cmap(colormap)
if len(zlimit) == 0:
    cnorm = mpl.colors.Normalize(np.min(data), np.max(data))
else:
    cnorm = mpl.colors.Normalize(zlimit[0], zlimit[1])
#define the edges of the map
xt = [xval[0], xval[-1], yval[0], yval[-1]]

#data map
ax[plotno].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')

#color bar
colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[barno].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = 'lower')
ax[barno].set_xticks([])
ax[barno].yaxis.tick_right()
ax[barno].yaxis.set_label_position('right')
ax[barno].set_yticks(zticks,zticks, color = 'w')

#Sets edges of image and labels
if xlimits == []:
    None
else:
    ax[plotno].set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[plotno].set_ylim(ylimits[0], ylimits[1])

ax[plotno].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[plotno].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax[barno].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = insetcolor, labelcolor = insetcolor)
ax[barno].yaxis.label.set_color('white')
ax[barno].spines['top'].set_color('white')
ax[barno].spines['bottom'].set_color('white')
ax[barno].spines['left'].set_color('white')
ax[barno].spines['right'].set_color('white')
ax[plotno].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[plotno].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax[barno].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = insetcolor)


'''
----------------------------------------------------        Fig 3c -- Photocurrent Map     -----------------------
'''

plotno = 4
barno = 5

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
pval = []
for i in range(powdata[0,0,:].size):
    pval.append(np.mean(powdata[:,:,i]))
pval = np.asanyarray(pval)

#load data map
data = np.load(join(path, datafile_name))
data = data * float(info['Pre-Amp Gain']) * 1e9
data = data[:,:,-14] - np.mean(data[:,:,0])


#                                     ~~ Data and Panel Parameters ~~


xlabel = r'$V_{\mathrm{sd}}$ (V)'
zlabel = "$I_{pc}$ (nA)"
ylabel = r'$V_{\mathrm{g}}$ (V)'
insetcolor = 'w'

linecolor = 'gist_rainbow_r'
mapcolor_name = 'ice'

hlines= np.linspace(-5.5, -7, 5)
vlines = []
hlcolor = ["blue", 'deepskyblue',  'orange', 'red','k']

smooth = 1.8

xlimits = [.4, 2.4]
ylimits = [-8, -4]
xticks = [.5, 1, 1.5, 2]
yticks = [-4, -5, -6, -7, -8]

zlimit = [-.1, .38]
zticks = [0, .15, .3]

#                                        ~~ Data Manipulations ~~

xlen = xval.size
ylen = yval.size

if smooth > 0:
    data = ndimage.gaussian_filter(data, smooth)


#                                         ~~ Plot Population ~~

for i in range(len(hlines)):
    ax[plotno].axhline(hlines[i], linestyle = "dashed", c = hlcolor[i + 0], linewidth = 1, alpha = .8)
for i in vlines:
    ax[plotno].axvline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .5)


#Get colormap and color normalization
cmap = pl.get_cmap(mapcolor_name)
if len(zlimit) == 0:
    cnorm = mpl.colors.Normalize(np.min(data), np.max(data))
else:
    cnorm = mpl.colors.Normalize(zlimit[0], zlimit[1])
#define the edges of the map
xt = [xval[0], xval[-1], yval[0], yval[-1]]

ax[plotno].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[barno].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = 'lower')
ax[barno].set_xticks([])
ax[barno].yaxis.tick_right()
ax[barno].yaxis.set_label_position('right')
ax[barno].set_yticks(zticks,zticks, color = 'w')


#Sets edges of image
if xlimits == []:
    None
else:
    ax[plotno].set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[plotno].set_ylim(ylimits[0], ylimits[1])

ax[plotno].set_xticks(xticks)
ax[plotno].set_yticks(yticks)
ax[plotno].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[plotno].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax[barno].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = insetcolor, labelcolor = insetcolor)
ax[barno].yaxis.label.set_color('white')
ax[barno].spines['top'].set_color('white')
ax[barno].spines['bottom'].set_color('white')
ax[barno].spines['left'].set_color('white')
ax[barno].spines['right'].set_color('white')
ax[plotno].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[plotno].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax[barno].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = insetcolor)



'''
----------------------------------------------------        Fig 3d -- Ipc-Vsd Curves     -----------------------
'''

plotno = 3
# barno = 6
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
pval = []
for i in range(powdata[0,0,:].size):
    pval.append(np.mean(powdata[:,:,i]))
pval = np.asanyarray(pval)

#load data map
data = np.load(join(path, datafile_name))
data = data * float(info['Pre-Amp Gain']) * 1e9


#                                     ~~ Data and Panel Parameters ~~


data = data[:,:,-14] - np.mean(data[:,:,0])
data = data + .008

xlabel = r'$V_{\mathrm{sd}}$ (V)'
zlabel = r'$V_{\mathrm{g}}$ (V)'
ylabel = r'$I_{\mathrm{pc}}$ (nA)'

save = True
show = True

mapcolor_name = 'OrRd' #Color of map
barcolor = 'OrRd_r'

colorpadding = 1
zerocolor = False
contours = False

hlines = []
vlines = []

smooth = 1.5



xlimits = [xval[0], xval[-1]]
xticks = [.5, 1, 1.5, 2 , 2.5]
ylimits = []
# yticks = np.arange(0, ylimits[1], 5)
yticks = []
# xlimits = [0]
# ylimits = [0]
zlimit = np.linspace(-5.5, -7, 5)
zticks = [0, -1, -2, -3, -4]
zticks = []
clr = ['k', 'b', 'r']
clr = []


#                                        ~~ Data Manipulations ~~


xlen = xval.size
ylen = yval.size

if smooth > 0:
    data = ndimage.gaussian_filter(data, smooth)


#                                         ~~ Plot Population ~~

for i in hlines:
    ax[plotno].axhline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .8)
for i in vlines:
    ax[plotno].axvline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .5)

if clr == []:
    cmap = pl.get_cmap(mapcolor_name, len(zlimit) + colorpadding)
    clr = cmap(range(len(zlimit) + colorpadding))

x = xval
for i in range(len(zlimit)):
    zindex = np.searchsorted(yval, zlimit[i])
    y = data[zindex,:]

    ax[plotno].plot(x, y, linewidth = linewidth, c = hlcolor[i+0])

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]

#Sets edges of image
if xlimits == []:
    None
else:
    ax[plotno].set_xlim(xlimits[0], xlimits[1])
    ax[plotno].set_xticks(xticks)
    ax[plotno].tick_params(axis='x', labelsize = ticksize)

if ylimits == []:
    None
else:
    ax[plotno].set_ylim(ylimits[0], ylimits[1])

ax[plotno].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'k')
ax[plotno].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'k')
# ax[-3].tick_params(axis='y', labelsize = insetfontsize, direction = 'in', color = 'k', labelcolor = 'k')
ax[plotno].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[plotno].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
# ax[-3].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = 'k')


'''
----------------------------------------------------        Fig 3e -- Interactions     -----------------------
'''

plotno = 6
barno = 7

#                                            ~~ Load Data ~~


path = join(dirname(__file__), 'dev2_photocurrent')
name = '2021_03_25_51_powercube_twobody3_'
# name = 'dev2_interaction-map_twobody'

savefile = join(path, name + ".npz")
if exists(savefile):
    f = np.load(savefile, allow_pickle=True)

    alpha = f['alpha']
    beta = f['beta']
    alphaerror = f['alphaerror']
    betaerror = f['betaerror']
    rsq = f['rsq']
    xval = f['xval']
    yval = f['yval']
    # rfi = f['rfi']
    # info = f['info']

data = alpha

#                                     ~~ Data and Panel Parameters ~~


xlabel = r'$V_{\mathrm{sd}}$ (V)'
zlabel = r'$ \gamma$ / $(\alpha + \beta)$ (arb. units)'
ylabel = r'$V_{\mathrm{g}}$ (V)'
ylabel = ""
insetcolor = 'k'

mapcolor_name = 'OrRd_r'
linecolor = 'gist_rainbow_r'

hlines = [ -5.5, -6, -6.5, -7]
vlines = []
colorpad = 0
# hlcolor = color_meline(linecolor, len(hlines) + colorpad)
hlcolor = ["blue", 'deepskyblue', 'orange', 'red']

smooth = .3
medianfilter = 2
errortrim = 0.65
upperlimit = 4e6

xlimits = [.4, 2.4]
ylimits = [-8, -4]
xticks = [ .5, 1, 1.5, 2 ]
yticks = [-8, -4, 0, 4, 8]
yticks = ([-4, -5, -6, -7, -8])
ticklabels = (['','','','',''])
zlimit = [3, 15]
zticks = ([np.log(10e2), np.log(10e4), np.log(10e6)], ['$10^{2}$', '$10^{4}$', '$10^{6}$'])

zlog = True

#                                        ~~ Data Manipulations ~~


xlen = xval.size
ylen = yval.size

if errortrim > 0:
    for i in range(xlen):
        for p in range(ylen):
            if rsq[p,i] < errortrim:
                data[p,i] = 0
            if data[p,i] < 0:
                data[p,i] = 0
            if data[p,i] == np.inf or data[p,i] == -np.inf:
                data[p,i] = 0

if upperlimit > 0:
    for i in range(xlen):
        for p in range(ylen):
            if np.abs(data[p,i]) > upperlimit:
                data[p,i] = 0

if np.min(xval) == np.max(xval):
    xval = np.linspace(xval[0] - 1, xval[0] + 1, xlen)
if np.min(yval) == np.max(yval):
    yval = np.linspace(yval[0] - 1, yval[0] + 1, ylen)

if smooth > 0:
    data = ndimage.gaussian_filter(data, smooth)

if medianfilter > 0:
    data = ndimage.median_filter(data, medianfilter)

xlen = xval.size
ylen = yval.size

if zlog:
    data = data + .01
    data = np.log(data)


#                                         ~~ Plot Population ~~

if zerocolor:
    zmax, zmin = np.max(data), np.min(data)
    if np.abs(zmax) > np.abs(zmin):
        zmin = -zmax
    else:
        zmax = -zmin

if zlimit == []:
    zlimit = [np.min(data), np.max(data)]

#Get colormap and color normalization
cmap = pl.get_cmap(mapcolor_name)
if len(zlimit) == 0:
    cnorm = mpl.colors.Normalize(np.min(data), np.max(data))
else:
    cnorm = mpl.colors.Normalize(zlimit[0], zlimit[1])
#define the edges of the map
xt = [xval[0], xval[-1], yval[0], yval[-1]]

ax[plotno].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[barno].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = 'lower')
ax[barno].set_xticks([])
ax[barno].set_yticks(zticks[0],zticks[1], color = 'k')


if contours:
    ax[plotno].contour( rf, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')

#Sets edges of image
if xlimits == []:
    None
else:
    ax[plotno].set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[plotno].set_ylim(ylimits[0], ylimits[1])

ax[plotno].set_xticks(xticks)
ax[plotno].set_yticks(yticks)
ax[plotno].set_yticklabels(ticklabels)
ax[plotno].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[plotno].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax[barno].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = insetcolor, labelcolor = insetcolor)
ax[barno].yaxis.label.set_color(insetcolor)
ax[barno].spines['top'].set_color(insetcolor)
ax[barno].spines['bottom'].set_color(insetcolor)
ax[barno].spines['left'].set_color(insetcolor)
ax[barno].spines['right'].set_color(insetcolor)

bbox = dict(boxstyle="square", ec=None, fc="k", alpha=0.5, pad = 0)

ax[plotno].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[plotno].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax[barno].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = insetcolor)

ax[0].set_xticks([])
ax[0].set_yticks([])

ax[0].spines['top'].set_color('w')
ax[0].spines['bottom'].set_color('w')
ax[0].spines['left'].set_color('w')
ax[0].spines['right'].set_color('w')


pl.show()