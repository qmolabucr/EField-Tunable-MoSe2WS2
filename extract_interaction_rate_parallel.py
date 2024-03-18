'''
Takes the photocurrent datacube and extracts the relative rate of two-body interaction (gamma/(alpha+beta)).
Runs in embarrassingly parallel but can be calculated linearly.

Fits data to 
   Ipc = (B * np.log(1+Power*A)) + D
and produces a .npz file with the following arrays:

alpha: 2d array of best fit for A
beta: 2d array of best fit for B
alphaerror: Error in A fit
betaerror: Error in B fit
rsq: R squared fit value
xval: X axis
yval: Y axis

The current settings are not the exact ones used for the map in figure 3e but they produce a similar map. Feel free to play around with it!
'''
# import matplotlib as mpl
# import matplotlib.pyplot as pl
import numpy as np
from scipy import ndimage
import scipy.optimize as op
from os import cpu_count
from os.path import join, dirname
from joblib import Parallel, delayed
import time

cores = cpu_count()

#The function we will fit the data to. Others can be used
# def twobody(x, a, b, d):
#     # i = ((b/a) * np.log(c+x*a)) +d
#     # i = ((b/a) * np.log(1+(x-c)*a)) +d
#     i = (b * np.log(1+x*a)) +d
#     return i

def twobody(x, a, b):
    # i = ((b/a) * np.log(c+x*a)) +d
    # i = ((b/a) * np.log(1+(x-c)*a)) +d
    i = (b * np.log(1+x*a))
    return i

#                                      ~~ Loads photocurrent data ~~

#Path to data directory and file name(s)
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

#load data map and calibrate to pA
data = np.load(join(path, datafile_name))
data = data * float(info['Pre-Amp Gain'])
data = (data - np.min(data)) * 1e12


#                                      ~~ Script Parameters ~~
#Fit function, others can be added
fitfunc = twobody

#Smooth? Not recommended for fitting unless the data is really bad
smooth = 0
smoothpower = 0

#Zero x/y data
zeroset_x = True
zeroset_y = False

#Offset x/y, neccessary for the log fit to work well
x_offset = 0.001
y_offset = 0

#Use index or searchsorted to trim data?
use_power_index = True
trimval = .16
trimindex = 26

#Fit parameters
p0 = [500, .01, .001]
p0 = [500, .01]

trackp0 = False
bounds_twobody2 = ([0, 0, 0], [np.inf, np.inf, np.mean(data[:,:,3])])
bounds_twobody2 = ([0, 0], [np.inf, np.inf])


#                                      ~~ Data Manipulation Prior to Fitting ~~

xlen = xval.size
ylen = yval.size
plen = cval.size


if use_power_index:
    ptrim = trimindex
else:
    ptrim = np.searchsorted(cval, trimval)
# lowptrim = np.searchsorted(cval, trimval)
noise = np.zeros(data[:,:,0].shape)
yrange = np.zeros(data[:,:,0].shape)
for i in range(ylen):
    for p in range(xlen):
        y = data[i, p, :ptrim]
        noise[i,p] = np.abs(np.max(y[-10:]) - np.min(y[-10:]))
        yrange[i,p] = np.abs(np.mean(y[-3:]) - np.mean(y[:3]))

if smooth > 0:
    for i in range(plen):
        data[:,:,i] = ndimage.gaussian_filter(data[:,:,i], smooth)

#                                        ~~ Fit Function We Will Use ~~

def Fit(x, y, fitfunc, p0, bounds):
    try:
        par, pcov = op.curve_fit(fitfunc, x, y, bounds = bounds, p0 = p0, maxfev = 3200)

        ss_residuals = np.sum( (y - fitfunc(x, *par))**2 )
        ss_squares = np.sum( (y - np.mean(y) )**2)
        r_sq = 1 - (ss_residuals/ss_squares)

        alpha = par[0]
        beta = par[1]
        alphaerror = np.sqrt(np.diag(pcov))[0]
        betaerror = np.sqrt(np.diag(pcov))[1]
    except RuntimeError:
        alpha = 0
        beta = 0
        alphaerror = -1
        betaerror = -1
        r_sq = 0

    return alpha, beta, alphaerror, betaerror, r_sq

#                                     ~~ Data Compilation/Preparation ~~

print("Starting fitting %i rows of points" % (ylen))
#compiles Ipc-Power data curves that will then be passed to the parallel fit function

failedfits = 0
rawdata = []

for i in range(ylen):
    for p in range(xlen):

        y = data[i, p, :ptrim]
        # x = powdata[i, p, :ptrim]
        x = cval[:ptrim]

        if zeroset_x:
            x = x - np.min(x) + x_offset
        if zeroset_y:
            y = y - np.min(y) + y_offset
        if smoothpower > 0:
            y = ndimage.gaussian_filter(y, smoothpower)

        rawdata.append([x,y])

#                                               ~~ Fit Proccess ~~

a = time.time()

#parallel process
x = Parallel(n_jobs=cores, backend='loky')(
    delayed(Fit)(rawd[0], rawd[1], fitfunc, p0, bounds_twobody2) for rawd in rawdata)

#Series process
# x = []
# cc = 0
# for rawd in rawdata:
#     x.append(Fit(rawd[0], rawd[1], fitfunc, p0, bounds_twobody2))
    # print(cc)
    # cc+=1

print('Fit %i curves in %.2f seconds' % (len(x), time.time() - a))

#                                          ~~ Tidying up and saving data ~~

alpha = np.zeros(data[:,:,0].shape)
beta = np.zeros(data[:,:,0].shape)
alphaerror = np.zeros(data[:,:,0].shape)
betaerror = np.zeros(data[:,:,0].shape)
rsq =  np.zeros(data[:,:,0].shape)

count = 0
for i in range(ylen):
    for p in range(xlen):
        alpha[i,p] = x[count][0]
        beta[i,p] = x[count][1]
        alphaerror[i,p] = x[count][2]
        betaerror[i,p] = x[count][3]
        rsq [i,p] = x[count][4]
        count+=1
        
print("  ~~ saving ~~  ")


savename = join(path, 'dev2_interaction-map' + "_" + str(fitfunc).split(" ")[1])
print(' Saving as  %s' % (savename))
np.savez(savename, alpha = alpha, beta = beta, alphaerror = alphaerror, betaerror = betaerror, rsq=rsq, xval = xval, yval = yval)

#To unpack the saved file use:
'''

f = np.load(savename, allow_pickle=True)

alpha = f['alpha']
beta = f['beta']
alphaerror = f['alphaerror']
betaerror = f['betaerror']
rsq = f['rsq']
xval1 = f['xval']
yval1 = f['yval']

'''
