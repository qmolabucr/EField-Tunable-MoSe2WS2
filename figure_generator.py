'''
Functions that generate the custom plots and colormaps for Figures 2-4
'''

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as pl


def fig_2():
    plotwidth_ = 2.78
    plotheight_ = plotwidth_
    smallplotwidth_ = 1.24
    smallplotheight_ = smallplotwidth_
    marginL_ = .5
    marginR_ = .5
    smallmarginoffset_ = .08
    marginT_ = .1
    marginB_ = .4
    plotsepw_ = 0
    plotseph_ = .35
    smallinsetbarh_ = smallplotheight_ * .45
    smallinsetbarw_ = smallplotwidth_ * .05
    smallinsetmarginB_ = .12
    smallinsetmarginR_ = .05
    maininsetbarh_ = plotheight_ * .36
    maininsetbarw_ = plotwidth_ * .042
    insetmarginB_ = .12
    insetmarginR_ = .12

   #calculates figuresize
    figuresize = ( marginL_ + marginR_ + plotwidth_, marginB_ + marginT_ + plotseph_ + plotheight_ + smallplotheight_) 

   #converts to percent
    marginL = marginL_ / figuresize[0]
    marginR = marginR_ / figuresize[0]
    smallmarginoffset = smallmarginoffset_ / figuresize[0]
    marginT = marginT_ / figuresize[1]
    marginB = marginB_ / figuresize[1]
    plotsepw = plotsepw_ / figuresize[0]
    plotseph = plotseph_ / figuresize[1]
    plotwidth = plotwidth_ / figuresize[0]
    plotheight = plotheight_ / figuresize[1]
    smallplotwidth = smallplotwidth_ / figuresize[0]
    smallplotheight = smallplotheight_ / figuresize[1]
    maininsetbarh = maininsetbarh_ / figuresize[1]
    maininsetbarw = maininsetbarw_ / figuresize[0]
    smallinsetbarh = smallinsetbarh_ / figuresize[1]
    smallinsetbarw = smallinsetbarw_ / figuresize[0]
    smallinsetmarginB = smallinsetmarginB_ / figuresize[1]
    smallinsetmarginR = smallinsetmarginR_ / figuresize[0]
    insetmarginB = insetmarginB_ / figuresize[1]
    insetmarginR = insetmarginR_ / figuresize[0]

    ax = []

    #top right axes
    left, bottom = marginL + plotwidth - smallplotwidth + smallmarginoffset, marginB + plotheight + plotseph
    ax.append([left, bottom, smallplotwidth, smallplotheight])

    #main axes
    left, bottom = marginL, marginB
    ax.append([left, bottom, plotwidth, plotheight])

    #top right colorbar
    left, bottom = marginL + plotwidth + smallmarginoffset - smallinsetmarginR - smallinsetbarw, marginB + plotheight + plotseph + smallinsetmarginB
    ax.append([left, bottom, smallinsetbarw, smallinsetbarh])

    #main colorbar
    left, bottom = marginL + plotwidth - insetmarginR - maininsetbarw, marginB + insetmarginB
    ax.append([left, bottom, maininsetbarw, maininsetbarh])

    #generates figure
    axs = []
    fig = pl.figure(num = 0, figsize=figuresize)
    for a in ax:
        axs.append(fig.add_axes(a))
    axs[0].set_xticks([])
    axs[0].set_yticks([])

    # return axes, figuresize
    return axs, fig


def fig_3():
    scalemod = 2
    totalwidth = 6.3 / scalemod
    marginL_ = .8/ scalemod
    marginR_ = .2/ scalemod
    marginT_ = .5/ scalemod
    marginB_ = .6/ scalemod
    smallplotsep_ = .2/ scalemod
    plotsepsmall_ = .1/ scalemod
    plotwidth_ = (totalwidth - plotsepsmall_) / 2
    plotheight_ = plotwidth_
    plotseph_ = .6/ scalemod
    smallplotheight_ = (plotheight_ - smallplotsep_)/2
    smallplotwidth_ = smallplotheight_
    smallmargin_ = ((plotwidth_ - smallplotwidth_)/2) - .2
    insetplotheight_ = smallplotheight_ * 1.2
    insetplotwidth_ = smallplotwidth_ * 1.2
    insetmargin_ = .1/ scalemod
    insetbarh_ = smallplotheight_ * .4
    insetbarw_ = smallplotwidth_ * .06
    maininsetbarh_ = plotheight_ * .4
    maininsetbarw_ = plotwidth_ * .06

   #calculates figuresize
    figuresize = ( (plotwidth_ * 2) + marginL_ + marginR_ + plotsepsmall_, (plotheight_ * 2) + (marginB_ + marginT_) + plotseph_ ) 

   #converts to percent
    marginL = marginL_ / figuresize[0]
    marginR = marginR_ / figuresize[0]
    marginT = marginT_ / figuresize[1]
    marginB = marginB_ / figuresize[1]
    plotseph = plotseph_ / figuresize[1]
    smallplotsep = smallplotsep_ / figuresize[1]
    plotsepsmall = plotsepsmall_ / figuresize[0]
    plotwidth = plotwidth_ / figuresize[0]
    plotheight = plotheight_ / figuresize[1]
    smallplotwidth = smallplotwidth_ / figuresize[0]
    smallplotheight = smallplotheight_ / figuresize[1]
    insetplotheight = insetplotheight_ / figuresize[1]
    insetplotwidth = insetplotwidth_ / figuresize[0]
    insetmarginw = insetmargin_ / figuresize[0]
    insetmarginh = insetmargin_ / figuresize[1]
    smallmargin = smallmargin_ / figuresize[0]
    insetbarh = insetbarh_ / figuresize[1]
    insetbarw = insetbarw_ / figuresize[0]
    maininsetbarh = maininsetbarh_ / figuresize[1]
    maininsetbarw = maininsetbarw_ / figuresize[0]

    ax = []

    #top left axes
    left, bottom = marginL + smallmargin, marginB + plotheight + plotseph + smallplotheight + smallplotsep
    ax.append([left, bottom, smallplotwidth, smallplotheight])

    #top left bottom axes
    left, bottom = marginL + smallmargin, marginB + plotheight + plotseph
    ax.append([left, bottom, smallplotwidth, smallplotheight])

    #top left bottom colorbar
    left, bottom = marginL + insetmarginw + smallmargin, marginB + insetmarginh + plotheight + plotseph
    ax.append([left, bottom, insetbarw, insetbarh])

    #top right
    left, bottom = marginL + plotwidth + plotsepsmall, marginB + plotheight + plotseph
    ax.append([left, bottom, plotwidth, plotheight])

    #bottom left axes
    left, bottom = marginL , marginB
    ax.append([left, bottom, plotwidth, plotheight])

    #bottom left colorbar
    left, bottom = marginL  + insetmarginw, marginB + plotheight - maininsetbarh - (1*insetmarginh)
    ax.append([left, bottom, maininsetbarw, maininsetbarh])

    #bottom right axes
    left, bottom = marginL +  plotwidth + plotsepsmall, marginB
    ax.append([left, bottom, plotwidth, plotheight])

    #bottom right colorbar
    left, bottom = marginL + plotwidth + plotsepsmall + plotwidth - insetmarginw - maininsetbarw, marginB + plotheight - maininsetbarh - (1*insetmarginh)
    ax.append([left, bottom, maininsetbarw, maininsetbarh])

   #generates figure
    axs = []
    fig = pl.figure(num = 0, figsize=figuresize)
    for a in ax:
        axs.append(fig.add_axes(a))
    axs[0].set_xticks([])
    axs[0].set_yticks([])

    # return ax, figuresize
    return axs, fig


def fig_4():
    #in inches
    plotwidth_ = 3
    plotheight_ = 3
    insetwidth_ = 1.2
    insetheight_ = 1.2
    barwidth_ = 1.4 * .75
    barheight_ = .22* .75
    barleftinset_ = 2.35* .75
    topinset_ = .05* .75
    leftinset_ = .05* .75
    figuremargin_ = 1.5* .75
    marginR_ = 0.05
    marginL_ = 0.4
    marginT_ = 0.4
    marginB_ = 0.4

    #Calculates figuresize
    figuresize = ((plotwidth_) + (marginR_ + marginL_), (plotheight_ *2) + (marginT_ + marginB_))

    #converting to percents
    plotwidth = plotwidth_ / figuresize[0]
    plotheight = plotheight_ / figuresize[1]
    insetwidth = insetwidth_ / figuresize[0]
    insetheight = insetheight_ / figuresize[1]
    barwidth = barwidth_ / figuresize[0]
    barheight = barheight_ / figuresize[1]
    barleftinset = barleftinset_ / figuresize[0]
    topinset = topinset_ / figuresize[1]
    leftinset = leftinset_ / figuresize[0]
    figuremarginR = marginR_ / figuresize[0]
    figuremarginL = marginL_ / figuresize[0]
    figuremarginT = marginT_ / figuresize[1]
    figuremarginB = marginB_ / figuresize[1]    

    ax = []

    #top axes
    left, bottom = figuremarginL, figuremarginB + plotheight
    ax.append([left, bottom, plotwidth, plotheight])

    #inset
    left, bottom = figuremarginL + leftinset, figuremarginB + plotheight + plotheight - insetheight - topinset
    ax.append([left, bottom, insetwidth, insetheight])

    #low axes
    left, bottom = figuremarginL, figuremarginB
    ax.append([left, bottom, plotwidth, plotheight])

    #colorbar
    left, bottom = figuremarginL + barleftinset, figuremarginB + plotheight + plotheight - barheight - topinset
    ax.append([left, bottom, barwidth, barheight])

   #generates figure
    axs = []
    fig = pl.figure(num = 0, figsize=figuresize)
    for a in ax:
        axs.append(fig.add_axes(a))
    axs[0].set_xticks([])
    axs[0].set_yticks([])

    # return ax, figuresize
    return axs, fig


#generates the custom colomaps used in figure 2 and 3

names = []
dicts = []

cdict = {'red':   ((0.0,  0.0, 0.0),
                   (0.1,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (0.9,  .9, .9),
                   (1.0,  1.0, 1.0)),

         'green': ((0.0,  0.0, 0.0),
                   (0.3, 0.0, 0.0),
                   (0.7, 1.0 , 1.0),
                   (1.0,  1.0, 1.0)),

         'blue':  ((0.0,  0.0, 0.0),
                   (0.1,  0.1, 0.1),
                   (0.5,  1.0, 1.0),
                   (0.9,  1.0, 1.0),
                   (1.0,  1.0, 1.0))}

names.append('ice_highconstrast')
dicts.append(cdict)

cdict = {'red':   ((0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (1.0,  1.0, 1.0)),

         'green': ((0.0,  0.0, 0.0),
                   (0.3, 0.0, 0.0),
                   (0.7, 1.0 , 1.0),
                   (1.0,  1.0, 1.0)),

         'blue':  ((0.0,  0.0, 0.0),
                   (0.5,  1.0, 1.0),
                   (1.0,  1.0, 1.0))}

names.append('ice')
dicts.append(cdict)

cdict = {'red':   ((0,  1.0, 1.0),
                   (0.25,  1.0, 1.0),
                  (0.5,  0.0, 0.0),
                   (0.75,  0.0, 0.0),
                   (1.0,  1.0, 1.0)),

         'green': ((0, 1.0, 1.0),
                   (0.15, 1.0, 1.0),
                   (0.35, 0.0, 0.0),
             
                  (0.5,  0.0, 0.0),
                   (0.65, 0.0, 0.0),
                   (0.85, 1.0 , 1.0),
                   (1.0,  1.0, 1.0)),

         'blue':  (
                  (0,  1.0, 1.0),
                  (0.25,  0.0, 0.0),
                  (0.5,  0.0, 0.0),
                   (0.75,  1.0, 1.0),
                   (1.0,  1.0, 1.0))}

names.append('double_ice')
dicts.append(cdict)

cdict = {'red':   ((0,  .95, .95),
                   (0.25,  .95, .95),
                  (0.5,  0.0, 0.0),
                   (0.75,  0.0, 0.0),
                   (1.0,  .95, .95)),

         'green': ((0, .95, .95),
                   (0.15, .95, .95),
                   (0.35, 0.0, 0.0),
                  (0.5,  0.0, 0.0),
                   (0.65, 0.0, 0.0),
                   (0.85, .95 , .95),
                   (1.0,  .95, .95)),

         'blue':  ((0,  .95, .95),
                  (0.25,  0.0, 0.0),
                  (0.5,  0.0, 0.0),
                   (0.75,  .95, .95),
                   (1.0,  .95, .95))}

names.append('double_ice2')
dicts.append(cdict)
sq_mod = .5
cdict = {'red':   ((0,  .2, .2),
                  (0.15*sq_mod,  0, 0),
                  (0.25*sq_mod,  0.0, 0.0),
                  (0.5*sq_mod,  1.0, 1.0),
                  (0.5,  0.5, 0.5),
                  (1.0,  1, 1)),

         'blue': ((0,  .6, .6),
                  (0.15*sq_mod,  1, 1),
                  (0.25*sq_mod,  0.4, 0.4),
                  (0.5*sq_mod,  .6, .6),
                  (0.5,  0.33, 0.33),
                  (1.0,  1, 1)),

         'green':  ((0,  .2, .2),
                  (0.25*sq_mod,  0.8, 0.8),
                  (0.5*sq_mod,  1, 1),
                  (0.5,  0.36, 0.36),
                  (1.0,  1, 1))}

names.append('terrain_spec')
dicts.append(cdict)

for i, name in enumerate(names):
    pl.register_cmap(cmap=LinearSegmentedColormap(name, dicts[i]))
