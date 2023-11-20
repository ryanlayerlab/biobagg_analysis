#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
import pylab
import random
from optparse import OptionParser


from matplotlib import rcParams
rcParams['font.family'] = 'Arial'

parser = OptionParser()

parser.add_option("--abvline",
                  dest="abvline",
                  help="Vertincal line")

parser.add_option("--noyaxis",
                  dest="noyaxis",
                  action="store_true",default=False,
                  help="Hide y-axsis")

parser.add_option("--tick_line_length",
                  dest="tick_line_length",
                  type=float,
                  default=2,
                  help="Tick line width")

parser.add_option("--tick_line_width",
                  dest="tick_line_width",
                  type=float,
                  default=0.5,
                  help="Tick line width")

parser.add_option("--axis_line_width",
                  dest="axis_line_width",
                  type=float,
                  default=0.5,
                  help="Axis line width")

parser.add_option("--axis_label_size",
                  dest="axis_label_size",
                  type=int,
                  default=8,
                  help="Axis label font size")

parser.add_option("--tick_label_size",
                  dest="tick_label_size",
                  type=int,
                  default=8,
                  help="Axis tick label font size")

parser.add_option("--density",
                  dest="density",
                  action="store_true",default=False,
                  help="Plot density")



parser.add_option("--xticks",
                  dest="xticks",
                  help="CSV ints to tick and label")

parser.add_option("--xtick_names",
                  dest="xtick_names",
                  help="CSV of xtick lables")


parser.add_option("--numyticks",
                  dest="numyticks",
                  help="Number of Y ticks")

parser.add_option("-t",
                  "--title",
                  dest="title",
                  help="Title")

parser.add_option("-x",
                  "--xlabel",
                  dest="xlabel",
                  help="X axis label")

parser.add_option("-y",
                  "--ylabel",
                  dest="ylabel",
                  help="Y axis label")

parser.add_option("-o",
                  "--output_file",
                  dest="output_file",
                  help="Data file")

parser.add_option("-b",
                  "--bins",
                  dest="bins",
                  default=10,
                  help="Number of bins or csv of bins")

parser.add_option("--x_max",
                  dest="max_x",
                  type="float",
                  help="Max x value")


parser.add_option("--x_min",
                  dest="min_x",
                  type="float",
                  help="Min x value")

parser.add_option("--y_max",
                  dest="max_y",
                  type="float",
                  help="Max y value")


parser.add_option("--y_min",
                  dest="min_y",
                  type="float",
                  help="Min y value")

parser.add_option("--color",
                  dest="color",
                  default=None,
                  help="Bar color")

parser.add_option("-c",
                  "--column",
                  dest="col",
                  type="int",
                  default="0",
                  help="Column in the data")

parser.add_option("-d",
                  "--delim",
                  dest="delim",
                  default="\t",
                  help="Field delimiter")

parser.add_option("-l",
                  "--ylog",
                  action="store_true",default=False,
                  dest="ylog",
                  help="Set y-axis to be log scale")

parser.add_option("--xlog",
                  action="store_true",default=False,
                  dest="xlog",
                  help="Set x-axis to be log scale")



parser.add_option("--x_sci",
                  action="store_true",default=False,
                  dest="x_sci",
                  help="Use scientific notation for x-axis")

parser.add_option("--y_sci",
                  action="store_true",default=False,
                  dest="y_sci",
                  help="Use scientific notation for y-axis")

parser.add_option("--width",
                  dest="width",
                  type="float",
                  default=5,
                  help="Figure width")

parser.add_option("--height",
                  dest="height",
                  type="float",
                  default=5,
                  help="Figure height")

parser.add_option("--black",
                  action="store_true", 
                  default=False,
                  dest="black",
                  help="black background")



(options, args) = parser.parse_args()
if not options.output_file:
    parser.error('Output file not given')

Y=[]
for l in sys.stdin:
    a = l.rstrip().split(options.delim)
    if len(a) == 1:
        if len(a[options.col]) != 0 :
            Y.append(float(a[options.col]))

matplotlib.rcParams.update({'font.size': 12})
#fig = matplotlib.pyplot.figure(figsize=(10,5),dpi=300)

#fig = matplotlib.pyplot.figure(figsize=(options.width,options.height),dpi=300)

if options.black:
    fig = matplotlib.pyplot.figure(\
            figsize=(options.width,options.height),\
            dpi=300,\
            facecolor='black')
else:
    fig = matplotlib.pyplot.figure(\
            figsize=(options.width,options.height),\
            dpi=300)



fig.subplots_adjust(wspace=.05,left=.01,bottom=.01)

x_max = max(Y)
x_min = min(Y)

if options.max_x:
    x_max = options.max_x
if options.min_x:
    x_min = options.min_x

if options.black:
    ax = fig.add_subplot(1,1,1,facecolor='k')
else:
    ax = fig.add_subplot(1,1,1)

bins = None
if ',' in str(options.bins):
    bins = [int(x) for x in options.bins.split(',')]
else:
    bins = int(options.bins)

h = ax.hist(Y, \
            bins,
            density=options.density,
            log=options.ylog, \
            histtype='bar', \
            rwidth=0.8, \
            color=options.color)

print(h)

if options.xlog:
    ax.set_xscale('log')

#labels, counts = np.unique(Y, return_counts=True)
#print labels
#ax.bar(labels, counts, align='center')
#plt.gca().set_xticks(labels)


if options.max_x:
    ax.set_xlim(xmax=options.max_x)
if options.min_x:
    ax.set_xlim(xmin=options.min_x)
if options.max_y:
    ax.set_ylim(ymax=options.max_y)
if options.min_y:
    ax.set_ylim(ymin=options.min_y)

if options.x_sci:
    formatter = matplotlib.ticker.ScalarFormatter()
    formatter.set_powerlimits((-2,2))
    ax.xaxis.set_major_formatter(formatter)

if options.y_sci:
    formatter = matplotlib.ticker.ScalarFormatter()
    formatter.set_powerlimits((-2,2))
    ax.yaxis.set_major_formatter(formatter)

if options.xlabel:
    ax.set_xlabel(options.xlabel, fontsize=options.axis_label_size)

if options.ylabel:
    ax.set_ylabel(options.ylabel, fontsize=options.axis_label_size)

if options.title:
    ax.set_title(options.title)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

if options.xticks:
    xticks = [int(x) for x in options.xticks.split(',')]
    xmajorlocator = matplotlib.ticker.FixedLocator(xticks)
    ax.xaxis.set_major_locator(xmajorlocator)

if options.xtick_names:
    xtick_locs = []
    xtick_names = []

    i = 0
    for xtick_name in options.xtick_names.split(','):
        if xtick_name != '':
            xtick_locs.append(i)
            xtick_names.append(xtick_name)
        i+=1

    xmajorlocator = matplotlib.ticker.FixedLocator(xtick_locs)
    ax.xaxis.set_major_locator(xmajorlocator)
    ax.set_xticklabels(xtick_names)

if options.black:
    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.title.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.xaxis.label.set_color('white')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')

ax.tick_params(axis='both',
               which='major',
               labelsize=options.axis_label_size,
               width=options.tick_line_width,
               length=options.tick_line_length)

ax.spines['bottom'].set_linewidth(options.axis_line_width)
ax.spines['left'].set_linewidth(options.axis_line_width)

if options.noyaxis:
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_ticklabels([])
    ax.set_yticks([], minor=True)
    ax.set_yticks([])

if options.abvline:
    #position, text, color, size = args.abvline.split(',')
    ax.axvline(float(options.abvline), color='red')
   # ax.text(float(position) + \
   #         (ax.get_xlim()[1]-ax.get_xlim()[0])*0.05,\
   #         ax.get_ylim()[1],
   #         text,
   #         va='top',\
   #         fontsize=float(size),\
   #         color=color)


if options.black:
    matplotlib.pyplot.savefig(options.output_file,bbox_inches='tight',\
            facecolor=fig.get_facecolor(),\
              transparent=True)
else:
    matplotlib.pyplot.savefig(options.output_file,bbox_inches='tight')
