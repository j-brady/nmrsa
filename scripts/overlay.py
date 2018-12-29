#!/Users/jacobbrady/virtual_envs/py35/bin/python
import sys
import os
import argparse
import nmrglue as ng
import numpy as np
import pandas as pd


if sys.version_info[0]==2:
    import matplotlib
    matplotlib.use("TkAgg")
    print("using python 2")

import matplotlib.pyplot as plt
from matplotlib.cm import Blues_r,get_cmap

from overlay_peaks import Hmqc


def plot_overlay(ax,hmqc,extent=[],cm=Blues_r,contour_start=1e7,
        contour_factor=1.2,contour_num=1,label="",**kwargs):
    """ Function to plot spectral overlays
        
        Function arguments:
        ax -- matplotlib axis object
        hmqc -- Hmqc class object
        peak -- pandas Series containing chem shifts
        extent -- list containing desired ppm limit eg. [x0,x1,y0,y1]
        cm -- matplotlib colormap
        contour_start -- start contouring at this intensity value
        contour_factor -- default 1.2
        contour_num -- default 1
        label -- label for plot

    """
    # convert ppm to points
    if len(extent)<4:
        ppm_13c_0, ppm_13c_1 = hmqc.uc_13c.ppm_limits()
        ppm_1h_0, ppm_1h_1 = hmqc.uc_1h.ppm_limits()
        ppm_1h = hmqc.uc_1h.ppm_scale()
        ppm_13c = hmqc.uc_13c.ppm_scale()
        slice = hmqc.data
        x1,x2,y1,y2 = ppm_1h_0,ppm_1h_1,ppm_13c_0,ppm_13c_1

#        minHpts = hmqc.uc_1h("%f ppm"%extent[0])
#        maxHpts = hmqc.uc_1h("%f ppm"%extent[1])
#        minCpts = hmqc.uc_13c("%f ppm"%extent[2])
#        maxCpts = hmqc.uc_13c("%f ppm"%extent[3])
#
#    if minHpts>maxHpts:
#        minHpts,maxHpts=maxHpts,minHpts
#    if minCpts>maxCpts:
#        minCpts,maxCpts=maxCpts,minCpts
#    slice = hmqc.data[minCpts:maxCpts+1,minHpts:maxHpts+1]
#
#    # calculate contour levels
#    x1,x2,y1,y2 = extent
#    # estimate start contour
#    contour_start = estimate_contour_start(hmqc,peak)
    cl = contour_start * contour_factor ** np.arange(contour_num) 
    if "colors" in kwargs.keys():
        colors = kwargs["colors"]
        ax.contour(slice, cl, 
                    extent=(x1,x2,y1,y2),colors=(colors,),
                    linewidths=0.5)
                    #,**kwargs)
 #       ax.plot(h1,c13,"o",color=kwargs["colors"])
        ax.plot([],[],"-",label=label,color=kwargs["colors"])
    else:
        ax.contour(slice, cl, cmap=cm, 
                    extent=(x2,x1,y2,y1),**kwargs)
        # hack for legend
#        ax.plot(h1,c13,"o",color=cm(1))
        ax.plot([],[],"-",label=label,color=cm(1))
    ax.set_ylim(y1,y2)
    ax.set_xlim(x1,x2)
    #ax.invert_yaxis()
    #ax.invert_xaxis()

def estimate_contour_start(hmqc,peak,extent=[5,5],cutoff=0.6):
    """ Estimate contour_start parameter for contour plot 
    
        Arguments:
        hmqc -- Hmqc class object
        peak -- pandas Series object containing peak information
        extent -- half width of box used to find maximum peak height
        cutoff -- fraction of peak height at which contours start

        Returns:
        cs -- contour start value
    """
    data = hmqc.data
    x,y = hmqc.uc_1h("%f ppm"%float(peak["1H"])),hmqc.uc_13c("%f ppm"%float(peak["13C"]))
    slice = data[y-extent[0]:y+extent[0]+1,x-extent[1]:x+extent[1]+1]
    maximum = slice.max()
    cs = maximum*cutoff
    #print(maximum,cs)
    return cs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot overlays of NMR (HMQC) spectra using matplotlib')
    parser.add_argument('-i','--input',nargs="+",help="Filenames of spectra - you can type as many as you like",type=str)
    parser.add_argument('-l','--labels',nargs="+",help="labels for spectra - you can type as many as you like",type=str)
    parser.add_argument('-o','--out',help="Name of output file",default="test.pdf",type=str)
    parser.add_argument('-c','--color',nargs="*",help="colors for spectra - you can type as many as you like",type=str)
    parser.add_argument('-s','--contourstart',nargs="*",default=1e6,help="contour start for spectra - you can type as many as you like",type=float)
    parser.add_argument('-n','--contournum',nargs="*",default=10,help="contour num for spectra - you can type as many as you like",type=int)
    parser.add_argument('-p','--plot',default=True,type=bool)

    args = parser.parse_args()
    
    outname = args.out
    spec_list = args.input
    lab_list = args.labels
    col_list = args.color
    plot = args.plot
    # contour number
    con_list = args.contournum
    # contour start
    cs_list = args.contourstart

    len_spec_list = len(spec_list)
    
    if not col_list or len(col_list)<len_spec_list:
        col_list = plt.cm.Set1.colors
        
    if not lab_list or len(lab_list)<len_spec_list:
        lab_list = spec_list

    if not con_list or len(con_list)<len_spec_list:
        con_list = [10 for _ in range(len_spec_list)]

    if not cs_list or len(cs_list)<len_spec_list:
        cs_list = [1e7 for _ in range(len_spec_list)]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for s,l,c,cn,cs in zip(spec_list,lab_list,col_list,con_list,cs_list):
        s = Hmqc(s)
        plot_overlay(ax,s,colors=c,
                     contour_start=cs,
                     contour_num=cn,
                     label=l,linewidths=0.1)

    ax.set_xlabel("1H ppm")
    ax.set_ylabel("15N ppm")
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3)
#    plt.tight_layout()
    plt.savefig(outname,bbox_inches='tight')
    if plot:
        plt.show()

    argfile = open("args.txt","w")
    argfile.write(" ".join(sys.argv))
    argfile.close()
