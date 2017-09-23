#!/Users/jacobbrady/virtual_envs/py35/bin/python
########################################################################
#                                                                      #
# Script for plotting overlays of spectra for inspection of            #
# chemical shift changes.                                              #
#                                                                      #
# Run using a yaml param file eg (params.yml):                         #
#                                                                      #
# general:                                                             #
#        # width and height of spectral window in ppm                  # 
#        width: 0.16                                                   #
#        height: 0.8                                                   #
#                                                                      #
# reference:                                                           #
#        spectrum: "../170914_WT_refold/added.ft2"                     #
#        type: pipe                                                    #
#        peaklist: WT_refold.tab                                       #
#        label: "WT-refold"                                            #
#        contour_start: 3e6                                            #
#        contour_num: 10                                               #
#        cm: Blues_r                                                   #
#                                                                      #
# comparison: # list of spectra to compare                             #
#        - spectrum: "170909_ILVM/bruker800/added.ft2"                 #
#          type: pipe                                                  #
#          peaklist: WT.tab                                            #
#          label: "WT"                                                 #
#          contour_start: 1e7                                          #
#          contour_num: 10                                             #
#          cm: Reds_r                                                  #
#          show_hz: False # whether or not to plot shift change        #
#                                                                      #
########################################################################
#             Written by Jacob Brady September 2017                    #
#                                                                      #
########################################################################
import sys
import yaml
import argparse
import nmrglue as ng
import numpy as np
import pandas as pd

if sys.version_info[0]==2:
    import matplotlib
    matplotlib.use("TkAgg")
    print("using python 2")

import matplotlib.pyplot as plt
from matplotlib.cm import Blues_r,Reds_r,Greens_r

cm_dict = {"Blues_r":Blues_r,"Reds_r":Reds_r,"Greens_r":Greens_r}

def clean_sort(df):
    # clean and sort data by Num
    df = df[df["Num"]<999]
    df = df.sort_values(by="Num")
    return df

class Hmqc:

    def __init__(self,path):

        self.dic, self.data = ng.pipe.read(path)
        self.uc_13c = ng.pipe.make_uc(self.dic, self.data, dim=0)
        # uc dict for 13c
        self.ppm_13c = self.uc_13c.ppm_scale()
        self.ppm_13c_0, self.ppm_13c_1 = self.uc_13c.ppm_limits()
        self.uc_1h = ng.pipe.make_uc(self.dic,self.data, dim=1)
        # uc dict for 1h
        self.ppm_1h = self.uc_1h.ppm_scale()
        self.ppm_1h_0, self.ppm_1h_1 = self.uc_1h.ppm_limits()

def plot_overlay(ax,hmqc,peak,extent,cm=Blues_r,contour_start=1e7,
        contour_factor=1.2,contour_num=1,label=""):
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
    # get peak vals
    h1 = float(peak["1H"])
    c13 = float(peak["13C"])
    ass = peak.H
    num = peak["Num"]
    # convert ppm to points
    minHpts = hmqc.uc_1h("%f ppm"%extent[0])
    maxHpts = hmqc.uc_1h("%f ppm"%extent[1])
    minCpts = hmqc.uc_13c("%f ppm"%extent[2])
    maxCpts = hmqc.uc_13c("%f ppm"%extent[3])
    if minHpts>maxHpts:
        minHpts,maxHpts=maxHpts,minHpts
    if minCpts>maxCpts:
        minCpts,maxCpts=maxCpts,minCpts
    slice = hmqc.data[minCpts:maxCpts+1,minHpts:maxHpts+1]

    # calculate contour levels
    x1,x2,y1,y2 = extent
    cl = contour_start * contour_factor ** np.arange(contour_num) 
    ax.contour(slice, cl, cmap=cm, 
		extent=(x2,x1,y2,y1))
    # hack for legend
    ax.plot(h1,c13,"o",color=cm(1))
    ax.plot([],[],"-",label=label,color=cm(1))
    ax.set_ylim(y2,y1)
    ax.set_xlim(x2,x1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot overlays of single peaks to compare chemical shifts')
    parser.add_argument('--yaml','-y',default="params.yml",help="YAML file containing parameters for running script. See script header for details.",type=str)
    args = parser.parse_args()
    yaml_file = open(args.yaml,"r")
    # extract params from yaml file
    params = yaml.load(yaml_file)
    yaml_file.close()
    general_params = params["general"]
    outname = general_params.get("outname","overlays.pdf")
    reference_params = params["reference"]
    comparison_list = params["comparison"]
    # read reference peak list
    ref_peaks = pd.read_table(reference_params["peaklist"],comment="#")
    ref_peaks = clean_sort(ref_peaks)
    # read reference spectrum 
    ref_spectrum = Hmqc(reference_params["spectrum"])
    # number of plots
    number = len(ref_peaks)
    plots_per_row = 6
    nrows = int(number/plots_per_row)
    if number%plots_per_row == 0:
        nrows = nrows
    else:
        nrows = nrows + 1
    size = 3 # inches
    fig, axes = plt.subplots(nrows, plots_per_row, figsize=(size*plots_per_row, size*nrows))
    axes = axes.ravel()
    # remove axes that are not used
    to_remove = axes[number:]
    [i.remove() for i in to_remove]
    height = general_params["height"]
    width = general_params["width"]
    # compare spectra and overlay
    contour_start = float(reference_params["contour_start"])
    contour_num = int(reference_params["contour_num"])
    label = reference_params["label"]
    cm = cm_dict[reference_params["cm"]]
    # comparison_peaks
    comp_peaks =[pd.read_table(i["peaklist"],comment="#") for i in comparison_list]
    comp_spectra = [Hmqc(i["spectrum"]) for i in comparison_list]
    comp_labels = [i["label"] for i in comparison_list]
    comp_cs = [float(i["contour_start"]) for i in comparison_list]
    comp_cnum = [int(i["contour_num"]) for i in comparison_list]
    comp_cms = [cm_dict[i["cm"]] for i in comparison_list]
    show_hzs = [i["show_hz"] for i in comparison_list]
    for i,ax in zip(ref_peaks.index,axes):
        i = ref_peaks.ix[i]
        i_num = i["Num"]
        i_H = i["H"]
        i_C = i["C"]
       
        c_ppm_1 = i["13C"]
        h_ppm_1 = i["1H"]
        # need to add check to make sure these values are sane since 
        # error will occour if outside bounds of spectrum
        minH,maxH = h_ppm_1-width/2.,h_ppm_1+width/2.
        minC,maxC = c_ppm_1-height/2.,c_ppm_1+height/2.
        extent = [minH,maxH,minC,maxC]
        # plot reference
        #print(type(float(reference_params["contour_start"])))
        plot_overlay(ax,ref_spectrum,i,extent,cm=cm,
                contour_start=contour_start,
                contour_num=contour_num,
                label=label)

        #for c_peak,c_spec,c_lab,c_cs,c_cnum,c_cm in zip_compare:
        for c_peak,c_spec,c_lab,c_cs,c_cnum,c_cm,show_hz in zip(comp_peaks,comp_spectra,comp_labels,comp_cs,comp_cnum,comp_cms,show_hzs):
            c_peak = clean_sort(c_peak)
            peak = c_peak[(c_peak["Num"]==i_num) & (c_peak.H==i_H) & (c_peak.C==i_C)]
            if peak.empty:
                print("No peak for %d"%i_num)
            elif len(peak)>1:
                print("Duplicated for %d"%i_num)
            else:
                dc = float(i["13C"] - peak["13C"])
                dh = float(i["1H"] - peak["1H"])
                dc_hz = float(dc*200)
                dh_hz = float(dh*800)
                c_ppm_2 = peak["13C"]
                h_ppm_2 = peak["1H"]
                comp = np.sqrt(dc_hz**2+dh_hz**2)
                print(peak)
                plot_overlay(ax,c_spec,peak,extent,cm=c_cm,contour_start=c_cs,
                        contour_num=c_cnum,label=c_lab)
            if show_hz:
                ax.text(0.1,0.1,"%.1f Hz"%(comp),transform=ax.transAxes,
                    bbox=dict(facecolor='yellow', alpha=0.5))
        ax.set_title("%d - %s"%(i["Num"],i.H))
        ax.set_ylabel("$^{13}$C ppm")
        ax.set_xlabel("$^{1}$H ppm")
        ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outname)
