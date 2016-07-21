#! /usr/bin/python
import os
import re
import argparse

import yaml
import subprocess as sp

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import numpy as np
import nmrglue as ng
from scipy.optimize import curve_fit

from jinja2 import Environment, PackageLoader

from nmrsa.ng import get_region, normalise, integrate
from nmrsa.fitting import Diffusion
from nmrsa.runNMRPipe import run_pipe_script

#sns.set_style("whitegrid")
sns.set_style("ticks")

def load_yaml(yaml_file):
    """ reads files containing YAML and converts to dictionary """
    f = open(yaml_file)
    y = f.read()
    return yaml.load(y)

def write_yaml(dict_obj,fname="output.yaml"):
    with open(fname,"a") as f:
        f.write( yaml.dump(f,default_flow_style=False))

def write_latex(string,fname="output.tex"):
    with open(fname,"a") as f:
        f.write("%s \n"% string)

def run_proc(yaml_dict,g2s,pdf,outfile):
    table = {}
    outfile.write("D\tErr\tT_diff\tDelta\tFile\tZGOPTNS\n")
    for k,v in yaml_dict.iteritems():
        k = k.replace("_","-")
        print "this is k %s" %k

        run_pipe_script(v["fid.com"],v["dirs"])
        run_pipe_script(v["ft.com"],v["dirs"])
        rows = []
        for ft in v["dirs"]:
            # for reading params
            param_dic,_data = ng.fileio.bruker.read(ft)
            ft = os.path.join(ft,v["filename"])

            dic,data = ng.pipe.read(ft)
            uc = ng.pipe.make_uc(dic,data)
            ppm_scale = uc.ppm_scale()
            start_ppm = v["start_ppm"]
            end_ppm = v["end_ppm"]
            start = uc(start_ppm,"ppm")
            end = uc(end_ppm,"ppm")
            regions,region_ppm = get_region(data,ppm_scale,start,end)
            
            areas = np.array([integrate(i) for i in regions])
            I0 = areas[0]
            print I0
            areas = areas/I0
            I0 = areas[0]
            D = 1e-10
            
            colormap = plt.cm.winter_r
            plt.rcParams['axes.color_cycle'] = [colormap(i) for i in np.linspace(0, 1, regions.shape[0])]
            fig = plt.figure(figsize=(12,6))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)

            #ax1 = plt.subplot2grid((4,4), (2, 0), colspan=2,rowspan=2)
            #ax2 = plt.subplot2grid((4,4), (0, 0), colspan=2,rowspan=2)
            #ax3 = plt.subplot2grid((4,4), (1, 0), colspan=1,rowspan=1)
            # x,y,wid,height
            ax3 = plt.axes([.27, .55, .20, .3])
            axes = [ax1,ax2,ax3]

            fit = Diffusion(T_diff=v["T_diff"],delta=v["delta"],Dtype=v["type"],bipolar=v["bipolar"])
            popt, pcov = curve_fit(fit.func, g2s, areas,[I0,D])
            perr = np.sqrt(np.diag(pcov))
            
            tex = " %.3e & $\pm$ %.3e & %.3f & %.3f & %s"%(popt[1],perr[1],v["T_diff"]*1000.,v["delta"]*1000,ft)
            out = " %.3e\t%.3e\t%.6f\t%.6f\t%s\t%s\n"%(popt[1],perr[1],v["T_diff"],v["delta"],ft,param_dic['acqus']['ZGOPTNS'])

            rows.append(tex)
            outfile.write(out)
            plot_fit(axes,g2s,areas,regions,region_ppm,fit.func,popt)

            plt.suptitle("%s:%s\nD=%8.3e $\pm$%.3e"%(k,ft,popt[1],perr[1]))
            pdf.savefig()
            plt.close()
            
            print "this is k NOW %s"%k
            table[k]=rows
    outfile.close()
    for k,rows in table.iteritems():
        table[k] = [val.replace("_","-") for val in rows]
    return table

def plot_fit(axes,g2s,areas,regions,region_ppm,func,popt):
    x = np.linspace(g2s.min(),g2s.max())
    ax1,ax2,ax3 = axes
    ax2.set_xlim(max(region_ppm),min(region_ppm))
    [ax2.plot(region_ppm,region,label="%d"%g2) for g2,region in zip(g2s,normalise(regions))]        
    ax2.set_xlabel("ppm")
    ax2.set_ylabel("normalised intensity")
    ax2.legend(title="$G^2$ - ($G^2 cm^{-2}$)",fontsize="small")

    ax1.plot(x,func(x,*popt),g2s,areas,"o")
    ax1.set_ylabel(r"I/I$_{0}$") 
    ax1.set_xlabel(r"G$^2$ $(G^2cm^{-2})$")

    
    ax3.plot(x,np.log(func(x,*popt)),g2s,np.log(areas),"o")
    ax3.set_ylabel(r"ln(I/I$_{0})$") 
    ax3.set_xlabel(r"G$^2$ $(G^2cm^{-2})$")
    
    sns.despine()

if __name__ == "__main__":
    """ Argument parser """
    parser = argparse.ArgumentParser(prog='PROG', usage='%(prog)s [options]',description="Script for processing psuedo 2D diffusion data written by Jacob Brady and Rui Huang.")

    parser.add_argument("-i","--params",
            type=str,help="Yaml file containing file names and parameters for processing.",
            default="proc.yaml")

    parser.add_argument("-g","--gradients",
            type=str,help="Yaml file containing gradient strengths.",
            default="gradients.yaml")
    
    #parser.add_argument("-dl","--difflist",
    #        type=str,help="Bruker difflist file containing gradient strengths.")

    parser.add_argument("--fitpeak",
            help="Fit peak volume to gaussian",
            action="store_true")

    parser.add_argument("-t","--title",
            type=str,help="title for result table",
            default="no title")
    
    parser.add_argument("-o","--outfile",
            type=str,help="name of output file",
            default="output.txt")

    args = parser.parse_args()
    yaml_file = args.params
    #difflist = args.difflist
    grad_file = args.gradients
    title = args.title
    outfile = open(args.outfile,"w")

    """ Getting gradients """
    grads = load_yaml(grad_file)
    max_strength = grads["max strength"]
    percent = np.array(grads["percentage"])
    gs = max_strength*percent
    g2s = np.square(gs)

    """ Getting params and processing """
    params = load_yaml(yaml_file)
    pdf = PdfPages("fits_summary.pdf")
    table = run_proc(params,g2s,pdf,outfile)
    pdf.close()
    
    """ Making results table """
    env = Environment(loader=PackageLoader('nmrsa', 'templates'))
    temp = env.get_template("diffusion_table.tex")
    oname = "summary.tex"
    out = open(oname,'w')
    out.write(temp.render(tables=table,title=title))
    out.close()
    sp.call("pdflatex %s"%oname,shell=True)
