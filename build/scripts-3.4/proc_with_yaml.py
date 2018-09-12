#!/home/jbrady/venvs/py34/bin/python
import os
import re
import sys
import argparse
from time import strftime, localtime

import yaml
import subprocess as sp

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import nmrglue as ng
from scipy.optimize import curve_fit
import pandas as pd

from jinja2 import Environment, PackageLoader

from nmrsa.ng import get_region, normalise, integrate
from nmrsa.fitting import Diffusion
from nmrsa.runNMRPipe import run_pipe_script


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

def drop_inds(array,inds):
    """ Function to remove data points by index
        
        Keyword arguments:
        array -- numpy array
        inds  -- list of indices of data points that you wish to drop
        Returns:
        np.array 
        Note:
        creates 1d boolean array and sets supplied list of indices to False
        then boolean index applied to array
    """
    bools = np.ones(len(array),dtype=np.bool)
    bools[inds]=False
    array = array[bools]
    return array

def run_proc(yaml_dict,g2s,pdf,outfile,read_params=False,delay=9,pulse=53):
    table = {}
    outfile.write("D\tErr\tT_diff\tDelta\tFile\tZGOPTNS\n")
    for k,v in yaml_dict.items():
        k = k.replace("_","-")
        #print("this is k %s" %k)

        run_pipe_script(v["fid.com"],v["dirs"])
        run_pipe_script(v["ft.com"],v["dirs"])
        rows = []
        integrals = pd.DataFrame()
        for ft in v["dirs"]:
            # for reading params
            param_dic,_data = ng.bruker.read(ft)
            ft = os.path.join(ft,v["filename"])

            dic,data = ng.pipe.read(ft)
            uc = ng.pipe.make_uc(dic,data)
            ppm_scale = uc.ppm_scale()
 
            if type(v["start_ppm"]) is list:
                """ If a list is given multiple peaks can be integrated and summed """
                areas = []
                for _start,_end in zip(v["start_ppm"],v["end_ppm"]):
                    if _start < _end:
                       _start, _end = _end, _start

                    start = uc(_start,"ppm")
                    end = uc(_end,"ppm")
                    regions,region_ppm = get_region(data,ppm_scale,start,end)
                    areas.append(np.array([integrate(i) for i in regions]))
                areas = np.array(areas)

            else:
                start_ppm = v["start_ppm"]
                end_ppm = v["end_ppm"]
                # Switch ppm values if wrong way round
                if start_ppm < end_ppm:
                   start_ppm, end_ppm = end_ppm, start_ppm

                start = uc(start_ppm,"ppm")
                end = uc(end_ppm,"ppm")
                regions,region_ppm = get_region(data,ppm_scale,start,end)
            
                areas = np.array([integrate(i) for i in regions])


            # Drop selected points here 
            drop_points = v.get("drop_points",[])
            if type(drop_points) is list and len(drop_points) >= 1:
                areas = drop_inds(areas,drop_points)
                _areas = areas
                _g2s = drop_inds(g2s,drop_points)
                #print(regions.shape,"regions")
                regions = drop_inds(regions,drop_points)
                #print(regions.shape,"regions")
                I0 = np.argmax(areas)
                I0 = areas[I0]
                #print(areas.shape,"SHAPE")
                areas = areas/I0
                #print(areas)
                I0 = areas[0]
            else:
                I0 = np.argmax(areas)
                I0 = areas[I0]
                #I0 = areas[0]
                _areas = areas
                _g2s = g2s
                #print(I0)
                areas = areas/I0
                #print(areas)
                I0 = areas[0]
                #print(I0)
            # initial    
            D = 1e-10
            
            colormap = plt.cm.winter_r
            #colormap = plt.cm.inferno
            plt.rcParams['axes.color_cycle'] = [colormap(i) for i in np.linspace(0, 1, regions.shape[0])]
            fig = plt.figure(figsize=(12,6))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)

            # x,y,wid,height
            ax3 = plt.axes([.27, .55, .20, .3])
            axes = [ax1,ax2,ax3]

            if read_params:
                print("Reading parameters from acqus file")
                T_diff = param_dic['acqus']['D'][delay]
                if T_diff > 1.0 or T_diff < 0.02:
                    raise ValueError("Are you sure you chose the right T_diff value? T_diff = %f s"%T_diff)
                delta  = param_dic['acqus']['P'][pulse]/1e6 # convert from us to s
                if delta > 0.005 or delta < 0.0001:
                    raise ValueError("Are you sure you chose the right delta value? delta = %f s"%delta)
                fit = Diffusion(T_diff=T_diff, delta=delta, Dtype=v["type"], bipolar=v["bipolar"])
            else:
                print("Reading parameters from yaml file")
                T_diff = v["T_diff"]
                delta  = v["delta"]
                fit = Diffusion(T_diff=T_diff, delta=delta, Dtype=v["type"], bipolar=v["bipolar"])
            #print(_g2s.shape,areas.shape,"g2s and areas")
            #print(_g2s)
            popt, pcov = curve_fit(fit.func, _g2s, areas,[I0,D])
            perr = np.sqrt(np.diag(pcov))
            
            ZGOPTNS = param_dic['acqus']['ZGOPTNS']
            tex = " %.3e & $\pm$ %.3e & %.3f & %.3f & %s"%(popt[1], perr[1], T_diff*1000., delta*1000, ft)
            out = " %.3e\t%.3e\t%.6f\t%.6f\t%s\t%s\n"%(popt[1], perr[1], T_diff, delta, ft, ZGOPTNS)
            data_dirs = [ft for _ in _areas]
            integral = pd.DataFrame({"Integral":_areas,"Normalised":areas,"G":np.sqrt(_g2s),"G^2":_g2s,"expt":data_dirs})
            integrals = integrals.append(integral,ignore_index=True)

            rows.append(tex)
            outfile.write(out)
            plot_fit(axes,_g2s,areas,regions,region_ppm,fit.func,popt)

            #plt.suptitle("%s:%s\nD=%8.3e $\pm$%.3e"%(k,ft,popt[1],perr[1]))
            plt.suptitle("T_diff = %f, delta = %f:%s\nD=%8.3e $\pm$%.3e: ZGOPTNS=%s"%(T_diff, delta, ft, popt[1], perr[1], ZGOPTNS))
            pdf.savefig()
            plt.close()
            
            #print("this is k NOW %s"%k)
            table[k]=rows
    outfile.close()
    
    integrals.to_csv("integrals.txt",index=False,sep="\t")
    integrals.to_pickle("integrals.pkl")
    for k,rows in table.items():
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

    #parser.add_argument("--fitpeak",
    #        help="Fit peak volume to gaussian",
    #        action="store_true")

    parser.add_argument("-r","--readAcqus",
            help="read parameters from acqus file",
            action="store_true")

    parser.add_argument("-nl","--nolatex",
            help="don't output latex table",
            action="store_true")

    parser.add_argument("-d","--delay",
            help="number of delay to read e.g. 9 for d9",
            default=9,type=int)

    parser.add_argument("-p","--pulse",
            help="number of pulse to read e.g. 53 for p53",
            default=53,type=int)

    parser.add_argument("-t","--title",
            type=str,help="title for result table",
            default="no title")
    
    parser.add_argument("-o","--outfile",
            type=str,help="name of output file",
            default="output.txt")

    parser.add_argument("-fs","--fitsummary",
            type=str,help="name of fit summary pdf",
            default="fits_summary.pdf")

    args = parser.parse_args()
    yaml_file = args.params
    #difflist = args.difflist
    grad_file = args.gradients
    title = args.title
    readAcqus = args.readAcqus
    delay = args.delay
    pulse = args.pulse
    nolatex = args.nolatex
    of = args.outfile
    outfile = open(of,"w")
    fs = args.fitsummary

    """ Getting gradients """
    grads = load_yaml(grad_file)
    max_strength = grads["max strength"]
    percent = np.array(grads["percentage"])
    gs = max_strength*percent
    g2s = np.square(gs)

    """ Getting params and processing """
    params = load_yaml(yaml_file)
    pdf = PdfPages(fs)
    table = run_proc(params,g2s,pdf,outfile,readAcqus,delay,pulse)
    pdf.close()
    
    """ Making results table """
    env = Environment(loader=PackageLoader('nmrsa', 'templates'))
    temp = env.get_template("diffusion_table.tex")
    oname = "summary.tex"
    out = open(oname,'w')
    out.write(temp.render(tables=table,title=title))
    out.close()

    """ logging """
    user = os.uname()[1]
    t = strftime("%a, %d %b %Y %H:%M:%S +0000", localtime())
    info = ["Run by %s\n"%user,t]
    log = open("proc.log","a")
    log.write("Ran %s\n using the following arguments:\n\n"%__file__)
    log.write("%s\n\n"%" ".join(sys.argv[1:]))
    log.write("proc.yaml:\n%s\n\n"%open(yaml_file,"r").read())
    log.write("gradients.yaml:\n%s\n\n"%open(grad_file,"r").read())
    log.write("output.txt:\n%s\n\n"%open(of,"r").read())
    for i in info:
        log.write(i)
    log.write("\n\n######################################################\n\n")
    log.close()

    """ latex """
    if nolatex:
        pass
    else:
        sp.call("pdflatex %s"%oname,shell=True)
