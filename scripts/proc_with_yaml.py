#! /usr/bin/python
import os
import re
import argparse

import yaml
import subprocess as sp

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import nmrglue as ng
from scipy.optimize import curve_fit

from jinja2 import Environment, PackageLoader

from NMR.diffusion import normalise, get_region, integrate, FitDiffusion 
from NMR.runNMRPipe import run_pipe_script


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

def run_proc(yaml_dict,g2s,pdf):
    table = {}
    for k,v in yaml_dict.iteritems():
        k = k.replace("_","-")
        print "this is k %s" %k

        run_pipe_script(v["fid.com"],v["dirs"])
        run_pipe_script(v["ft.com"],v["dirs"])
        rows = []
        for ft in v["dirs"]:
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
            fig, axes = plt.subplots(1, 2,figsize=(12,6))

            fit = FitDiffusion(T_diff=v["T_diff"],delta=v["delta"],Dtype=v["type"])
            popt, pcov = curve_fit(fit.func, g2s, areas,[I0,D])
            perr = np.sqrt(np.diag(pcov))
            
            tex = " %.3e & $\pm$ %.3e & %.3f & %.3f & %s"%(popt[1],perr[1],v["T_diff"]*1000.,v["delta"]*1000,ft)
            rows.append(tex)
            plot_fit(axes,g2s,areas,regions,region_ppm,fit.func,popt)

            plt.suptitle("%s:%s\nD=%8.3e $\pm$%.3e"%(k,ft,popt[1],perr[1]))
            pdf.savefig()
            plt.close()
            print "this is k NOW %s"%k
            table[k]=rows
    for k,rows in table.iteritems():
        table[k] = [val.replace("_","-") for val in rows]
    return table

def plot_fit(axes,g2s,areas,regions,region_ppm,func,popt):
    x = np.linspace(g2s.min(),g2s.max())
    ax1,ax2 = axes
    ax2.set_xlim(max(region_ppm),min(region_ppm))
    [ax2.plot(region_ppm,region,label="%d"%g2) for g2,region in zip(g2s,normalise(regions))]        
    ax2.set_xlabel("ppm")
    ax2.set_ylabel("normalised intensity")
    ax2.legend(title="Gradient strength")

    ax1.plot(x,func(x,*popt),g2s,areas,"o")
    ax1.set_ylabel(r"I/I$_{0}$") 
    ax1.set_xlabel(r"G$^2$ $(G^2cm^{-2})$")



if __name__ == "__main__":
    """ Argument parser """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--params",
            type=str,help="Yaml file containing file names and parameters for processing.",
            default="proc.yaml")

    parser.add_argument("-g","--gradients",
            type=str,help="Yaml file containing gradient strengths.",
            default="gradients.yaml")

    parser.add_argument("-t","--title",
            type=str,help="title for result table",
            default="no title")

    args = parser.parse_args()
    yaml_file = args.params
    grad_file = args.gradients
    title = args.title

    """ Getting gradients """
    grads = load_yaml(grad_file)
    max_strength = grads["max strength"]
    percent = np.array(grads["percentage"])
    gs = max_strength*percent
    g2s = np.square(gs)

    """ Getting params and processing """
    params = load_yaml(yaml_file)
    pdf = PdfPages("fits_summary.pdf")
    table = run_proc(params,g2s,pdf)
    pdf.close()
    
    """ Making results table """
    env = Environment(loader=PackageLoader('NMR', 'templates'))
    temp = env.get_template("diffusion_table.tex")
    oname = "summary.tex"
    out = open(oname,'w')
    out.write(temp.render(tables=table,title=title))
    out.close()
    sp.call("pdflatex %s"%oname,shell=True)
