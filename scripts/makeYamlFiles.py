#!/usr/bin/python
import argparse
import os
from time import strftime, localtime

from numpy import arange
from jinja2 import Environment, PackageLoader


user = os.uname()[1]
t = strftime("%a, %d %b %Y %H:%M:%S +0000", localtime())
info = ["File created by %s"%user,t]

PATH = os.getcwd()
DIRS = filter(os.path.isdir, os.listdir(PATH))
DIRS = [i for i in DIRS if "ser" in os.listdir(os.path.join(PATH,i))]
# Need to figure out how to deal with 1d data set
if len(DIRS)<1:
    raise Exception("No ser files found! Are you sure you downloaded the data?")

""" argparser stuff """
parser = argparse.ArgumentParser(description="Script to generate yaml files for processing diffusion data")
parser.add_argument("-p","--proc",
                    type=str,help="Name for yaml file containing file names and parameters for processing.",
                    default="proc.yaml")
parser.add_argument("-g","--gradients",
                    type=str,help="Name for yaml file containing file names and parameters for processing.",
                    default="gradients.yaml")
parser.add_argument("-ng","--nogradients",
                    help="Do not generate yaml file containing gradients",
                    action="store_true")
parser.add_argument("-d","--dtype",type=str,help="Type of diffusion experiment: single (1Q) or triple (3Q) quantum.",
                    default="1Q",choices=set(("1Q","2Q","3Q")))
parser.add_argument("-b","--bipolar", help="Use flag if bipolar gradients are used",
                    action="store_true")
args = parser.parse_args()

""" jinja2 setup """
env = Environment(loader=PackageLoader('nmrsa', 'templates'))
gradients_temp = env.get_template("gradients.yaml")
proc_temp = env.get_template("proc.yaml")

""" proc.yaml """ 
filename = "test.ft"
fid_com = "fid.com"
ft_com = "ft.com"
dtype = args.dtype
bipolar = args.bipolar
start_ppm = 2.1
end_ppm = 1.9

outname = args.proc 
if outname in os.listdir("./"):
    os.system("cp %s %s.backup"%(outname,outname))
    print("Moved %s to %s.backup"%(outname,outname))

outfile = open(outname,"w")
outfile.write(proc_temp.render(info=info,directories=DIRS,filename=filename,fid_com=fid_com,
        ft_com=ft_com,start_ppm=start_ppm,end_ppm=end_ppm,bipolar=bipolar,type=dtype))

""" gradients.yaml """ 
no_gradients = args.nogradients
outname = args.gradients
if no_gradients:
    pass

else:
    if outname in os.listdir("./"):
        os.system("cp %s %s.backup"%(outname,outname))
        print("Moved %s to %s.backup"%(outname,outname))

    max_strength = 45.0
    percentages = arange(0.1,1,0.1)
    outfile = open(outname,"w")
    outfile.write(gradients_temp.render(info=info,max_strength=max_strength,percentages=percentages))

