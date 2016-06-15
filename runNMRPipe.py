import os
import subprocess as sp


def print_script(name):
    with open(name,'r') as f:
        f = f.read()
        print(f) 

def run_pipe_script(script,dirs):
    """ Run i.communicate() for each element in proc to see shell output """
    for d in dirs:
        print("Processing: %s"%d)
        sp.Popen("%s %s"%(script,d),shell=True).wait()

def run_pipe_2args(script,dirs,arg2):
    """ Run i.communicate() for each element in proc to see shell output """
    for d in dirs:
        print("Processing: %s"%d)
        sp.Popen("%s %s %s"%(script,d,arg2),shell=True).wait()
        
#os.path.join ...
""" |& will pass both stdout and stderr through pipe. Normally stderr is not piped. """ 
if __name__  == "__main__":

    cwd = os.getcwd()
    fid_com = "fid.com"
    xy_com = "xy.com"
    fid = os.path.join(cwd,fid_com)
    xy = os.path.join(cwd,xy_com)
    print(fid)
    print_script(fid)
    print(xy)
    print_script(xy)

    dirs = ["aqueous/913","915","916","917","918","919"]
    run_pipe_script(fid,dirs)
    run_pipe_script(xy,dirs)
    print("All done!")
    addNMR = "comb_spectra.sh"
    print("Adding spectra using %s" % addNMR)
    sp.Popen("./%s" % addNMR,shell=True)
