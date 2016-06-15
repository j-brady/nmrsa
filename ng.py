import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt
from colormap import viridis,magma,plasma 

""" My nmrglue functions """
class OneD:
    
    def __init__(self,fname="test.ft",spectype="pipe"):
        if spectype == "pipe":

            self.dic,self.data = ng.pipe.read(fname) 
            self.uc = ng.pipe.make_uc(self.dic, self.data)
            self.ppm_scale = self.uc.ppm_scale()

    def get_region(self,start_ppm,end_ppm):

        if start_ppm < end_ppm:
            start_ppm, end_ppm = end_ppm, start_ppm

        start_pts = self.uc(start_ppm,"ppm")
        end_pts = self.uc(end_ppm,"ppm")
        self.region = self.data[start_pts:end_pts+1]
        self.region_ppm = self.ppm_scale[start_pts:end_pts+1]

        return self.region,self.region_ppm


class Pseudo2D:


    def __init__(self,fname="test.ft",spectype="pipe"):

        if spectype == "pipe":

            self.dic,self.data = ng.pipe.read(fname) 
            self.uc = ng.pipe.make_uc(self.dic, self.data)
            self.ppm_scale = self.uc.ppm_scale()

    def get_region(self,start_ppm,end_ppm):

        if start_ppm < end_ppm:
            start_ppm, end_ppm = end_ppm, start_ppm

        #print start_ppm,end_ppm
        start_pts = self.uc(start_ppm,"ppm")
        end_pts = self.uc(end_ppm,"ppm")
        #print start_pts,end_pts
        self.region = self.data[:,start_pts:end_pts+1]
        self.region_ppm = self.ppm_scale[start_pts:end_pts+1]

        #colormap = plt.cm.winter_r
        colormap = viridis
        #colormap = magma
        #colormap = plasma
        plt.rcParams['axes.color_cycle'] = [colormap(k) for k in np.linspace(0, 1, self.region.shape[0])]

        return self.region,self.region_ppm
            
    def getData(self):
        return self.data

    def getDic(self):
        return self.dic
    
    def getPPMscale(self):
        return self.ppm_scale


class Pseudo3D:


    def __init__(self,fname="test.ft",spectype="pipe"):

        if spectype == "pipe":
            self.dic,self.data = ng.pipe.read(fname)
            self.uc_w0 = ng.pipe.make_uc(self.dic, self.data,1)
            self.uc_w1 = ng.pipe.make_uc(self.dic, self.data,2)
            self.ppm_scale_w0 = self.uc_w0.ppm_scale()
            self.ppm_scale_w1 = self.uc_w1.ppm_scale()

    def get_region(self,w0,w1):
        """ 
            ppm values
            w1 = [x1,x2]
            w0 = [y1,y2]
        """
        y1,y2 = w0
        x1,x2 = w1
        if y1<y2:
            y1,y2=y2,y1
        if x1<x2:
            x1,x2=x2,x1

        #print start_ppm,end_ppm
        start_pts_w0,end_pts_w0 = self.uc_w0(y1,"ppm"),self.uc_w0(y2,"ppm")
        start_pts_w1,end_pts_w1 = self.uc_w1(x1,"ppm"),self.uc_w1(x2,"ppm")
        #print start_pts,end_pts
        self.region = self.data[:,start_pts_w0:end_pts_w0+1,
                                  start_pts_w1:end_pts_w1+1]

        self.region_ppm = [self.ppm_scale_w0[start_pts_w0:end_pts_w0+1],
                           self.ppm_scale_w1[start_pts_w1:end_pts_w1+1]]

        return self.region,self.region_ppm

    def get_data(self):
        return self.data


def normalise(x):
    return (x-np.min(x))/(np.max(x)-np.min(x))

def get_region(data,ppm_scale,start,end):
    region = data[:,start:end+1]
    region_ppm = ppm_scale[start:end+1]
    return region,region_ppm
    
def integrate(region):
    """ Data is 1d array """
    area = region.sum()
    return area

    regions,region_ppm = get_region(data,ppm_scale,start,end)

    """ Integration """ 
    areas = np.array([integrate(i) for i in regions])
    """ normalise integrals """
    I0 = areas[0]
    areas = areas/I0
    """ Start params """    
    I0 = areas[0]
    D = 1e-10
    """ Gradient strengths """
    #procpar = ng.bruker.read_procpar("procpar")
    #gs = np.array(procpar["gzlvl1"]["values"],dtype="float")
    #g2s = np.square(gs*0.00179)
    grad_max = 45. # Gcm-1
    gs = np.array([.1,.2,.3,.4,.5,.6,.7,.8,.9])
    gs = gs*grad_max
    g2s = np.square(gs)
    """ Fitting """
    popt, pcov = curve_fit(func, g2s[:-1], areas[:-1],[I0,D])
    #class_fit = FitDiffusion()
    popt, pcov = curve_fit(func, g2s, areas,[I0,D])
    result = "Fitting params\n"+r"I$_0$ = %8.3f"%popt[0]+"\n"+r"D = %8.3e $cm^2s^{-1}$"%popt[1] #+"\n"+r"%8.3e $m^2s^{-1}$"%(popt[1]/10000.)
    
    """ Plotting fits """
    x = np.linspace(g2s.min(),g2s.max())
    #colormap = plt.cm.autumn
    colormap = plt.cm.winter_r
    plt.rcParams['axes.color_cycle'] = [colormap(k) for k in np.linspace(0, 1, regions.shape[0])]
    fig = plt.figure(figsize=(16,8))
    ax1 = fig.add_subplot(121) 
    ax1.plot(x,func(x,*popt),"k--") 
    ax1.plot(g2s,areas,"o")
    ax1.text(g2s.max()*.7,0.97,result)
    ax1.set_ylabel(r"I/I$_{0}$")
    ax1.set_xlabel(r"G$^2$ $(G^2cm^{-2})$")
    """ Plotting spectra """
    ax2 = fig.add_subplot(122)
    #ax2 = fig.add_subplot(122,projection='3d')
    ax2.set_xlim(max(region_ppm),min(region_ppm))
    [ax2.plot(region_ppm,region,label="%d"%g2) for g2,region in zip(g2s,normalise(regions))]
    #[ax2.plot(region_ppm,np.zeros(region_ppm.shape[0])+grad,zs=region) for grad,region in zip(grads,normalise(regions))]
    ax2.set_xlabel("ppm")
    ax2.set_ylabel("normalised intensity")
    #ax2.set_ylabel(r"G$^2$")
    #ax2.set_zlabel("normalised intensity")
    ax2.legend(title="Gradient strength")
    plt.tight_layout()
    figname="result%.1fto%.1f.pdf"%(start_ppm,end_ppm)
    #figname="result%.1fto%.1f_wo_lastpoint.pdf"%(start_ppm,end_ppm)
    plt.savefig(figname)
