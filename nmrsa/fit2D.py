import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import pandas as pd
import nmrglue as ng
from nmrglue.analysis import linesh,peakpick
from lmfit import Model
from lmfit.models import LorentzianModel, GaussianModel, VoigtModel

from scipy.interpolate import griddata

from ng import Pseudo3D


model = {"G":GaussianModel(),
         "L":LorentzianModel(),
         "V":VoigtModel()}


def get_peak_bounds(peak,F1_width=0.20,F2_width=0.80):
    f1 = [peak.F1 + F1_width/2. , peak.F1 - F1_width/2]
    f2 = [peak.F2 + F2_width/2. , peak.F2 - F2_width/2]
    return f1,f2

def peak_mask(shape,center,radius):
    f2,f1 = np.ogrid[:shape[0],:shape[1]]
    cf2,cf1 = center



if __name__ == "__main__":
    spectrum = Pseudo3D("spectra/pseudo3d/buffer_1432.ft") 
    peaks= pd.read_table("spectra/pseudo3d/buffer_peaks1432.txt",comment="#",names=["F1","F2","Assignment"],delim_whitespace=True)
    peak = peaks.ix[1]
    print(peaks)
    f1,f2 = get_peak_bounds(peak)
    region,region_ppm = spectrum.get_region(f2,f1)

    #region_ppm_f1 = np.vstack([region_ppm[0] for _ in region[0,:,0]])
    #region_ppm_f2 = np.vstack([region_ppm[1] for _ in region[0,0,:]])

    
    #X = np.array(region_ppm[0])
    #print(X.shape)
    #Y = np.array(region_ppm[1])
    #print(Y.shape)
    #Z = region[0].T
    
    #X,Y = np.meshgrid(*region_ppm)
    #print(X.shape,Y.shape,Z.shape)
    #print(Z.ndim)
    #peaks = peakpick.pick(region, 1e5)
    #print("peaks")
    #print(peaks)
    #params = [[(x, x_lw),(y, y_lw)] for x, x_lw, y, y_lw in zip(peaks['X_AXIS'], peaks['X_LW'], peaks['Y_AXIS'], peaks['Y_LW'])]
    #print(params)
    #lineshapes = ["l","l"]
    #params_best, amp_best, iers = linesh.fit_NDregion(Z,lineshapes,params,peaks["VOL"])
    #print(params_best)

    # simulate the spectrum
    #sdata = linesh.sim_NDregion(Z.shape, lineshapes, params_best, amp_best)
    #print(sdata.shape)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(X,Y,Z,lw="2")
    
    #ax.plot_wireframe(X,Y,sdata,color="r")
#    ax.scatter(peaks['X_AXIS'],peaks['Y_AXIS'])
    plt.show()

    
    #pars = model.guess(region,x=region_ppm)

    #result = model.fit(region,x=region_ppm,amplitude=region.max())
    #print result.fit_report()

    #plt.plot(region_ppm, region,         'bo')
    #plt.plot(region_ppm, result.init_fit, 'k--')
    #plt.plot(region_ppm, result.best_fit, 'r-')
    #plt.show()
