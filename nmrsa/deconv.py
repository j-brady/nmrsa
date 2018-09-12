from lmfit.models import (GaussianModel, ExponentialModel,
        LorentzianModel, VoigtModel, PseudoVoigtModel)
from lmfit import Parameters, Model
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from nmrsa.global_fitting import residual

np.random.seed(1987)

def test_exp_dec(Z,d,npoints=6):
    t = np.arange(npoints)
    ints = []
    for t_i in t:
        ints.append(Z * np.exp(-d*t_i))
    return np.array(ints)

def gauss2D(xy,amplitude,sigma_x,sigma_y,center_x,center_y):
    """ p is lmfit Parameter object """
    x, y = xy
    g = amplitude * np.exp(-((x-center_x)**2./2.*sigma_x**2. + (y-center_y)**2./2.*sigma_y**2.))
    # Convert to 1D for fitting
    return g.ravel()

def update_params(params,values):
    for k,v in zip(params,values):
        params[k].value = v

def make_models(model,peaks):
    """ Make composite models for multiple peaks 
        
        Arguments:
            -- models
            -- peaks: dict of {prefix:"peak_name_", model_params:()}
    """
    if len(peaks)<2:
        mod = Model(model,prefix=peaks[0][0])
        p_guess = mod.make_params()
        values = peaks[0][1:]
        update_params(p_guess,values)

    elif len(peaks)>1:

        mod = Model(model,prefix=peaks[0][0])
        
        values = peaks[0][1:]
        print(values)
        for i in peaks[1:]:
            mod += Model(model,prefix=i[0])
            values.extend(i[1:])

        p_guess = mod.make_params()
        
        update_params(p_guess,values)

        return mod, p_guess

def make_test_peaks(model,xy,param_list,sigma=0.1):
    mu = 0
    data = model(xy,*param_list[0])
    data_noise = data + data*np.random.normal(mu,sigma,data.shape)
    for i in param_list[1:]:
        d = model(xy,*i)
        data+=d
        data_noise+=d + d*np.random.normal(mu,sigma,d.shape)

    return data, data_noise

def fix_params(params,to_fix):
    """ Set parameters to fix """ 
    for k in params:
        for p in to_fix:
            if p in k:
                params[k].vary = False

def get_params(params,name):
    ps = []
    ps_err = []
    for k in params:
        if name in k:
            ps.append(params[k].value)
            ps_err.append(params[k].stderr)
    return ps, ps_err


if __name__ == "__main__":
    
    #mod = GaussianModel(prefix="g1_")
    #mod2 = GaussianModel(prefix="g2_")
    #x = np.linspace(-5,5)
    #sim_y = mod.func(x,center=-2.,amplitude=10,sigma=1)
    #sim_y2 = mod2.func(x,center=2.,amplitude=5,sigma=0.75)
    ##sim_y = sim_y + sim_y*np.random.normal(mu,sigma,x.shape)
    ##sim_y2 = sim_y2 + sim_y2*np.random.normal(mu,sigma,x.shape)
    ## add noise
    #conv = sim_y+sim_y2 + (sim_y+sim_y2)*np.random.normal(mu,sigma,x.shape)
    #print(sim_y)
    ## now fit
    #mod = mod + mod2
    #params = mod.make_params()
    #out = mod.fit(conv,params,x=x)    
    #print(out.fit_report(min_correl=0.5))
    #plt.plot(x, out.best_fit, 'r-')
    ##plt.plot(x,sim_y)
    ##plt.plot(x,sim_y2)
    ##plt.plot(x,sim_y+sim_y2)
    #plt.plot(x,conv)

    #plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x = np.linspace(-5,5)
    y = np.linspace(-5,5)
    X, Y = np.meshgrid(x,y) 

    test_peaks = [[10,1,1,1,1],
                  [7,1,1,-2,-2],
                  [7,1,1,3,3],
                  [7,1,1,2,-2],
                  [4,1,1,-3,3],
                  [6,1,1,3,-3],
                  [12,1,1,-2,2],
                  [12,1,1,0,0]]

    test_data, test_data_noise = make_test_peaks(gauss2D,(X,Y),test_peaks,sigma=0.3)
    exp_test_data = test_exp_dec(test_data_noise,.5,6)
    print(exp_test_data.shape)

    peaks = [["g1_",10,1,1,1,1],
             ["g2_",7,1,1,-2,-2],
             ["g3_",7,2,1,3,3],
             ["g4_",7,1,1,2,-2],
             ["g5_",7,1,1,-3,3],
             ["g6_",7,1,1,3,-3],
             ["g7_",7,1,1,-2,2],
             ["g8_",7,1,1,0,0]]

    mod, p_guess = make_models(gauss2D,peaks)
    print(p_guess)

    out = mod.fit(test_data_noise,p_guess,xy=(X,Y))
    to_fix = ["sigma","center"]
    fix_params(out.params,to_fix)
    for k in out.params:
        print(out.params[k].vary)
    amps = []
    amps_std = []
    for d_i in exp_test_data:
        o_i = mod.fit(d_i,out.params,xy=(X,Y))
        amp, amp_err = get_params(o_i.params,"amplitude")
        amps.append(np.array(amp))
        amps_std.append(np.array(amp_err))
        print(o_i.fit_report(min_correl=0.5))
    amps = np.vstack(amps)
    amps_std = np.vstack(amps_std)

    #print(out.fit_report(min_correl=0.5))
    ##ax.plot_wireframe(X, Y, Z_noise.reshape(50,50),color="r")
    #ax.plot_wireframe(X, Y, Z_conv.reshape(50,50),color="r")
    ax.plot_wireframe(X, Y, test_data_noise.reshape(50,50),color="r")
    colors = ["r","g","b","orange","k","c"]
    #for i,c in zip(exp_test_data,colors):
    #    ax.plot_wireframe(X, Y, i.reshape(50,50),color=c)

    ax.plot_wireframe(X,Y,out.best_fit.reshape(50,50),color="b")
    ##ax.contour(X,Y,Z)
    ax.set_ylabel("y")
    ax.set_xlabel("x")
    plt.show()
    t = np.arange(6)

    for i in range(len(exp_test_data)):
        mod = ExponentialModel()
        pars = mod.guess(amps[:,i],x=t)
        fit = mod.fit(amps[:,i],pars,x=t)

        plt.errorbar(t,amps[:,i],yerr=amps_std[:,i],fmt="ro")
        plt.plot(t,fit.best_fit,"k--")

        plt.show()
