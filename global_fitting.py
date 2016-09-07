import random

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq, curve_fit
from lmfit import minimize


""" For debugging """
#np.random.seed(1987)

def model(x,params):
    return None

def gaussian(x,mu,sigma):
    a = 1./(sigma*np.sqrt(2*np.pi))
    return a*np.exp(-np.square(x-mu)/(2.*sigma**2.))

def residual(params,y,x,function,global_fit=False,sigma=None,lmfit=False):
    
    """ Function for calculating residuals for leastsq fit.
        
        Args:

            params (list): your list of parameters in same order as arguments for your function/fitting model e.g. [mu,sigma]
            y: your y values. This can be a list of y datasets if global_fit is True
            x: your x values. This can be a list of x datasets if global_fit is True
            function: your fitting model. Called like so, function(x,*params)
            global_fit: If set to True you can supply the function with a list of x's,y's and functions for global fitting.

        Returns:

            y - model(x,*params)

        Example:

            for a global fit:
                x = [x1,x2,x3]
                y = [y1,y2,y3]
                f = [f1,f2,f3]
                global_fit=True
                pO = [list of starting params] # same order as for *params in fitting function
                result, pcov, infodict, errmsg, success = leastsq(residual,p0,args=(y,x,f,global_fit),full_output=True)


    """

    if sigma is None:

        sigma = [1. for _ in range(len(x))]

    else:

        sigma = sigma


    if global_fit:

        errs = []

        for y_i, x_i, f_i, s_i in zip(y,x,function,sigma):
            if lmfit:
                err = (y_i - f_i(x_i, params)) / s_i
		errs.extend(err)

            else:
                err = (y_i - f_i(x_i, *params)) / s_i
		errs.extend(err)

        errs = np.array(errs).flatten()
        return errs
# This commented code does not work with pandas data Series ...
#        err = [y_i - f_i(x, *params) for y_i, x_i, f_i in zip(y,x,function)]
#        err = np.array(err).flatten()
#        return err

    else:
        if lmfit:
            model = function(x,params)
        else:
            model = function(x,*params)
        return (y-model)/sigma

def resample(x,y,yerr=None,max_replacement=100):
    """ Resample dataset with replacement using numpy.random.choice for bootstrap fit.
        
        Args:
            x: set of x values
            y: set of y values
            yerr: set of y errors
            max_replacement: maximum percentage of data points that can be replaced. Important for small datasets.

        Returns:
            if yerrs:
                resample_x, resample_y
            else:
                resample_x, resample_y, resample_yerr
                
        Example:
            
    """ 

    if len(x) == len(y):

        length = len(x)
        #inds = sorted([random.choice(range(length)) for _ in range(length) ])
        #print inds
        #resample_x = x[inds]
        #resample_y = y[inds]

        # Code to use numpy.random.choice
        data_indices = np.arange(length)
        choice = np.random.choice(data_indices,int(length*max_replacement/100.))
        replace = np.random.choice(data_indices,int(length*max_replacement/100.))
        data_indices[choice] = data_indices[replace]
        inds = np.sort(data_indices)

        #inds = np.sort(np.random.choice(data_indices,length))
        print inds
        resample_x = x[inds]
        resample_y = y[inds]
        if yerr is None:
            return resample_x, resample_y

        elif len(yerr) == len(y):
            resample_yerr = yerr[inds]
            return resample_x, resample_y, resample_yerr

        else:
            raise TypeError("yerr should be a list of length y containing the errors of y.")
    else:
        print("len(x) != len(y)!")
        raise TypeError("len(x) should be the same as len(y).")


def bootstrap(params, y, x, function, global_fit=False, lmfit=False, yerr=None, iterations=100, montecarloX=False, mc_kwargs={"std_x":0.1,"iterations":200}, **kwargs):
    if montecarloX:
        from errors import monte_carlox
    else:
        pass
    results = []
    
    for _ in range(iterations):

        if global_fit:
            if yerr is None:    
                data = np.array([resample(x_i,y_i,**kwargs) for x_i, y_i in zip(x,y)])
                x_rs = data[:,0]
                y_rs = data[:,1]
                print x_rs, y_rs
            else:
                data = np.array([resample(x_i,y_i,yerr_i,**kwargs) for x_i, y_i, yerr_i in zip(x,y,yerr)])
                x_rs = data[:,0]
                y_rs = data[:,1]
                yerr_rs = data[:,2]
                print x_rs, y_rs, yerr_rs
             
        else:
            if yerr is None:
                x_rs, y_rs = resample(x,y,**kwargs)
                print x_rs,y_rs
                yerr_rs = None
            else:
                x_rs, y_rs, yerr_rs = resample(x,y,yerr,**kwargs)
                print x_rs,y_rs,yerr_rs
             
#        print data
        #if yerrs is None:
        #    yerrs = [1. for i in y_rs]
        #else
        if montecarloX:
            """ Do montecarlo simulation of x for every y bootstrap """
            std_x = mc_kwargs["std_x"]
            mc_iterations = mc_kwargs["iterations"]
            result = monte_carlox(function, x_rs, params, y_rs, std_x, yerr_rs, global_fit, lmfit, iterations=mc_iterations)
            results.append(result)
        else:
            if lmfit:     
                result = minimize(residual,params,args=(y_rs, x_rs, function, global_fit, yerr_rs, lmfit))
            else:
                result, pcov, infodict, errmsg, success = leastsq(residual, params, args=(y_rs, x_rs, function, global_fit, yerr_rs),full_output=True)
            print result
            results.append(result)
#    print results
    if lmfit:
        return np.array(results).ravel()
    else:
        return np.array(results)

if __name__ == "__main__":

    mu = 10.
    sigma = 1.

    d1 = np.random.normal(mu, sigma, 100)
    d2 = np.random.normal(mu, sigma, 100)
    d3 = np.random.normal(mu, sigma, 100)
    d4 = np.random.normal(mu, sigma, 100)
    d5 = np.random.normal(mu, sigma, 100)
    d6 = np.random.normal(mu, sigma, 100)
    data = np.array([d1,d2,d3,d4,d5,d6])

    hist, bin_edges = np.histogram(d1,20,density=True)
    p0 = [mu,sigma]
    x = bin_edges[:-1] 
    y = hist
    result, pcov, infodict, errmsg, success = leastsq(residual,p0,args=(y,x,gaussian),full_output=True)
    print result
    bs = bootstrap(p0, y, x, gaussian, global_fit=False,iterations=1000)
    av_mu = np.mean(bs[:,0])
    std_mu = np.std(bs[:,0])
    print av_mu,std_mu
    

    """ Global fit """
#    bins = np.min(data),np.max(data)
#    print bins
#    data = np.array([np.histogram(d,20,density=True) for d in data])
#    x = data[:,1]
#    x = [i[:-1] for i in x]
#    y = data[:,0]
#    funcs = [gaussian for _ in x]
#    global_fit=True
#
#    result, pcov, infodict, errmsg, success = leastsq(residual,p0,args=(y,x,funcs,global_fit),full_output=True)
#    print result
#    bs = bootstrap(p0, y, x, funcs, global_fit=True,iterations=100)
#    av_mu = np.mean(bs[:,0])
#    std_mu = np.std(bs[:,0])
#    av_sigma = np.mean(bs[:,1])
#    std_sigma = np.std(bs[:,1])
#    print result
#    print av_mu,std_mu,av_sigma,std_sigma
#
#    c = ["r","g","b","y","orange","pink"]
#    #plt.clf()
#    for i,j,k in zip(x,y,c):
#        plt.bar(i,j,width=0.75,color=k,alpha=0.5)#plt.plot(bins,gaussian(bins,mu,sigma))
#    sim_x = np.linspace(min(x[0]),max(x[0]),100)
#    plt.plot(sim_x,gaussian(sim_x,av_mu,av_sigma),label=r"$ \mu = %.3f \pm %.3f, \sigma = %.3f \pm %.3f$" % (av_mu,std_mu,av_sigma,std_sigma) )
#    #plt.plot(x,y)
#    plt.legend(frameon = False)
#    plt.title("Test global bootstrap fit for Gaussian data.")
#    #print result
#    plt.show()
