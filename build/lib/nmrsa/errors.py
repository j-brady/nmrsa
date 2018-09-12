""" Functions for estimating errors in fitting nmr data """
import sys
import numpy as np
#import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from .global_fitting import residual, resample
from lmfit import minimize

#np.random.seed(1987)
def multiplicative(sigmas):
    """ list of arrays of sigma values """
    return np.sqrt(sum([np.square(s) for s in sigmas])) 

def logarithmic(xs):
    """ List of tuples containing [(array of x,array of sigmas )]"""
    l = []
    for x,s in xs:
        sl = 0.434(s/x)
        l.append(sl)
    return np.array(l)

def gen_data(array,sigma=0.1):
    """ 
        Generate array with noise added using np.random.normal
    """
    print("sigma")
    print(sigma)
    array_len = len(array)
    t = np.linspace(0,1,array_len)
    noise = np.random.normal(0,sigma,array_len)
    data_w_noise = array+noise
    return t,data_w_noise


def redchi(f,x,popt,raw_y_vals,yerrs):
    residuals = (f(x,*popt) - raw_y_vals)/yerrs
    N = len(x)
    n = len(popt)
    return sum(residuals**2)/(N-n)

def monte_carlo(func,x,popt,raw_y_vals,iterations=200,red_chi_sqrs=False,yerrs=None):
    """ 
        Monte Carlo error estimation. Standard deviation of raw_y_vals from y_vals generated using 
        optimised fitting parameters is used to add noise to the simulated y_vals which are then fitted
        using curve_fit.

        parameters:
            func         -- function used for fitting with curve_fit
            x            -- x values (np.array)
            popt         -- optimised parameters from curve_fit (list)
            raw_y_vals   -- raw y data values (np.array)  
            iterations   -- number of iterations (int)
        
        returns:
            np.array of simulated fitting parameters --> np.array([popt,popt,etc])

    """

    fit_y_vals = func(x,*popt)
    #plt.plot(x,fit_y_vals)
    std = sum(np.abs(raw_y_vals-fit_y_vals))/len(x)
    fit_params = []

    if red_chi_sqrs and yerrs is not None:
        red_chis = []
    if red_chi_sqrs and yerrs is None:
        print("You need to have the errors in y to calculate chis", sys.exc_info()[0])
        raise

    for _ in range(iterations):
        t,sim_y_vals = gen_data(fit_y_vals,sigma=std)
        popt,pcov = curve_fit(func,x,sim_y_vals,popt)
        #plt.plot(x,sim_y_vals,"o")
        fit_params.append(popt)

        # Calucalate chi_squares
        if red_chi_sqrs:
            red_chis.append(redchi(func,x,popt,raw_y_vals,yerrs))


    if red_chi_sqrs:                                            
        return np.array(fit_params), np.array(red_chis)

    else:
        return np.array(fit_params)

def monte_carlox(func,x,p,y,std_x,y_err=None,global_fit=False,lmfit=True,iterations=200):
    """ 
        Monte Carlo error estimation in x values. Standard deviation of raw_y_vals from y_vals generated using 
        optimised fitting parameters is used to add noise to the simulated y_vals which are then fitted
        using curve_fit.

        parameters:
            func         -- function used for fitting with curve_fit. If lmfit is used then function should be lmfit compatible.
            x            -- x values (np.array)
            p            -- optimised parameters from curve_fit (list) or lmfit parameter object
            y            -- raw y data values (np.array)  
            std_x        -- estimated standard error of x values to use in simulation (e.g. 0.1 for 10%)
            y_err        -- errors in y points 
            lmfit        -- if True lmfit minimize module is used for fitting, else curvefit used
            iterations   -- number of iterations (int)
        
        returns:
            if lmfit:
                numpy array of MinimizerResult objects
            else:
                numpy array of simulated fitting parameters --> np.array([popt,popt,etc])

    """
    if y_err is None:
        y_err = None
    elif len(y_err) == len(y):
        print("Using y_errs to weight fit")
        y_err = y_err
    else:
        print("len(y_err) != len(y)")
        raise TypeError("y_err should be list (or np.array) of errors in y with length y")

    fit_results = []
    
    for _ in range(iterations):
        dist = np.random.normal(0.0,std_x,len(x))
        sim_x = x+(x*dist)

        if lmfit:
            result = minimize(residual,p,args=(y,sim_x,func,global_fit,y_err,lmfit))
            fit_results.append(result)
        else:
            popt,pcov = curve_fit(func,sim_x,y,popt)
            fit_results.append(popt)
                                                
    return np.array(fit_results)

def tmean_std(array,percentiles):
    """ Calculates the trimmed mean and the std of the trimmed data
    
        Arguments:
        array -- numpy array of data
        percentiles -- tuple of lower and upper percentiles to use for tmean

        Returns:
        tmean -- trimmed mean
        std -- standard deviation of trimmed region
        indices -- indices of data used for tmean calculation

    """
    # percentiles
    lp = np.percentile(array,percentiles[0])
    up = np.percentile(array,percentiles[1])
    mask = np.ma.masked_inside(array, lp, up)
    tarray = array[mask.mask]
    tmean = tarray.mean()
    std = tarray.std()
    mask = mask.mask
    return tmean, std, mask

def bootstrapYmontecarloX(n_bs=100,n_mc=100):
    """ 
        Generate n data sets by randomly replacing a subset of points.
        
    """
