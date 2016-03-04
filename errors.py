""" Functions for estimating errors in fitting nmr data """
import numpy as np
#import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gen_data(array,sigma=0.1):
    """ 
        Generate array with noise added using np.random.normal
    """
    array_len = len(array)
    t = np.linspace(0,1,array_len)
    noise = np.random.normal(0,sigma,array_len)
    data_w_noise = array+noise
    return t,data_w_noise


def monte_carlo(func,x,popt,raw_y_vals,iterations=200):
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
    std = np.sqrt(np.square(raw_y_vals-fit_y_vals)/len(x))
    fit_params = []

    for _ in range(iterations):
        t,sim_y_vals = gen_data(fit_y_vals,sigma=std)
        popt,pcov = curve_fit(func,x,sim_y_vals,popt)
        #plt.plot(x,sim_y_vals,"o")
        fit_params.append(popt)
                                                
    return np.array(fit_params)
