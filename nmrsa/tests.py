from numpy import arange, exp, linspace
from numpy.random import normal
from lmfit import Parameters

from errors import monte_carlox
from global_fitting import resample

def strExp(x,p):
    amplitude = p["amplitude"].value
    decay = p["decay"].value
    nu = p["nu"].value
    return amplitude * exp(-1./decay * x**nu)


def test_monte_carlox():
    p = Parameters()
    #           (Name,  Value,  Vary,   Min,  Max,  Expr)
    p.add_many(('amplitude',    100,  True, None, None,  None),
               ('decay'    ,   6e-7,  True,  None, None,  None),
               ('nu'       ,  1.0, True, None, None, None))

    x = linspace(0,500,50)
    y = strExp(x,p)
    y_err = abs(y - (y * normal(0.0,0.1,len(y))))
    lmfit = True
    global_fit = False
    func = strExp
    std_x = 0.1
    return monte_carlox(func,x,p,y,std_x,y_err,global_fit,lmfit,100)


def test_resample():
    p = Parameters()
    #           (Name,  Value,  Vary,   Min,  Max,  Expr)
    p.add_many(('amplitude',    100,  True, None, None,  None),
               ('decay'    ,   6e-7,  True,  None, None,  None),
               ('nu'       ,  1.0, True, None, None, None))
    x = linspace(0,500,50)
    y = strExp(x,p)
    y_err = abs(y - (y * normal(0.0,0.1,len(y))))

    return resample(x,y,yerr=None,max_replacement=100)
