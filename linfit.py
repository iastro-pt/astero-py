import math
import numpy as np
from scipy import optimize
from pylab import *

def linfit(x, y, err=0.0, use_err=False, full_output=False, plot=False):
    """
    Finds the least-squares fit to y = a + bx
     
    inputs:
       x            sequence of independent variable data
       y            sequence of dependent variable data
       err          sequence of uncertainties in dependent variable (by default,
                     if it is not provided, a standard unweigthed fit is made)
       full_output  (optional, boolean) return dictionary with all results and stats
       
    outputs (not using full_output):
       a      y intercept of regression line
       b      slope of regression line
       
    full output = True:
       also returns dictionary containing statistics of fit
       key    value
       ---    -----------------------
       a      a in: y = a + bx
       b      a in: y = a + bx
       ap     a' in: x = a' + b'y
       bp     b' in: x = a' + b'y
       r2     correlation coeffecient
       var_x  variance of x (sigma**2)
       var_y  variance of y (sigma**2)
       cov    covariance
       SEa    standard error for a
       SEb    standard error for b
    """

    # notation as in
    # http://mathworld.wolfram.com/CorrelationCoefficient.html
    # http://mathworld.wolfram.com/LeastSquaresFitting.html

    x=np.asarray(x)
    y=np.asarray(y)

    n = len(x)
    assert n == len(y)

    if (use_err): aa, bb, berr, corr_coef = wls(x, y, err)

    mean_x = np.sum(x.astype(np.float))/n
    mean_y = np.sum(y.astype(np.float))/n

    SSxx = np.sum( (x-mean_x)**2 )
    SSyy = np.sum( (y-mean_y)**2 )
    SSxy = np.sum( x*y ) - n*mean_x*mean_y

    # y = a + b x
    # x = ap + bp y
    b = SSxy/SSxx
    bp = SSxy/SSyy

    a = mean_y - b*mean_x

    if not full_output:
        return a, b

    ap = mean_x - bp*mean_y

    s2 = (SSyy - b*SSxy)/(n-2)
    s = math.sqrt(s2)

    SEa = s*math.sqrt( 1.0/n + mean_x**2/SSxx )
    if(use_err): 
        SEb = berr
        r2 = corr_coef
    else: 
        SEb = s/math.sqrt(SSxx)
        r2 = b*bp

    stats = dict(
        r2 = r2,             # correlation coefficient
        var_x = SSxx/n, # variance of x (sigma**2)
        var_y = SSyy/n, # variance of y (sigma**2)
        cov = SSxy/n,   # covariance
        SEa = SEa,      # standard error for a
        SEb = SEb,      # standard error for b
        a = a,          # a in: y = a + bx
        b = b,          # b in: y = a + bx
        ap = ap,        # a' in: x = a' + b'y
        bp = bp,        # b' in: x = a' + b'y
     )
    
    # plot the result?
    if(plot and use_err): plotlin(x, y, aa, bb)
    if(plot and not use_err): plotlin(x, y, a, b)

    if(use_err): return aa, bb, stats
    else: return a, b, stats


def wls(x, y, err):
    """
    Performs a weigthed least-squares fit to the data
    """

    fit = optimize.leastsq

    # fit functions
    fitfunc = lambda p,x: p[0] + p[1]*x
    errfunc = lambda p,x,y,errors : (y-fitfunc(p,x))/errors

    p0 = [1.0,-1.0]
    out = fit(errfunc, p0, args=(x, y, err), full_output=1)
    
    ##
    one_err = np.ones(len(x))
    out_inverted = fit(errfunc, p0, args=(y, x, one_err))
    bp = out_inverted[0][1]
    ##

    best_params = out[0]    # fitted parameters
    cov = out[1]            # covariance matrix

    a = best_params[0]
    b = best_params[1]

    berr = sqrt(diag(cov))[0]
    r2 = b*bp
#    Sw = np.sum(w)
#    Swy = np.sum(w*y)
#    Swx = np.sum(w*x)

#    Swx2 = np.sum(w*x*x)
#    Swxy = np.sum(w*x*y)
#    Swxwx = np.sum((w*x)**2)

#    # y = a + b x
#    a = (Swy*Swx2 - Swx*Swxy) / (Sw*Swx2 - Swxwx)
#    b = (Sw*Swxy - Swx*Swy) / (Sw*Swx2 - Swxwx)
    
    return a, b, berr, r2 


def plotlin(x, y, a, b):
    """
    Plots the best fit line
    """

    xfit = np.linspace(min(x), max(x), 100)
    yfit = a + b*xfit

    plot(x, y, 'ko', xfit, yfit, 'b-')
    show()

    return
    






