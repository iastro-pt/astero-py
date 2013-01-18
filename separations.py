import numpy as np
import uncertainties as un


###############################################################################
# Large separation
#  Dl(n) =   nu     -  nu
#              n,l       n-1,l
###############################################################################
def large_sep ( mx , mxErr ):
    """ 
    Given a frequency matrix and errors calculates the large 
    separations for each degree l and propagates errors.
    """

    (Nn, Nl) = np.shape(mx)
    Dl    = np.zeros((Nn-1, Nl))
    DlErr = np.zeros((Nn-1, Nl))

    for l in range(Nl):
        for n in range(1,Nn):
            if (mx[n,l] and mx[n-1,l]):
                a = un.ufloat( (mx[n,l], mxErr[n,l]) )
                b = un.ufloat( (mx[n-1,l], mxErr[n-1,l]) )
                result = a-b
                Dl[n-1,l]    = un.nominal_value(result)
                DlErr[n-1,l] = un.std_dev(result)	 
    Dl[Dl<0.] = 0.
    DlErr[Dl<0.] = 0.
    return Dl, DlErr

###############################################################################
# Small separation 02 (scaled)
#  d02(n) =    nu     -  nu
#                n,0       n-1,2
###############################################################################
def small_sep02 ( mx , mxErr, lag=0):
	""" 
	Given a frequency matrix and errors calculates the (scaled) small 
	separation d02 and propagates errors.
	
	Notes
	-----
	The parameter lag is the difference between the radial orders of
	the first modes of l=0 and l=2.  
	"""

	(Nn, Nl) = np.shape(mx)
	d02    = np.zeros((1,Nn))
	d02Err = np.zeros((1,Nn))
	d02.fill(None)      # values that can't be calculated are NaN

	
	for n in range(1-lag,Nn):
		if (mx[n,0] != 0. and mx[n-1+lag,2] != 0.):
		    a = un.ufloat( (mx[n,0], mxErr[n,0]) )
		    b = un.ufloat( (mx[n-1+lag,2], mxErr[n-1+lag,2]) )
		    result = (a-b) / 3.	
		    d02[0,n-1+lag]    = un.nominal_value(result)
		    d02Err[0,n-1+lag] = un.std_dev(result)
		
	return d02, d02Err


###############################################################################
# Small separation 13 (scaled)
#  d13(n) =    nu     -  nu
#                n,1       n-1,3
###############################################################################
def small_sep13 ( mx , mxErr, lag=0):
	""" 
	Given a frequency matrix and errors calculates the (scaled) small 
	separation d13 and propagates errors.
	
	Notes
	-----
	The parameter lag is the difference between the radial orders of
	the first modes of l=1 and l=3.  
	"""

	(Nn, Nl) = np.shape(mx)
	d13    = np.zeros((1,Nn))
	d13Err = np.zeros((1,Nn))
	d13.fill(None)      # values that can't be calculated remain NaN

	for n in range(1-lag,Nn):
		if (mx[n,1] != 0. and mx[n-1+lag,3] != 0.):     
        		a = un.ufloat( (mx[n,1], mxErr[n,1]) )
        		b = un.ufloat( (mx[n-1+lag,3], mxErr[n-1+lag,3]) )
        		result = (a-b) / 5.
        		d13[0,n-1+lag]    = un.nominal_value(result)
        		d13Err[0,n-1+lag] = un.std_dev(result)

	return d13, d13Err


###############################################################################
# Separation d01 - 3 point 
# -> Roxburgh (1993, 2009)
#
#  d01(n) =   nu     - ( nu      + nu      )   /
#               n,0    (   n-1,1     n,1 | )  /   2
###############################################################################
def sep01_3(mx, mxErr, lag=0):
	""" 
	Given a frequency matrix and errors calculates the 3-point 
	separation d01 and propagates errors. 
	"""

	(Nn, Nl) = np.shape(mx)
	d01    = np.zeros((1,Nn))
	d01Err = np.zeros((1,Nn))
	d01.fill(None)      # values that can't be calculated are NaN

	for n in range(0,Nn-lag):
		if (mx[n,0] != 0. and mx[n-1+lag,1] != 0. and mx[n+lag,1] != 0. and (n-1+lag) >= 0):
			#print mx[n,0], mx[n-1+lag,1], mx[n+lag,1], n-1+lag, n
			a = un.ufloat( (mx[n,0],mxErr[n,0]) )
			b = un.ufloat( (mx[n-1+lag,1],mxErr[n-1+lag,1]) )
			c = un.ufloat( (mx[n+lag,1],mxErr[n+lag,1]) )
			result = a - (b+c) / 2.
			d01[0,n]    = un.nominal_value(result)
			#print d01[0,n]
			d01Err[0,n] = un.std_dev(result) 

	return d01, d01Err


###############################################################################
# Separation d10 - 3 point 
# -> Roxburgh (1993) , Roxburgh & Vorontsov (2003)
#
#  d10(n) = | nu    + nu      |  /     - nu   
#           |   n,0     n+1,0 | /  2       n,1
###############################################################################
def sep10_3(mx, mxErr, lag=0):
	""" 
	Given a frequency matrix and errors calculates the 3-point 
	separation d10 and propagates errors. 
	"""

	(Nn, Nl) = np.shape(mx)
	d10    = np.zeros((1,Nn-1))
	d10Err = np.zeros((1,Nn-1))
	d10.fill(None)      # values that can't be calculated are NaN

	for n in range(0,Nn-1):
		if (mx[n+lag,1] != 0. and mx[n,0] != 0. and mx[n+1,0] != 0.):      
        		a = un.ufloat( (mx[n+lag,1], mxErr[n+lag,1]) )
        		b = un.ufloat( (mx[n,0], mxErr[n,0]) )
        		c = un.ufloat( (mx[n+1,0], mxErr[n+1,0]) )
        		result = -a + (b+c) / 2.
        		d10[0,n]    = un.nominal_value(result)
        		d10Err[0,n] = un.std_dev(result)

	return d10, d10Err



################################################################################
## Separation d01 - 5 point 
## -> Roxburgh & Vorontsov (2003)
##
## d01(n) = 1/8 * |nu      - 2nu      + 6nu    - 4nu    + nu      |
##                |  n-1,0      n-1,1      n,0      n,1     n+1,0 |
################################################################################
#def sep01_5 ( mx ):
#	""" given a frequency matrix 'mx', calculates the 5-point separation d01 """

#	(Nn, Nl) = np.shape(mx)
#	d01 = np.zeros((1,Nn-1))
#	for n in range(1,Nn-1):
#		d01[0,n-1]  = mx[n-1,0] - 4*mx[n-1,1] + 6*mx[n,0] - 4*mx[n,1] + mx[n+1,0]
#		d01[0,n-1]  = (1./8.) * d01[0,n-1]

#	return d01



#############################################
##
## Separation d10 - 5 point 
## Roxburgh&Vorontsov (2003)
##
## d10(n)= -1/8 * nu      - 4nu    + 6nu    - 4nu      + nu
##                  n-1,1      n,0      n,1      n+1,0     n+1,1
##
#############################################
#def sep10_5 ( mx ):
#	""" given a frequency matrix 'mx', calculates the 5-point separation d10 """

#	(Nn, Nl) = np.shape(mx)
#	d10 = np.zeros((1,Nn-1))
#	for n in range(1,Nn-1):
#		d10[0,n-1]  = mx[n-1,1] - 4*mx[n,0] + 6*mx[n,1] - 4*mx[n+1,0] + mx[n+1,1]
#		d10[0,n-1]  = -(1./8.) * d01[0,n-1]

#	return d10



###############################################################################
# Ratio r02 
# -> Roxburgh & Vorontsov (2003), Cunha et al. (2007)
#
# r02 = d02/D1 = | nu    -  nu      |  / | nu    - nu      |
#                |   n,0      n-1,2 | /  |   n,1     n-1,1 |
#
###############################################################################
def ratio02 ( mx , mxErr, lag02=0, lag01=0, dnu_mean=0):
	""" 
	Given a frequency matrix and errors calculates the ratio r02
	and propagates errors.
    """

	(Nn, Nl) = np.shape(mx)
	r02    = np.zeros((1,Nn-1))
	r02Err = np.zeros((1,Nn-1))
	r02.fill(None)      # values that can't be calculated are NaN

	for n in range(0,Nn-1):
		if (mx[n,0] != 0. and mx[n-1+lag02,2] != 0. and mx[n,1] != 0. and mx[n-1,1] != 0.
	    and n-1 >= 0 and n-1+lag02 >= 0):
        		a = un.ufloat( (mx[n,0], mxErr[n,0]) )
        		b = un.ufloat( (mx[n-1+lag02,2], mxErr[n-1+lag02,2]) )
        		c = un.ufloat( (mx[n,1], mxErr[n,1]) )
        		d = un.ufloat( (mx[n-1,1], mxErr[n-1,1]) )
        		result = (a-b) / (6.*(c-d))
        		r02[0,n]    = un.nominal_values(result)
        		r02Err[0,n] = un.std_dev(result)

	return r02, r02Err


###############################################################################
# Ratio r13 
# Roxburgh & Vorontsov (2003)
#
# r13 = d13/D0 = | nu    -  nu      |  / | nu      - nu    |
#                |   n,1      n-1,3 | /  |   n+1,0     n,0 |
###############################################################################
def ratio13 ( mx, mxErr, lag=0, dnu_mean=0):
	""" 
	Given a frequency matrix and errors calculates the ratio r13
	and propagates errors.
    """

	(Nn, Nl) = np.shape(mx)
	r13    = np.zeros((1,Nn-1))
	r13Err = np.zeros((1,Nn-1))
	r13.fill(None)      # values that can't be calculated are NaN

	for n in range(0,Nn-1):
		if (mx[n,1] and mx[n-1+lag,3] and mx[n+1,0] and mx[n,0]
		    and (mx[n+1,0]-mx[n,0] < 1.5*dnu_mean) ):     
        		a = un.ufloat( (mx[n,1], mxErr[n,1]) )
        		b = un.ufloat( (mx[n-1+lag,3], mxErr[n-1+lag,3]) )
        		c = un.ufloat( (mx[n+1,0], mxErr[n+1,0]) )
        		d = un.ufloat( (mx[n,0], mxErr[n,0]) )
        
        		result = (a-b) / (10.*(c-d))
        		r13[0,n]    = un.nominal_values(result)
        		r13Err[0,n] = un.std_dev(result)

	return r13, r13Err


############################################################
#
# Ratio r01 
# Roxburgh & Vorontsov (2003)
#
# r01 = d01/D1 = || nu    - | nu      + nu      |   /    ||  / | nu    - nu      |
#                ||   n,0   |   n-1,1     n,1   |  /   2 || /  |   n,1     n-1,1 |


############################################################
def ratio01 ( mx, mxErr, lag=0 ):
	""" 
	Given a frequency matrix and errors calculates the ratio r01
	and propagates errors.
    """

	(Nn, Nl) = np.shape(mx)
	r01    = np.zeros((1,Nn))
	r01Err = np.zeros((1,Nn))
	r01.fill(None)      # values that can't be calculated are NaN

	for n in range(0,Nn-lag):
		if (mx[n,0] and mx[n-1+lag,1] and mx[n+lag,1] and (n-1+lag) >= 0):
			#print mx[n,0], mx[n-1+lag,1], mx[n+lag,1], mx[n+lag,1], mx[n-1+lag,1]
			a = un.ufloat( (mx[n,0], mxErr[n,0]) )
			b = un.ufloat( (mx[n-1+lag,1], mxErr[n-1+lag,1]) )
			c = un.ufloat( (mx[n+lag,1], mxErr[n+lag,1]) )

			result = (a - (b+c)/2.) / (c-b)
			r01[0,n]    = un.nominal_value(result)
			r01Err[0,n] = un.std_dev(result)

	return r01, r01Err


############################################################
#
# Ratio r10 
# Roxburgh & Vorontsov (2003)
#
# r10 = d10/D0 = || nu    + nu      |  /     - nu    |  / | nu      - nu    |
#                ||   n,0     n+1,0 | /  2       n,1 | /  |   n+1,0     n,0 |

#
############################################################
def ratio10 ( mx, mxErr, lag=0):
	""" 
	Given a frequency matrix and errors calculates the ratio r01
	and propagates errors.
    """

	(Nn, Nl) = np.shape(mx)
	r10    = np.zeros((1,Nn-1))
	r10Err = np.zeros((1,Nn-1))
	r10.fill(None)      # values that can't be calculated are NaN
	
	for n in range(0,Nn-1):
		if (mx[n+lag,1] and mx[n,0] and mx[n+1,0]):
			#print mx[n,0], mx[n+1,0], mx[n+lag,1], mx[n+1,0], mx[n,0]
			a = un.ufloat( (mx[n+lag,1], mxErr[n+lag,1]) )
			b = un.ufloat( (mx[n,0], mxErr[n,0]) )
			c = un.ufloat( (mx[n+1,0], mxErr[n+1,0]) )

			result = (-a + (b+c)/2.) / (c-b)
			r10[0,n]    = un.nominal_value(result)
			r10Err[0,n] = un.std_dev(result)

	return r10, r10Err



############################################################
#
# Second differences D2nu
# Houdek & Gough (2007)
#
# D2l = nu      - 2nu    + nu
#         n-1,l      n,l     n+1,l
#
############################################################
def second_dif( mx, mxErr ):
	""" given a frequency matrix 'mx', calculates the second differences """
	
	(Nn, Nl) = np.shape(mx)
	D2    = np.zeros((Nn-1,Nl))
	D2Err = np.zeros((Nn-1,Nl))

	# second differences for l=(0,1,2,...Nl)
	for l in range(Nl):
		for n in range(1,Nn-1):
			a = un.ufloat( (mx[n-1,l], mxErr[n-1,l]) )
			b = un.ufloat( (mx[n,l], mxErr[n,l]) )
			c = un.ufloat( (mx[n+1,l], mxErr[n+1,l]) )

			result = a - 2*b + c
			D2[n-1,l]    = un.nominal_value(result)
			D2Err[n-1,l] = un.std_dev(result)			

	return D2, D2Err



