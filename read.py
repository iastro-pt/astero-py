import numpy as np

###############################################################################
# Read model frequencies
# type of file:
#	l   n    nu_model   nu_corrected
#	... ...  ...        ...
#
###############################################################################
def readModel( filename ):
	""" 
	Reads frequencies from 'filename' and returns frequency matrix.
	Type of file:
	  l   n    nu_model   nu_corrected
	  ... ...  ...        ...
	"""

	l,n,fmod,fcorr = np.loadtxt(filename, unpack=True)

	lUnq, indUnq = np.unique(l, return_index=True)
	nmax = max(sum(l==0), sum(l==1), sum(l==2), sum(l==3))
	lmax = int(max(l))

	# build the frequency matrices
	mx_mod = (nmax, lmax+1)
	mx_mod = np.zeros(mx_mod)
	mx_corr = (nmax, lmax+1)
	mx_corr = np.zeros(mx_corr)

	# how many n for each l?
	ind = ind=np.zeros((1,4))[0]
	ind[0] = sum(l==0)
	ind[1] = sum(l==1)
	ind[2] = sum(l==2)
	ind[3] = sum(l==3)
	# vector 'ind' contains the number of radial orders n for each l

	# for each l (each column) populate the matrix with 
	# as many n as there are
	k = 0
	for j in range(lmax+1):
		for i in range(nmax):
			if (i < ind[j]):
				mx_mod[i, j]  = fmod[k]
				mx_corr[i, j] = fcorr[k]
				k += 1
			else:
				mx_mod[i, j]  = 0.
				mx_corr[i, j] = 0.

	# returns both matrices
	return mx_mod, mx_corr, n
		
	
###############################################################################
# Read model frequencies
# type of file:
#	l   n   nu   inertia  beta   splitting
#	... ... ...  ...      ...    ...
#
###############################################################################
def readModel2( filename ):
	""" 
	Reads frequencies from 'filename' and returns frequency matrix.
	Type of file:
      l    n    nu   inertia  beta   splitting
      ...  ...  ...  ...      ...    ...
	"""

	l, n, fmod, inert, beta, split = np.loadtxt(filename, unpack=True)

	lUnq, indUnq = np.unique(l, return_index=True)
	nmax = max(sum(l==0), sum(l==1), sum(l==2), sum(l==3))
	lmax = int(max(l))

	# build the frequency matrices
	mx_mod = (nmax, lmax+1)
	mx_mod = np.zeros(mx_mod)

	# how many n for each l?
	ind = (1, lmax+1);	ind = np.zeros(ind);	ind=ind[0]
	i = 1
	ind[0] = np.argmax(n[0:])
	start = ind[0]
	while (start+1 < len(n)):
		ind[i] = np.argmax(n[start+1:])
		start += ind[i] + 1
		i += 1
	# vector 'ind' contains the number of radial orders n for each l

	# for each l (each column) populate the matrix with 
	# as many n as there are
	k = 0
	for j in range(lmax+1):
		for i in range(nmax):
			if (i <= ind[j]):
				mx_mod[i, j]  = fmod[k]
				k += 1
			else:
				mx_mod[i, j]  = 0.

	# returns both matrices
	return mx_mod, n



###############################################################################
# Read model frequencies
# type of file :
#	l    nu_obs   err_obs   l     nu_model
#	...  ...      ...       ...   ...
#
###############################################################################
def readModel3( filename ):
	""" 
	Reads frequencies from 'filename' and returns frequency matrix
	"""

	lobs, fobs, err, l, fmod = np.loadtxt(filename, unpack=True)

	lUnq, indUnq = np.unique(l, return_index=True)
	nmax = max(sum(l==0), sum(l==1), sum(l==2), sum(l==3))
	lmax = int(max(l))

	# build the frequency matrices
	mx_mod = (nmax, lmax+1)
	mx_mod = np.zeros(mx_mod)

	# how many n for each l?
	ind = ind=np.zeros((1,4))[0]
	ind[0] = sum(l==0)
	ind[1] = sum(l==1)
	ind[2] = sum(l==2)
	ind[3] = sum(l==3)
	# vector 'ind' contains the number of radial orders n for each l

	# for each l (each column) populate the matrix with 
	# as many n as there are

	# l=0
	mx_mod[0:ind[0],0] = fmod[0:ind[0]]
	# l=1	
	mx_mod[0:ind[1],1] = fmod[ind[0]: ind[0]+ind[1]]
	# l=2
	mx_mod[0:ind[2],2] = fmod[ind[0]+ind[1]: ind[0]+ind[1]+ind[2]]
	# l=3
	mx_mod[0:ind[3],3] = fmod[ind[0]+ind[1]+ind[2]: ind[0]+ind[1]+ind[2]+ind[3]]

	# return matrix
	return mx_mod

	
###############################################################################
# Read model frequencies
# type of file:
#	l   n   nu   inertia 
#	... ... ...  ...     
#
###############################################################################
def readModel4( filename ):
	""" 
	Reads frequencies from 'filename' and returns frequency matrix
	"""

	l, n, fmod, inert = np.loadtxt(filename, unpack=True)

	nmax = max(sum(l==0), sum(l==1), sum(l==2), sum(l==3))
	lmax = int(max(l))

	# build the frequency matrices
	mx_mod = (nmax, lmax+1)
	mx_mod = np.zeros(mx_mod)

	# how many n for each l?
	ind = ind=np.zeros((1,4))[0]
	ind[0] = sum(l==0)
	ind[1] = sum(l==1)
	ind[2] = sum(l==2)
	ind[3] = sum(l==3)
	# vector 'ind' contains the number of radial orders n for each l

	# for each l (each column) populate the matrix with 
	# as many n as there are
	k = 0
	for j in range(lmax+1):
		for i in range(nmax):
			if (i < ind[j]):
				mx_mod[i, j]  = fmod[k]
				k += 1
			else:
				mx_mod[i, j]  = 0.

	# returns matrix
	return mx_mod, n



###############################################################################
# Read frequencies from observations
# type of file:
#	l   n   nu   err
#	... ... ...  ...
#
###############################################################################
def readObs( filename ):
	""" 
	Reads frequencies from 'filename' and returns frequency matrices
	"""

	l,n,fobs,err = np.loadtxt(filename, unpack=True)

	lUnq, indUnq = np.unique(l, return_index=True)
	nmax = max(sum(l==0), sum(l==1), sum(l==2), sum(l==3))
	lmax = int(max(l))

	# build the frequency matrix
	mx_obs = (nmax, lmax+1)
	mx_obs = np.zeros(mx_obs)
	# and the error matrix
	mx_err = (nmax, lmax+1)
	mx_err = np.zeros(mx_err)

	# how many n for each l?
	ind = ind=np.zeros((1,4))[0]
	ind[0] = sum(l==0)
	ind[1] = sum(l==1)
	ind[2] = sum(l==2)
	ind[3] = sum(l==3)
	# vector 'ind' contains the number of radial orders n for each l

	# for each l (each column) populate the matrix with 
	# as many n as there are
	k = 0
	for j in range(lmax+1):
		for i in range(nmax):
			if (i < ind[j]):
				mx_obs[i, j]  = fobs[k]
				mx_err[i,j] = err[k]
				k += 1
			else:
				mx_obs[i, j]  = 0.
				mx_err[i,j] = 0.
	# returns both matrices
	return mx_obs, mx_err, n


###############################################################################
# Read frequencies from observations
# type of file:
#	l    nu   err  flag1  flag2
#	...  ...  ...  ...    ...
#
###############################################################################
def readObs2( filename ):
	""" 
	Reads frequencies from 'filename' and returns frequency matrix
	"""

	l,fobs,err,flg1,flg2 = np.loadtxt(filename, unpack=True)

	lUnq, indUnq = np.unique(l, return_index=True)
	nmax = max(sum(l==0), sum(l==1), sum(l==2), sum(l==3))
	lmax = int(max(l))

	# build the frequency matrix
	mx_obs = (nmax, lmax+1)
	mx_obs = np.zeros(mx_obs)
	# and the error matrix
	mx_err = (nmax, lmax+1)
	mx_err = np.zeros(mx_err)

	# how many n for each l?
	ind = ind=np.zeros((1,4))[0]
	ind[0] = sum(l==0)
	ind[1] = sum(l==1)
	ind[2] = sum(l==2)
	ind[3] = sum(l==3) 
	# vector 'ind' contains the number of radial orders n for each l

	# for each l (each column) populate the matrix with 
	# as many n as there are
	k = 0
	for j in range(lmax+1):
		for i in range(nmax):
			if (i < ind[j]):
				mx_obs[i,j] = fobs[k]
				mx_err[i,j] = err[k]
				k += 1
			else:
				mx_obs[i,j] = 0.
				mx_err[i,j] = 0.

	# returns both matrices
	return mx_obs, mx_err
