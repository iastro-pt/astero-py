import numpy as np
from pylab import *

class _astero:
    """
    Common class for a set of oscillation frequencies
    """

    def __init__(self, filename, filetag, name='undefined'):
        """ 
        Constructs an object that contains the given set of oscillation 
        frequencies, related quantities and methods.

        Parameters
        ----------
        filename : string
                   Name of the file with the frequencies.
                   
        filetag  : integer
                   Specifies the layout of the file. Range 1-6
                   1 - 	l, n, nu_model, nu_corrected
                   2 -  l, n, nu_model, inertia, beta, splitting
                   3 -  l, nu_obs, err_obs, l, nu_model
                   4 -  l, n, nu_model, inertia
                   
                   5 -  l, n, nu_obs, err_obs
                   6 -  l, nu_obs, err_obs, flag1, flag2
                   
        name     : string, optional
                   Name of the star
        """

        
        import read

        self.name = name

        if (filetag == 1): 
	        self.mx, self.mxCorr, self.nvector = read.readModel(filename)
	        self.mxErr = np.zeros(np.shape(self.mx))
	        self.type = 'model'
	        (self.Nn, self.Nl) = np.shape(self.mx)
	
        elif (filetag == 2):
	        # see docstring for read.readModel2 for type of file
	        self.mx, self.nvector = read.readModel2(filename)
	        self.mxErr = np.zeros(np.shape(self.mx))
	        self.type = 'model'
	        (self.Nn, self.Nl) = np.shape(self.mx)

        elif (filetag == 3):
			# see docstring for read.readModel3 for type of file
			self.mx = read.readModel3(filename)
			self.mxErr = np.zeros(np.shape(self.mx))
			self.type = 'model'
			(self.Nn, self.Nl) = np.shape(self.mx)
			
        elif (filetag == 4):
			# see docstring for read.readModel4 for type of file
			self.mx, self.nvector = read.readModel4(filename)
			self.mxErr = np.zeros(np.shape(self.mx))
			self.type = 'model'
			(self.Nn, self.Nl) = np.shape(self.mx)
			
        elif (filetag == 5):
			# see docstring for read.readObs for type of file
			self.mx, self.mxErr, self.nvector = read.readObs(filename)
			self.type = 'observational'
			(self.Nn, self.Nl) = np.shape(self.mx)
			
        elif (filetag == 6):
			# see docstring for read.readObs2 for type of file
			self.mx, self.mxErr = read.readObs2(filename)
			self.type = 'observational'
			(self.Nn, self.Nl) = np.shape(self.mx)

        else:
			print 'The filetag argument must be in the range 1-6.'


    ## Calculate delta_nu
    ## (based on method of White et al. 2011)
    ##
    def delta_nu(self, use_errors=False, plot=False, quiet=False):
        """
        To obtain delta_nu and epsilon, we perform a leastsquares fit to the
        radial (l = 0) frequencies as a function of n. 
        By the asymptotic relation, the gradient of this fit is delta_nu and
        the intercept is delta_nu * epsilon. 
        """

        from linfit import linfit

        f = self.mx[:,0]   # l=0 frequencies
        f = f[f != 0]   # only non-zero ones
        err = self.mxErr[:,0]   # l=0 frequency errors
        err = err[err != 0]
        try:
            n = self.nvector[:len(f)]
        except AttributeError:
            print 'There is no information on the radial orders of the frequencies\n'
            print 'Fit is made assuming consecutive values'
            n = range(1,len(f)+1)

        if(use_errors): a, dnu, dic = linfit(n, f, err=err, use_err=True, full_output=True, plot=plot)
        else: a, dnu, dic = linfit(n,f,full_output=True, plot=plot)

        if (not quiet):
            print 'Fit statistics (nu = a + b*n):'
            print '   b (=dnu) = ' + str(dnu) + ' +/- ' + str(dic['SEb'])
            print '   R^2 = ' + str(dic['r2'])
        self.dnu = dnu



    ## Calculate large separations
    ##
    def large_sep(self):
        """ 
        Calculate large separations, propagating errors 
        results: Dl, DlErr
        """
        import separations as diff
        self.Dl, self.DlErr = diff.large_sep(self.mx, self.mxErr)


    ## Calculate small separations
    ##
    def small_sep(self, lag02=0, lag13=0, lag01=0):
        """ 
        Calculate small separations, propagating errors
        results: d01, d10, d02, d13
        """
        import separations as diff
        self.d02, self.d02Err = diff.small_sep02(self.mx, self.mxErr, lag=lag02)  # (scaled)
        self.d02 = self.d02[0];	self.d02Err = self.d02Err[0]

        try:
            self.d13, self.d13Err    = diff.small_sep13(self.mx, self.mxErr, lag=lag13)  # (scaled)
            self.d13 = self.d13[0];	self.d13Err = self.d13Err[0]
        except (IndexError):
            print 'd13 can''t be calculated'
            pass

        self.d01, self.d01Err = diff.sep01_3(self.mx, self.mxErr, lag=lag01)
        self.d01 = self.d01[0];	self.d01Err = self.d01Err[0]

        self.d10, self.d10Err = diff.sep10_3(self.mx, self.mxErr, lag=lag01)
        self.d10 = self.d10[0];	self.d10Err = self.d10Err[0]


    ##
    ## Calculate ratios
    ##
    def ratios(self, lag02=0, lag13=0, lag01=0):
        import separations as diff
        
        try:
            dnu = self.dnu
        except AttributeError:
            self.delta_nu(quiet=True)
            dnu = self.dnu
        
        self.r02, self.r02Err = diff.ratio02(self.mx, self.mxErr, lag02=lag02, dnu_mean=dnu)
        self.r02 = self.r02[0];	self.r02Err = self.r02Err[0]
        try:
            self.r13, self.r13Err = diff.ratio13(self.mx, self.mxErr, lag=lag13, dnu_mean=dnu)
            self.r13 = self.r13[0];	self.r13Err = self.r13Err[0]
        except (IndexError):
            print 'r13 can''t be calculated'
            pass    
        self.r01, self.r01Err = diff.ratio01(self.mx, self.mxErr, lag=lag01)
        self.r01 = self.r01[0];	self.r01Err = self.r01Err[0]
        self.r10, self.r10Err = diff.ratio10(self.mx, self.mxErr, lag=lag01)
        self.r10 = self.r10[0];	self.r10Err = self.r10Err[0]



    ##
    ## Calculate the second differences
    ##
    def sec_dif(self):
    	import separations as diff
    	
    	try:
            dnu = self.dnu
        except AttributeError:
            self.delta_nu(quiet=True)
            dnu = self.dnu 
        
        self.D2, self.D2Err = diff.second_dif(self.mx, self.mxErr, dnu_mean=dnu)
        #print 'Second differences done'


	##
	## Plot the large separations as a function of frequency
	##
    def plotls(self, xlim=[None,None], ylim=[None,None], fig=1):

        try: 
            self.Dl
        except AttributeError:
            print 'need to run "large_sep()" before'
            return
            
        fig = figure(fig)
        ax = subplot(111)
        for i in range(self.Nl):
            lbl = 'l='+str(i)
            ax.errorbar(self.mx[0:-1,i], self.Dl[:,i], yerr=self.DlErr[:,i], fmt='-o', label=lbl)
        leg=ax.legend()
        ax.grid(False)
        ax.set_xlabel('Frequency')
        ax.set_ylabel(r'$\Delta_{\ell} \, (\mu Hz)$')

        ax.set_title(r'Large Separation $\Delta_l$')
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        show()				

	##
	## Plot the (scaled) small separations as a function of frequency
	##
    def plotss(self, xlim=[None,None], ylim=[None,None], fig=1, only='all'):

        oplot = False
        if (fignum_exists(fig)):
            oplot = True
            symbols = ['d', '+', '*']
            index = np.random.randint(0,len(symbols))
            symbol = symbols[index]
        else:
            symbol = 'o'
            
        fig = figure(fig)
        ax = subplot(111)
        
        if (only == 'all'):
		    ax.errorbar(self.mx[:,0],self.d02, yerr=self.d02Err, fmt='r-'+symbol)
		    ax.errorbar(self.mx[:,1],self.d13, yerr=self.d13Err, fmt='b-'+symbol)
		    ax.errorbar(self.mx[:,0],self.d01, yerr=self.d01Err, fmt='g-'+symbol)
		    ax.errorbar(self.mx[0:-1,1],self.d10, yerr=self.d10Err, fmt='k-'+symbol)
		    leg = ax.legend((r'$d_{02}/3$',r'$d_{13}/5$',r'$d_{01}$',r'$d_{10}$'), 'upper right')
        elif (only == 'd02'):
            ax.errorbar(self.mx[:,0],self.d02, yerr=self.d02Err, fmt='r-'+symbol)
            leg = ax.legend('d02 / 3', 'upper right')
        elif(only == 'd13'):
            ax.errorbar(self.mx[:,1],self.d13, yerr=self.d13Err, fmt='b-'+symbol)
            leg = ax.legend('d13 / 5', 'upper right')
        elif(only == 'd010'):
            mfc = 'orange' if oplot else None
            ms = 2 if oplot else None
            ax.errorbar(self.mx[:,0],self.d01, yerr=self.d01Err, fmt='g'+symbol, mfc=mfc)
            ax.errorbar(self.mx[0:-1,1],self.d10, yerr=self.d10Err, fmt='k'+symbol, mfc=mfc)
            leg = ax.legend((r'$d_{01}$',r'$d_{10}$'), 'upper right')
        
        ax.grid(False)
        ax.set_xlabel('Frequency')
        ax.set_ylabel(r'$d \, (\mu Hz)$')

        ax.set_title(r'Small Separations')
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        plt.show()	

	##
	## Plot the ratios as a function of frequency
	##
    def plotr(self, xlim=[None,None], ylim=[None,None], fig=1, only='all'):

        oplot = False
        if (fignum_exists(fig)):
            oplot = True
            symbols = ['d', '+', '*']
            index = np.random.randint(0,len(symbols))
            symbol = symbols[index]
        else:
            symbol = 'o'
            
        fig = figure(fig)
        ax = subplot(111)
        
        if (only == 'all'):
            ax.errorbar(self.mx[0:-1,0],self.r02, yerr=self.r02Err, fmt='r-'+symbol)
            ax.errorbar(self.mx[0:-1,1],self.r13, yerr=self.r13Err, fmt='b-'+symbol)
            ax.errorbar(self.mx[:,0],self.r01, yerr=self.r01Err, fmt='g-'+symbol) 
            ax.errorbar(self.mx[0:-1,1],self.r10, yerr=self.r10Err, fmt='k-'+symbol) 
            leg = ax.legend((r'$r_{02}$',r'$r_{13}$',r'$r_{01}$',r'$r_{10}$'))
        elif (only == 'r02'):
            ax.errorbar(self.mx[0:-1,0],self.r02, yerr=self.r02Err, fmt='r-'+symbol)
            leg = ax.legend((r'$r_{02}$',r'$r_{13}$',r'$r_{01}$',r'$r_{10}$'))
        elif(only == 'r13'):
            ax.errorbar(self.mx[0:-1,1],self.r13, yerr=self.r13Err, fmt='b-'+symbol)
            leg = ax.legend((r'$r_{02}$',r'$r_{13}$',r'$r_{01}$',r'$r_{10}$'))
        elif(only == 'r010'):
            mfc = 'orange' if oplot else None
            ms = 2 if oplot else None
            ax.errorbar(self.mx[:,0],self.r01, yerr=self.r01Err, fmt='g'+symbol, mfc=mfc)
            ax.errorbar(self.mx[0:-1,1],self.r10, yerr=self.r10Err, fmt='k'+symbol, mfc=mfc)
            leg = ax.legend((r'$r_{02}$',r'$r_{13}$',r'$r_{01}$',r'$r_{10}$'))
            
        ax.grid(False)
        ax.set_xlabel('Frequency')
        ax.set_ylabel('ratio')

        ax.set_title(r'Ratios')
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        plt.show()


	##
	## Plot the dr as a function of frequency (stub!)
	##
	def plotdr(self, xlim=[None,None], ylim=[None,None]):
		if (self.type == 'm'):
			if (self.star=='A'): fig = plt.figure(5)
			else: fig = plt.figure(9)
			
			ax = fig.add_subplot(211)
			ax.plot(self.mx[0:-1,0],(self.r02 - self.r13),'r-o')
			leg = ax.legend((r'$dr$'), 'upper left')
			ax.grid(False)
			ax.set_xlabel('Frequency')
			ax.set_ylabel(r'$dr \, (\mu Hz)$')

			if (self.star=='A'): ax.set_title(r'(modelo;*A)')
			else: ax.set_title(r'(modelo;*B)')

			ax.set_ylim(ylim)
			ax.set_xlim(xlim)

			plt.show()

		else:
			if (self.star=='A'): fig = plt.figure(5)
			else: fig = plt.figure(9)
			
			ax = fig.add_subplot(212)
			ax.plot(self.mx[0:-1,0],self.r02-self.r13,'r-o')
			leg = ax.legend((r'$dr$'), 'upper left')
			ax.grid(False)
			ax.set_xlabel('Frequency')
			ax.set_ylabel(r'$dr \, (\mu Hz)$')
			
			if (self.star=='A'): ax.set_title(r'(observacoes;*A)')
			else: ax.set_title(r'(observacoes;*B)')
			
			ax.set_ylim(ylim)
			ax.set_xlim(xlim)
			plt.show()



	##
	## Plot the second differences as a function of frequency (stub!)
	##
	def plotsd(self, xlim=[None,None], ylim=[None,None]):
		
		if (self.type == 'm'):
			if (self.star=='A'): fig = plt.figure(6)
			else: fig = plt.figure(10)

			ax = fig.add_subplot(111)
			ax.errorbar(self.mx[0:-1,:], self.D2[:,:], fmt='-o')
			leg = ax.legend(('l=0', 'l=1', 'l=2', 'l=3'), 'lower right')
			ax.grid(False)
			ax.set_xlabel('Frequency')
			ax.set_ylabel(r'$\Delta_{\ell} \, (\mu Hz)$')

			if (self.star=='A'): ax.set_title(r'Second differences $\Delta_2$ (modelo;*A)')
			else: ax.set_title(r'Second differences $\Delta_2$ (modelo;*B)')

			ax.set_ylim(ylim)
			ax.set_xlim(xlim)
			plt.show()
		else:
			if (self.star=='A'): fig = plt.figure(6)
			else: fig = plt.figure(10)

			ax = fig.add_subplot(111)
			for i in range(self.Nl):
				lbl = 'l='+str(i)
				ax.errorbar(self.mx[0:-2,i], self.D2[0:-1,i], yerr=self.D2Err[0:-1,i], fmt='-o', label=lbl)
			leg=ax.legend(loc='lower center')
			ax.grid(False)
			ax.set_xlabel('Frequency')
			ax.set_ylabel(r'$\Delta_{\ell} \, (\mu Hz)$')

			if (self.star=='A'): ax.set_title(r'Large Separation $\Delta_l$ (observacoes;*A)')
			else:  ax.set_title(r'Large Separation $\Delta_l$ (observacoes;*B)')

			ax.set_ylim(ylim)
			ax.set_xlim(xlim)
			plt.show()	

	##
	## Relative differences to observations (stub!)
	##
	def relative (self, obs, plot=True, xlim=[None,None], ylim=[None,None]):

		self.rel1 = (self.Dl- obs.Dl) / obs.Dl
		self.rel2 = (self.d02- obs.d02) / obs.d02
		self.rel3 = (self.d13- obs.d13) / obs.d13

		if (plot):
			#for i in range(3):
				fig = plt.figure(i)
				ax = fig.add_subplot(111)
				ax.plot(self.rel, self.D2[:,:], fmt='-o')
				#leg = ax.legend(('l=0', 'l=1', 'l=2', 'l=3'), 'lower right')
				ax.grid(False)
				#ax.set_xlabel('Frequency')
				#ax.set_ylabel(r'$\Delta_{\ell} \, (\mu Hz)$')
				#ax.set_title(r'Second differences $\Delta_2$ (modelo;*A)')
				ax.set_ylim(ylim)
				ax.set_xlim(xlim)
				plt.show()


    ##
    ## Plot the echelle diagram
    ##
    def echelle (self, dnu=None, nu0=0.0, fig=None):
        """
        Plot the echelle diagram of the object frequencies.
        
        Parameters
        ----------
        dnu : mean large separation. By default, delta_nu() is run to calculate
              it from l=0 frequencies. The user can input any other value.
        nu0 : reference frequency
        fig : identification of the Figure window for overplotting (experimental)
        """
        
        import numpy as np
        from pylab import figure, subplot, errorbar, minorticks_on, show
        from pylab import setp, gca
        from pylab import fignum_exists
        if (dnu == None):
            try:
                dnu = self.dnu
            except AttributeError:
                self.delta_nu()
                dnu = self.dnu
        
        if (fignum_exists(fig)):
            colors = ['b', 'r', 'g', 'y']
            index = np.random.randint(0,len(colors))
            color = colors[index]
        else:
            color = 'k'
        face = 'w' if color=='k' else color
            
        figure(fig)
        ax = subplot(111)
        nuMODdnu = self.mx % (dnu - nu0)
        
        l0,cl,bl = errorbar(nuMODdnu[:,0], self.mx[:,0], 
                    xerr=self.mxErr[:,0], fmt=color +'*', ms=8, capsize=2)
        l1,cl,bl = errorbar(nuMODdnu[:,1], self.mx[:,1], 
                    xerr=self.mxErr[:,1],fmt=color +'o', mfc=face, ms=6)
        l2,cl,bl = errorbar(nuMODdnu[:,2], self.mx[:,2], 
                    xerr=self.mxErr[:,2],fmt=color +'d', mfc=face, ms=7,capsize=2)
        l3,cl,bl = errorbar(nuMODdnu[:,3], self.mx[:,3], 
                    xerr=self.mxErr[:,3],fmt=color +'s', mfc=face, ms=6,capsize=3)

        ax.legend((r'$\ell=0$', r'$\ell=1$', r'$\ell=2$', r'$\ell=3$'),
                  frameon=False, numpoints=1)
        ax.set_xlabel(r'$\nu$  ' + 'mod' + r' $\Delta\nu \, (\mu Hz)$',
                      fontsize=20, labelpad=15)
        ax.set_ylabel(r'$\nu \, (\mu Hz)$', fontsize=20, labelpad=15)
        ax.set_title(r'$\Delta \nu =$' + str(dnu)[0:6], fontsize=22)
        #ax.set_ylim([1600,3500])
        #ax.set_xlim([20,140])

        setp(ax.get_xticklabels(), fontsize=18)
        setp(ax.get_yticklabels(), fontsize=18)
        leg = gca().get_legend()
        ltext  = leg.get_texts()
        setp(ltext, fontsize=18) 
        minorticks_on()
        show()

    
    ##
    ## Output frequencies to MESA inlist form
    ##
    def mesa_create_inlist(self, max_error=None):
    
        print 'nl0 = ' + str( (self.mx[:,0] != 0.).sum() )
        for i in range(len(self.mx[:,0])):
            if (self.mx[i,0] != 0):
                print 'l0_obs(' + str(i+1) + ') = ' + str(self.mx[i,0]) + 'd0'
                print 'l0_obs_sigma(' + str(i+1) + ') = ' + str(self.mxErr[i,0]) + 'd0'
        
        print '\n'
        
        print 'nl1 = ' + str( (self.mx[:,1] != 0.).sum() )
        for i in range(len(self.mx[:,1])):
            if (self.mx[i,1] != 0):
                print 'l1_obs(' + str(i+1) + ') = ' + str(self.mx[i,1]) + 'd0'
                print 'l1_obs_sigma(' + str(i+1) + ') = ' + str(self.mxErr[i,1]) + 'd0'
        
        print '\n'
        
        print 'nl2 = ' + str( (self.mx[:,2] != 0.).sum() )
        for i in range(len(self.mx[:,2])):
            if (self.mx[i,2] != 0):
                print 'l2_obs(' + str(i+1) + ') = ' + str(self.mx[i,2]) + 'd0'
                print 'l2_obs_sigma(' + str(i+1) + ') = ' + str(self.mxErr[i,2]) + 'd0'
    
        print '\n'
        
        print 'nl3 = ' + str( (self.mx[:,3] != 0.).sum() )
        for i in range(len(self.mx[:,3])):
            if (self.mx[i,3] != 0):
                print 'l3_obs(' + str(i+1) + ') = ' + str(self.mx[i,3]) + 'd0'
                print 'l3_obs_sigma(' + str(i+1) + ') = ' + str(self.mxErr[i,3]) + 'd0'








