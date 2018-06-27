def calc_area(redshift, mu_min, if_plot=False):
    
    import magnif 
    import pylab
    from astropy.io import fits
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.cosmology import FlatLambdaCDM
    
    xdata = []
    ydata = []
    n1 = []

    file1 = 'outmagnif0.5-'+str(redshift)+'.fits'
    print('file', file1)  

    magnif.savemagnif('hlsp_frontier_model_macs1149_sharon_v4_kappa.fits', 'hlsp_frontier_model_macs1149_sharon_v4_gamma1.fits', 0.5, file1, redshift)
    print('done saving')

  
    #Open data file 
    hdulist=fits.open(file1)
    
    (n,bins,params) = pylab.hist(hdulist[0].data[:,:].flatten(), histtype='step', bins=np.arange(10,1000,2))
   
 
    #Build arrays to plot data on log scale
    
    for i in range(len(n)):
        if n[i] == 0:
            a = [0]
        else:
            a = [np.log10(n[i])]
        b = [np.log10(bins[i])]
        ydata = np.concatenate((ydata,a))
        xdata = np.concatenate((xdata,b))
    
    
    #Fit a line to log-scale data
    
    def f(x,k,b):
        return -1*k*x+b
    
    popt = curve_fit(f,xdata,ydata)
    #print(popt)
    array = popt[0]
    if if_plot == True:    
        plt.clf()
        plt.plot(xdata, ydata)
        plt.plot(xdata, f(xdata, array[0], array[1]))
        plt.xlabel('Magnification')
        plt.ylabel('Pixels')
    
    #To plot area vs. magnification
    
    cosmo = FlatLambdaCDM(70,0.3)
    dist, other = str(cosmo.angular_diameter_distance(redshift)).split(' ')
    dist = float(dist)
    
    cd1_1 = hdulist[0].header['CD1_1'] 
    cd2_1 = hdulist[0].header['CD2_1']
    ang = (cd1_1**2 + cd2_1**2)**(0.5)
    
    area = (dist*0.01745329*ang*1000*2)**2 #kpc**2/pixel
    
    for i in range(len(n)):
        a1 = [np.log10(n[i]*area/bins[i])]
        n1 = np.concatenate((n1,a1))
    
    if if_plot == True:
        plt.figure()
        plt.plot(xdata, n1)
        plt.xlabel('Magnification')
        plt.ylabel('Area (Kpc^2)')
        plt.show()
    
    #Calculate integral
    
    mu_min = np.log10(mu_min)
    #print('mu', mu_min)
    sum_mu = 0.
    
    n10 = n1
    for i in range(len(n1)):
        n10[i] = 10**n1[i]
    
    
    if mu_min < abs(xdata[0]):
        print('out of range')
    elif mu_min > abs(xdata[-1]):
        print('out of range')
    else:
        for i in range(len(n)):
            if xdata[i] >= mu_min:
                sum_mu = sum_mu + (n10[i])*((bins[i+1])-(bins[i]))
    
        print('integral:',sum_mu)
    
    
        def func(x,c,d,e):
            return c/((x+d)**e)
    
        popt, other = curve_fit(func,bins[:-1],n10)
        #print(popt)
        const = popt[0]
        offset = popt[1]
        exponent = popt[2]
        #print(const, offset, exponent)
    
        int2 = (const*(bins[-1]+offset)**(-exponent+1))/(-exponent+1) - (const*(10**mu_min + offset)**(-exponent+1))/(-exponent+1)
        print('integral2:',int2)
    
    
        if if_plot == True:
            plt.clf()
            plt.plot(bins[:-1],n10)
            plt.plot(bins[:-1],func(bins[:-1],const,offset,exponent))
            plt.xlabel('Magnification')
            plt.ylabel('Area (Kpc^2)')
            plt.show()

    convert, other = str(cosmo.arcsec_per_kpc_comoving(redshift)).split('a')
    convert = float(convert)
    int2 = (convert**2)*int2
    return int2 
   


 
if __name__ == "__main__":
    import sys
    import numpy as np
    #file1=sys.argv[1]
    redshift = np.float(sys.argv[1])
    mu_min = np.float(sys.argv[2])
    calc_area(redshift, mu_min, if_plot=True)

