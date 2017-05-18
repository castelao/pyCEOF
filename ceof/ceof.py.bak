#!/usr/bin/env python
# -*- coding: Latin-1 -*-

""" Class to deal with Complex EOF
"""

from UserDict import UserDict

import numpy
import numpy as np
from numpy import ma
import scipy.fftpack

from pyclimate.svdeofs import svdeofs, getvariancefraction

from utils import scaleEOF


def ceof_scalar2D(data):
    """ Estimate the complex EOF on a 2D array.

        Time should be the first dimension, so that the PC (eigenvalues) will
          be in respect to the first dimension.

    """
    assert type(data) is np.ndarray, \
        "ceof_scalar2D requires an ndarray but got: %s" % type(data)
    assert np.isfinite(data).all(), \
        "ceof_scalar2D requires a full valid values array"

    # ---- Creating the complex field using Hilbert transform
    input_H = numpy.empty(data.shape, dtype=data.dtype)
    for i in range(data.shape[1]):
        input_H[:,i] = scipy.fftpack.hilbert(data[:,i])

    U = data + 1j*input_H

    pcs, lambdas, eofs = svdeofs(U)

    return pcs, lambdas, eofs


def CEOF_2D(data, cfg=None):
    """ Complex EOF of a scalar 2D array
    """
    assert data.ndim == 2, "CEOF_2D requires a 2D ndarray input"

    if cfg is None:
        cfg = {'cumvar':1, 'normalize':'pc_median'}
    assert type(cfg) is dict, "cfg must be a dictionary"

    # ---- Creating the complex field using Hilbert transform
    #if self['input'].mask.any():
    #    print "There are masked values in U at CEOF_2D()"

    pcs, lambdas, eofs = ceof_scalar2D(data)

    expvar = getvariancefraction(lambdas)

    # Define how many modes will be returned by the explainned variance.
    # cumvar = 1 means 100%, i.e. all modes
    if cfg['cumvar'] == 1:
        nmodes = len(lambdas)
    else:
        # This doesn't work. Need to improve it.
        nmodes = (np.arange(len(expvar))[np.cumsum(expvar)>=cfg['cumvar']])[0]

    if 'maxnmodes' in cfg:
        nmodes = min(nmodes, cfg['maxnmodes'])

    print "Considering the first %s of %s modes." % (nmodes,len(lambdas))

    # ---- Normalize -----------------------------------------------------
    if 'normalize' in  cfg:
        pcs, eofs = scaleEOF(pcs, eofs, scaletype=cfg['normalize'])

    output = {}
    output['eofs'] = eofs[:,:nmodes]
    output['pcs'] = pcs[:,:nmodes]
    output['lambdas'] = lambdas[:nmodes]
    output['variancefraction'] = expvar[:nmodes]

    return output


def make_animation(data, eofdata, t, lat, lon, outputfilename, limits = None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt   # For plotting graphs.
    import pylab
    from mpl_toolkits.basemap import Basemap

    #import numpy as np
    #import subprocess                 # For issuing commands to the OS.
    import os
    import sys                        # For determining the Python version.

    not_found_msg = """
    The mencoder command was not found;
    mencoder is used by this script to make an avi file from a set of pngs.
    It is typically not installed by default on linux distros because of
    legal restrictions, but it is widely available.
    """

    if limits == None:
        LatIni = lat.min()
        LatFin = lat.max()
        LonIni = lon.min()
        LonFin = lon.max()
    else:
        LatIni = limits['LatIni']
        LatFin = limits['LatFin']
        LonIni = limits['LonIni']
        LonFin = limits['LonFin']

    #try:
    #    subprocess.check_call(['mencoder'])
    #except subprocess.CalledProcessError:
    #    print "mencoder command was found"
    #    pass # mencoder is found, but returns non-zero exit as expected
    #         # This is a quick and dirty check; it leaves some spurious output
    #	     # for the user to puzzle over.
    #except OSError:
    #	 print not_found_msg
    #     sys.exit("quitting\n")

    parallels = numpy.arange(-5,20.1,5)
    meridians = numpy.arange(300,340,10)
    V = range(-20, 21, 1)
    for i in range(len(t)) :
        #fig = plt.figure(figsize=(14,10.5), dpi=100)
        pylab.subplot(211)
        map = Basemap(projection='merc',lat_ts=0,llcrnrlon=LonIni,llcrnrlat=LatIni, urcrnrlon=LonFin, urcrnrlat=LatFin,resolution='l',area_thresh=1000.)
        X, Y = map(*pylab.meshgrid(lon, lat))
        map.contourf(X,Y, data[i], V)
        plt.colorbar(shrink=0.8)
        map.drawcoastlines()
        map.fillcontinents(color='0.0')
        map.drawparallels(parallels,labels=[1,0,0,1])
        map.drawmeridians(meridians,labels=[1,0,0,1])
        plt.title('%s' % (t[i].strftime('%Y-%m-%d')))
        #
        pylab.subplot(212)
        map = Basemap(projection='merc',lat_ts=0,llcrnrlon=LonIni,llcrnrlat=LatIni, urcrnrlon=LonFin, urcrnrlat=LatFin,resolution='l',area_thresh=1000.)
        X, Y = map(*pylab.meshgrid(lon, lat))
        map.contourf(X,Y, eofdata[i], V)
        plt.colorbar(shrink=0.8)
        map.drawcoastlines()
        map.fillcontinents(color='0.0')
        map.drawparallels(parallels,labels=[1,0,0,1])
        map.drawmeridians(meridians,labels=[1,0,0,1])
        plt.title('%s (EOF reconstructed)' % (t[i].strftime('%Y-%m-%d')))
        # ----
        filename = str('../tmp/%04d' % i) + '.png'
        plt.savefig(filename, dpi=100)
        print 'Wrote file', filename
        #
        # Clear the figure to make way for the next image.
        #
        plt.clf()

    command = ('mencoder',
               'mf://../tmp/*.png',
               '-mf',
               'type=png:w=800:h=600:fps=2',
               '-ovc',
               'lavc',
               '-lavcopts',
               'vcodec=mpeg4',
               '-oac',
               'copy',
               '-o',
               outputfilename)
    
    os.spawnvp(os.P_WAIT, 'mencoder', command)

    #print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
    #subprocess.check_call(command)

    print "\n\n The movie was written to 'output.avi'"

    print "\n\n You may want to delete *.png now.\n\n"
    return


class CEOF_Filter():
    """ Unfinished
    
        This will be a class to filter using EOF. One define the number of modes
	  or the desired variability to be explainned. The field is decomposed
	  by CEOF, than the field is reconstructed considering only the n first
	  modes.



    """
    def __init__(self,input,metadata={'variancefraction_explainned':0.95}):
        """
	"""
	N = pcs.shape[1]    # Now it is all the modes.
        T = pcs.shape[0]
	I = eofs.shape[0]
	J = eofs.shape[1]
        data_filtered = numpy.zeros((T, I, J))

        eof_amp=(eof.real**2+eof.imag**2)**0.5
        eof_phase=numpy.arctan2(eof.imag,eof.real)
        pc_amp = (numpy.real(pc)**2+numpy.imag(pc)**2)**0.5
        pc_phase = numpy.arctan2(numpy.imag(pc),numpy.real(pc))	

	for t in range(T):
            for n in range(N):
	        data_filtered[t] = data_filtered[t] + eof_amp[:,:,n]*pc_amp[t,n]*numpy.cos(eof_phase[:,:,n]+pc_phase[t,n])
        return

def CEOF_2D_limited(input,metadata={'variancefraction_explainned':0.95}):
        #self.data['input']=input.copy
        ceof = CEOF_2D(input)
        if 'nmodes' not in metadata:
            nmodes=len(ceof['variancefraction'])
        if 'variancefraction_explainned' in metadata:
            nmodes=(numpy.ones(ceof['variancefraction'].shape)[numpy.cumsum(ceof['variancefraction'])<=metadata['variancefraction_explainned']]).sum().astype('i')
        if 'nmodes_max' in metadata:
            if metadata['nmodes_max']<nmodes:
                nmodes=metadata['nmodes_max']
        nt,ni,nj=input.shape
        filtred=numpy.zeros((nt,ni,nj))
        for t in range(nt):
            #for n in range(20):
            for n in range(ceof['pcs'].shape[1]):
                filtred[t,:,:]+=ceof['eofs'][:,:,n].real*x['pcs'][t,n].real
        #self.data['filtred']=filtred
        return filtred


class CEOF(UserDict):
    """
    """
    def __init__(self, input, metadata={}, logger=None, **keywords):
        """
            Time should be the first dimension, i.e. axis=0
        """
        self.input = input.copy()
        self.data = input.copy()
        self.metadata = metadata

        self.go()
        return

    def filter(self,var,l,type,l2=None):
        #from maud import window_mean
        from maud import window_1Dmean_grid, get_halfpower_period
        from datetime import timedelta
        
        if len((set(numpy.diff(self.data['datetime'])))) !=1:
            print "Class incomplete. Can't deal with a non regular time series"
            return

        dt=self.data['datetime'][1]-self.data['datetime'][0]
        if type == 'bandpass':
            tscale = dt.days+dt.seconds/86400.
            ll = (l.days+l.seconds/86400.)/tscale
            ll2 = (l2.days+l2.seconds/86400.)/tscale
            #lowpass = window_mean.window_1Dmean_grid(self.data[var], ll/2., method='hann', axis=0)
            lowpass = window_1Dmean_grid(self.data[var], ll/2., method='hann', axis=0)
            output = window_1Dmean_grid(lowpass, ll2/2., method='hann', axis=0)
            output = lowpass - output

            print "ATENTION!!!! Improve this here!!!"
            self.halfpower_period = "20-120"
        else:
            ll=(l.days+l.seconds/86400.)/(dt.days+dt.seconds/86400.)

            if ll<1:
                print "This filter will have no effect. Data series have not enough resolution."
                return

            lowpass = window_1Dmean_grid(self.data[var],ll/2.,method='hann',axis=0)
            #lowpass=window_mean.window_1Dmean(self.data[var],ll,method='hanning',axis=0)
            if type=='lowpass':
                output=lowpass
            elif type=='highpass':
                #self.data[var]=x_highpass=self.data[var]-(lowpass-lowpass.mean())
                output=self.data[var]-(lowpass)
                #output=self.data[var]-lowpass
            else:
                print "On function filter, type must be lowpass or highpass"

            halfpower_period = get_halfpower_period(self.data[var],
                    output, dt=dt)
            print "Filter half window size: %s" % l
            print "Half Power Period: %s" % halfpower_period
            self.halfpower_period = halfpower_period
        
        # ----
	## I should move this to inside the window_mean_1D_grid
	#nt,ni,nj = self.data[var].shape
	#gain = ma.masked_all((nt,ni,nj))
	#for i in range(ni):
	#    for j in range(nj):
	#        if output[:,i,j].mask.all()==False:
        #            gain[:,i,j] = numpy.absolute(numpy.fft.fft(output[:,i,j]-output[:,i,j].mean())) / numpy.absolute(numpy.fft.fft(self.data[var][:,i,j]-self.data[var][:,i,j].mean()))
	#gain_median = ma.masked_all(nt)
	#for t in range(nt):
	#    gain_median[t] = numpy.median(gain[t,:,:].compressed()[numpy.isfinite(gain[t,:,:].compressed())])
	#freq=numpy.fft.fftfreq(nt)/dt.days
	#import rpy2.robjects as robjects
	#smooth = robjects.r['smooth.spline'](robjects.FloatVector(gain_median[numpy.ceil(nt/2.):]),robjects.FloatVector(-freq[numpy.ceil(nt/2.):]),spar=.4)
	##smooth = robjects.r['smooth.spline'](robjects.FloatVector(-freq[numpy.ceil(nt/2.):]),robjects.FloatVector(gain_median[numpy.ceil(nt/2.):]),spar=.4)
	#s_interp = robjects.r['predict'](smooth,x=0.5)
	#halfpower_period = 1./s_interp.rx2['y'][0]

        # ----
        self.data[var]=output

        return

    def select_data(self, var, polygon_coordinates):
        """
        """
        #var = 'ssh'
        T,I,J = self.data[var].shape
        tmp=numpy.ones((J,K))==1
        for i in range(I):
            for j in range(J):
                   tmp[i,j] = ((self.data[var].mask)[:,i,j]).all()==False


        from shapely.geometry import Polygon
        from shapely.geometry import Point
        polygon = Polygon(polygon_coordinates)

        ind = ind&tmp
        return ind



    def set_wavelenght(self):
        """ Estimate the wavelenghts from the gradient of the EOF


	"""

        eof_phase=numpy.arctan2(self['eofs'].imag, self['eofs'].real)

        eof_phase_360 = eof_phase.copy()
        eof_phase_360[eof_phase<0] = 2*numpy.pi+eof_phase[eof_phase<0]

        from fluid.common.common import _diff_centred
        dx_eof_phase = _diff_centred(eof_phase,dim=1)
        dx_eof_phase_360 = _diff_centred(eof_phase_360, dim=1)

        ind = abs(dx_eof_phase)>abs(dx_eof_phase_360)
        dx_eof_phase[ind] = dx_eof_phase_360[ind]

	self.data['dx_eof_phase'] = dx_eof_phase

        #from scipy.interpolate import bisplrep, bisplev
        #tck = bisplrep(x['Lon'], x['Lat'], eof_phase)
        #dx_eof_phase_spline = bisplev(x['Lon'][0,:], x['Lat'][:,0],tck,dx=1)#/self.data['dX']

        from fluid.common.common import lonlat2dxdy
        #dX, dY = lonlat2dxdy( x['Lon'][0,:], self['Lat'][:,0])
        dX, dY = lonlat2dxdy( self['lon'], self['lat'])

        L_x = ma.masked_all(dx_eof_phase.shape)
	for n in range(dx_eof_phase.shape[-1]):
	    L_x[:,:,n] = dX/dx_eof_phase[:,:,n]*2*numpy.pi*1e-3

        #self.data['L_x'] = dX/dx_eof_phase*2*numpy.pi*1e-3
        self.data['L_x'] = L_x


        #from fluid.common.common import _diff_centred
        #eof_phase=numpy.arctan2(x['eofs'].imag,x['eofs'].real)
        #deof_x=_diff_centred(eof_phase,dim=1)
        #pylab.contourf(eof_phase[:,:,0])
        #pylab.figure()
        #pylab.contourf(deof_x[:,:,0],numpy.linspace(-0.5,0.5,20))
        #pylab.colorbar()
        #pylab.show()

        return


    def go(self):
        var = self.metadata['ceof']['var']

        if ('Lat' not in self.keys()) or ('Lon' not in self.keys()):
            self.data['Lon'], self.data['Lat'] = numpy.meshgrid(self.data['lon'],self.data['lat'])
        if 'prefilter' in self.metadata:
            print "Filtering in time"
            if self.metadata['prefilter']['type'] == 'bandpass':
                self.filter(var,l=self.metadata['prefilter']['l'],type=self.metadata['prefilter']['type'], l2=self.metadata['prefilter']['l2'],)
            else:
                self.filter(var,l=self.metadata['prefilter']['l'],type=self.metadata['prefilter']['type'])
        # ---- Normalize -----------------------------------------------------
        #self.data['ssh']=self.data['ssh']-self.data['ssh'].mean()
        # --------------------------------------------------------------------
        ind = ma.getmaskarray(self.data[var]).any(axis=0) == False

        I, J, K = self.data[var].shape
    
        if 'ceof_coord' in self.metadata:
            coord = self.metadata['ceof_coord']
            assert type(coord) == list
            from shapely.geometry import Point, Polygon
            polygon = Polygon(coord)
            for j in range(J):
                for k in range(K):
                        ind[j,k] = ind[j,k] & polygon.intersects(
                                Point(self.data['Lon'][j,k],
                                    self.data['Lat'][j,k]))
  
        #ind = ind&tmp
    
        N = int((numpy.ones(ind.shape)[ind]).sum())
        grid_index = ma.masked_all((N,2), dtype='i')
        n = -1
        for j in range(J):
            for k in range(K):
                if ind[j, k]:
                    n += 1
                    grid_index[n, :] = [j, k]
    
        self.grid_index = grid_index
        data2D = numpy.zeros((I,N), dtype=self.data[var].dtype)
        for n, ind in enumerate(self.grid_index):
            data2D[:, n] = self.data[var][:,ind[0], ind[1]]

        print "Running CEOF_2D()"
        output = CEOF_2D(data2D, cfg=self.metadata['ceof'])

        nmodes = len(output['lambdas'])

        for k in [k for k in output.keys() if k is not 'ceof']:
            self.data[k] = output[k]

        self.data['eofs'] = ma.masked_all((J, K, nmodes), dtype='complex128')
        for n,ind in enumerate(self.grid_index):
            self.data['eofs'][ind[0],ind[1],:] = output['eofs'][n,:]

	# ----
	self.set_wavelenght()

        #for k in self.data.keys():
        #    print k, self.data[k].shape, type(self.data[k])

        if 'figs' in self.metadata:
	    print "Creating figures for %s modes" % nmodes
            #for nmode in range(5):
            #for n in range(10):
            for n in range(nmodes):
                #eof2D=ceof['eofs'][:,nmode]
                #pc2D=ceof['pcs'][:,nmode]
                ##
                ##eof=ma.masked_all(self.data['ssh'].shape[1:],dtype=eofs.dtype)
                #eof=ma.masked_all(self.data['ssh'].shape[1:],dtype='complex128')
                ##pc=ma.masked_all(self.data['ssh'].shape[1:],dtype=pcs.dtype)
                #pc=ma.masked_all(self.data['ssh'].shape[1:],dtype='complex128')
                ##
                ##for n,i in enumerate(ind):
                ##
                ## ---- 2D back to grid ----
                #for n,ind in enumerate(self.data2D['grid_index']):
                #    #print n,ind
                #    eof[ind[0],ind[1]] = eof2D[n]
                #    #pc[ind[0],ind[1]] = pc2D[n]
                #pc = pc2D
                #print pc2D.shape
                #varfrac = round(getvariancefraction(lambdas)[nmode]*1e2)
                #varfrac = ceof['variancefraction'][nmode]
                #fig = self.plot(eof_amp,eof_phase,pc_amp,pc_phase,nmode,varfrac)
                if 'suffix' in self.metadata['figs']:
                    filename="../figs/CEOF_%s_mode%s.eps" % (self.metadata['figs']['suffix'],(n+1))
                else:
                    filename="../figs/CEOF_mode%s.eps" % (n+1)
                limits={'LatIni':-5, 'LatFin':15, 'LonIni':-60, 'LonFin':-25}
                #self.plot(self['eofs'][:,:,n],self['pcs'][:,n],(n+1),self['variancefraction'][n],filename=filename,limits=limits,cumvarfrac=self['variancefraction'][:(n+1)].sum())
                import graphics
                graphics.plot(self['eofs'][:,:,n], self['pcs'][:,n],
                        (n+1),self['variancefraction'][n],
                        filename=filename,
                        data = self.data,
                        limits=limits,
                        cumvarfrac=self['variancefraction'][:(n+1)].sum())
                #fig.show()

                #import pylab
                #pylab.savefig("../fig/CEOF_mode%s.eps" % nmode)
                ##print "dir(fig)",dir(fig)
                ##fig.close()
                #pylab.close()



        # --------------------------------------------------------------------
        #y=x['eofs'][:,:,0]
        #
        #eofs = x['eofs'][:,:,3]
        #
        #
        #eofs_amp = numpy.absolute(eofs)
        #eofs_phase=numpy.arctan2(eofs.imag,eofs.real)
        ##pcs_phase=numpy.arctan2(pcs.imag,pcs.real)
        #
        #dx_eofs_phase=ma.masked_all(eofs_phase.shape)
        #dy_eofs_phase=ma.masked_all(eofs_phase.shape)
        #
        ##dx_eofs_phase[:,0,:]=eofs_phase[:,1,:]-eofs_phase[:,0,:]
        ##dx_eofs_phase[:,:,1:-1]=eofs_phase[:,:,2:]-eofs_phase[:,:,:-2]
        ##dx_eofs_phase[:,-1,:]=eofs_phase[:,-1,:]-eofs_phase[:,-2,:]
        ##ind=(dx_eofs_phase>3)|(dx_eofs_phase<-3)
        ##dx_eofs_phase.mask[ind]=True
        #
        #
        #dx_eofs_phase[:,1:-1] = (eofs_phase[:,2:]-eofs_phase[:,:-2])/2.
        #ind = ((numpy.sign(eofs_phase[:,2:])*numpy.sign(eofs_phase[:,:-2]))<0) & (numpy.absolute(dx_eofs_phase[:,1:-1])>2)
        #dx_eofs_phase.mask[ind] = True
        #
        #dy_eofs_phase[1:-1,:] = (eofs_phase[2:,:]-eofs_phase[:-2,:])/2.
        #ind = ((numpy.sign(eofs_phase[2:,:])*numpy.sign(eofs_phase[:-2,:]))<0) & (numpy.absolute(dy_eofs_phase[1:-1,:])>2)
        #dy_eofs_phase.mask[ind] = True
        # --------------------------------------------------------------------
        #eofs_phase=numpy.arctan2(eofs.imag,eofs.real)
        #pcs_phase=numpy.arctan2(pcs.imag,pcs.real)

        #dx_eofs_phase=ma.masked_all(eofs_phase.shape)

        #dx_eofs_phase[:,0,:]=eofs_phase[:,1,:]-eofs_phase[:,0,:]
        #dx_eofs_phase[:,1:-1,:]=eofs_phase[:,2:,:]-eofs_phase[:,:-2,:]
        #dx_eofs_phase[:,-1,:]=eofs_phase[:,-1,:]-eofs_phase[:,-2,:]
        #ind=(dx_eofs_phase>3)|(dx_eofs_phase<-3)
        #dx_eofs_phase.mask[ind]=True



        #dphase=ma.masked_all(pcs_phase.shape)


        #dphase[0,:]=(pcs_phase[1,:]-pcs_phase[0,:])
        #dphase[1:-1,:]=(pcs_phase[2:,:]-pcs_phase[:-2,:])
        #dphase[-1,:]=(pcs_phase[-1,:]-pcs_phase[-2,:])

        #ind = (numpy.absolute(dphase[1:-1,:])>2) & (pcs_phase[:-2,:]<0) & (pcs_phase[2:,:]>0)
        #dphase[1:-1,:][ind]=(pcs_phase[2:,:][ind]-(2*numpy.pi+pcs_phase[:-2,:][ind]))
        ##ind = (numpy.absolute(dphase[1:-1,:])>2) & (pcs_phase[:-2,:]>0) & (pcs_phase[2:,:]<0)
        ##dphase[1:-1,:][ind]=(2*numpy.pi+pcs_phase[2:,:][ind])-pcs_phase[:-2,:][ind]

        #ind = (numpy.absolute(dphase[1:-1,:])>3) & (pcs_phase[:-2,:]<0) & (pcs_phase[2:,:]<0)
        #ind = (numpy.absolute(dphase[1:-1,:])>3) & (pcs_phase[:-2,:]>0) & (pcs_phase[2:,:]<0)
        # --------------------------------------------------------------------
