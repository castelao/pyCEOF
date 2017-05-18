#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def plot(eof, pc, nmode, varfrac, filename, data, limits=None, cumvarfrac=None):
    """ Plot one mode of the CEOF
    """
    import pylab
    import matplotlib

    if limits == None:
        #LatIni = self['Lat'].min()
        #LatFin = self['Lat'].max()
        #LonIni = self['Lon'].min()
        #LonFin = self['Lon'].max()
        pass
    else:
        LatIni = limits['LatIni']
        LatFin = limits['LatFin']
        LonIni = limits['LonIni']
        LonFin = limits['LonFin']

    # ----
    cdict = {'red': ((0.0, 0.0, 0.0),
             (0.5, 0.879, 0.879),
             (1.0, 0.0, 0.0)),
     'green': ((0.0, 0.281, 0.281),
               (0.5, 0.418, 0.418),
               (1.0, 0.281, 0.281)),
     'blue': ((0.0, 0.195, 0.195),
              (0.5, 0.184, 0.184),
              (1.0, 0.281, 0.281))}
    um_sym_cmap = matplotlib.colors.LinearSegmentedColormap('um_sym_colormap',cdict,256)

    cdict = {'red': ((0.0, 0.879, 0.879),
             (1.0, 0.0, 0.0)),
     'green': ((0.0, 0.418, 0.418),
               (1.0, 0.281, 1.0)),
     'blue': ((0.0, 0.184, 0.184),
              (1.0, 0.281, 0.281))}
    um_cmap = matplotlib.colors.LinearSegmentedColormap('um_colormap',cdict,256)

    #cdict = {'red': ((0.0, 0.879, 0.879),
    #         (0.5, 0.0, 0.0),
    #         (1.0, 0.879, 0.879)),
    # 'green': ((0.0, 0.418, 0.418),
    #           (0.5, 0.281, 0.281),
    #           (1.0, 0.418, 0.418)),
    # 'blue': ((0.0, 0.184, 0.184),
    #          (0.5, 0.195, 0.195),
    #          (1.0, 0.184, 0.184))}
    #um_cmap = matplotlib.colors.LinearSegmentedColormap('um_colormap',cdict,256)

    cdict = {'red': ((0.0, 0.0, 0.0),
             (0.5, 1.0, 1.0),
             (1.0, 0.0, 0.0)),
     'green': ((0.0, 0.0, 0.0),
               (0.5, 1.0, 1.0),
               (1.0, 0.0, 0.0)),
     'blue': ((0.0, 0.0, 0.0),
              (0.5, 1.0, 1.0),
              (1.0, 0.0, 0.0))}
    bw_cmap = matplotlib.colors.LinearSegmentedColormap('bw_colormap',cdict,256)

    cdict = {'red': ((0.0, 0.45, 0.45),
             (0.5, 0.95, 0.95),
             (1.0, 0.45, 0.45)),
     'green': ((0.0, 0.45, 0.45),
               (0.5, 0.95, 0.95),
               (1.0, 0.45, 0.45)),
     'blue': ((0.0, 0.45, 0.45),
              (0.5, .95, 0.95),
              (1.0, 0.45, 0.45))}
    grey_cmap = matplotlib.colors.LinearSegmentedColormap('grey_colormap',cdict,256)


    # ----

    parallels = np.arange(-5,20.1,5)
    meridians = np.arange(300,340,10)

    margin=0.08
    left=margin
    bottom=margin
    height_eof = (1-4*margin)*.44
    width_eof =  (1-3*margin)*.5

    height_pc = (1-4*margin)*.28
    width_pc = (1-2*margin)*1
    # ----
    eof_amp=(eof.real**2+eof.imag**2)**0.5
    eof_phase=np.arctan2(eof.imag,eof.real)
    pc_amp = (np.real(pc)**2+np.imag(pc)**2)**0.5
    pc_phase = np.arctan2(np.imag(pc),np.real(pc))

    fig = pylab.figure(figsize=(14,10.5), dpi=300)
    cumvarfrac = None
    if cumvarfrac != None:
        title = "Mode: %i (%i%%) (cumulative %i%%)" % (nmode,varfrac*1e2,cumvarfrac*1e2)
    else:
        title = "Mode: %i (%i%%)" % (nmode,varfrac*1e2)

    #if 'halfpower_period' in self:
    #if 'prefilter' in self.metadata:
    #    if type(self.halfpower_period) == str:
    #        halfpower_period = self.halfpower_period
    #    else:
    #        halfpower_period = round(self.halfpower_period)
    #    title = "%s (%s half power:%s days)" % (title,self.metadata['prefilter']['type'], halfpower_period)
    fig.text(.5, .95, title, horizontalalignment='center',fontsize=16)
    #
    pylab.axes([left, bottom + 2*height_pc + 2*margin, width_eof, height_eof])
    map = Basemap(projection='merc', lat_ts=0, llcrnrlon=LonIni,
            llcrnrlat=LatIni, urcrnrlon=LonFin, urcrnrlat=LatFin,
            resolution='l', area_thresh=1000.)
    X, Y = map(*pylab.meshgrid(data['lon'], data['lat']))
    map.contourf(X, Y, eof_amp*1e2)
    pylab.title("CEOF amplitude")
    cbar = pylab.colorbar()
    cbar.set_label('[cm]')
    map.drawcoastlines()
    map.fillcontinents(color='0.0')
    map.drawparallels(parallels,labels=[1,0,0,1])
    map.drawmeridians(meridians,labels=[1,0,0,1])

    pylab.axes([left+width_eof+margin, bottom + 2*height_pc + 2*margin, width_eof, height_eof])
    map = Basemap(projection='merc',lat_ts=0,llcrnrlon=LonIni,llcrnrlat=LatIni, urcrnrlon=LonFin, urcrnrlat=LatFin,resolution='l',area_thresh=1000.)
    X, Y = map(*pylab.meshgrid(data['lon'],data['lat']))
    V=[-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180]
    #V = range(-180,181,20)
    #import matplotlib.cm as cm
    #map.pcolor(X[0,:],Y[:,0],eof_phase*180/np.pi,cmap=um_sym_cmap)
    #map.contourf(X,Y,eof_phase*180/np.pi,V,cmap=um_sym_cmap)
    #from scipy.stats import scoreatpercentile
    #scoreatpercentile(x.flatten(),15)

    from numpy import ma
    #ind_sig = eof_amp<0.01
    #eof_phase_deg = eof_phase*180/np.pi
    map.contourf(X, Y, eof_phase*180/np.pi, V, cmap=grey_cmap)
    map.contourf(X,Y,ma.masked_array(eof_phase*180/np.pi, mask=eof_amp<0.01),V,cmap=um_sym_cmap)

    cbar = pylab.colorbar()
    cbar.set_label('[degrees]')
    pylab.title("CEOF phase")
    map.drawcoastlines()
    map.fillcontinents(color='0.0')
    map.drawparallels(parallels,labels=[1,0,0,1])
    map.drawmeridians(meridians,labels=[1,0,0,1])
    # ----
    #pylab.subplot(2,2,2)
    pylab.axes([left, bottom+margin+height_pc, width_pc, height_pc])
    pylab.plot_date(pylab.date2num(data['datetime']),pc_amp,'-')
    fig.autofmt_xdate()
    pylab.title("PC amplitude")
    pylab.ylabel('[dimensionless]')
    pylab.grid()
    # ----
    #pylab.subplot(2,2,4)
    pylab.axes([left, bottom, width_pc, height_pc])
    pylab.plot_date(pylab.date2num(data['datetime']),pc_phase*180/np.pi,'.')
    fig.autofmt_xdate()
    v = pylab.axis()
    pylab.axis((v[0],v[1],-181,181))
    pylab.title("PC phase")
    pylab.ylabel('[degrees]')
    pylab.grid()
    # ----
    #pylab.subplot(2,2,4)
    #pylab.axes([left + margin + width_l, bottom, width_r, height_r])
    #pylab.plot(np.absolute(scipy.fftpack.fft(pc_amp))[1:pc_amp.shape[0]/2])
    #pylab.title("PC FFT")
    #pylab.grid()
    # ----
    #pylab.show()
    print "Saving figure %s" % filename
    #fig.savefig(filename)
    pylab.savefig(filename)
    pylab.close()
