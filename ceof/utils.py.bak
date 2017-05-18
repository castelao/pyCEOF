import numpy as np

def scaleEOF(pcs, eofs, scaletype):
    """ Scale the EOFS and PCS preserving the mode

        If scaled by the pcmedian:
            pc = pc/median(pc)
            eof = eof*median(pc)
        This is the scale that I most use. On this case, the EOF structure
          has that magnitude at least half of the timeseries.
    """

    assert pcs.ndim == 2
    assert eofs.ndim == 2, "Expected a 2D EOF array."
    assert pcs.shape[-1] == eofs.shape[-1]

    nmodes = pcs.shape[-1]

    if scaletype == 'pc_std':
        print "Normalizing by the pc_std"
        for n in range(nmodes):
            fac = (np.absolute(pcs[:,n])).std()
            pcs[:,n] = pcs[:,n]/fac
            eofs[:,n] = eofs[:,n]*fac

    elif scaletype == 'pc_median':
        print "Normalizing by the pc_median"
        for n in range(nmodes):
            fac = np.median((np.absolute(pcs[:,n])))
            pcs[:,n] = pcs[:,n]/fac
            eofs[:,n] = eofs[:,n]*fac

    elif scaletype == 'pc_max':
        print "Normalizing by the pc_max"
        for n in range(nmodes):
            fac = (np.absolute(pcs[:,n])).max()
            pcs[:,n] = pcs[:,n]/fac
            eofs[:,n] = eofs[:,n]*fac

    elif scaletype == 'eof_std':
        print "Normalizing by the eof_std"
        for n in range(nmodes):
            fac = (np.absolute(eofs[:,n])).std()
            eofs[:,n] = eofs[:,n]/fac
            pcs[:,n] = pcs[:,n]*fac

    elif scaletype == 'eof_max':
        print "Normalizing by the eof_max"
        for n in range(nmodes):
            fac = (np.absolute(eofs[:,n])).max()
            eofs[:,n] = eofs[:,n]/fac
            pcs[:,n] = pcs[:,n]*fac

    else:
        print "Don't understand the normalization config: %s" % scaletype
        return

    return pcs, eofs


def ceof_reconstruct(eofs, pcs, nmodes=None):
    """ Reconstruct the dataset from the sum of eofs*pcs

        If modes is 'all', uses all modes to reconstruct, otherwise it
          is expected an integer, and will be the number of modes used.
    """
    assert eofs.shape[-1] == pcs.shape[-1], \
            "Last dimension of eofs must be equal to last dimension of pcs."
    assert (nmodes == None) or (type(nmodes) == int), \
            "Valid values for modes are: 'all' or an integer."

    if nmodes is None:
        nmodes = range(pcs.shape[-1])
    else:
        nmodes = range(nmodes)

    print "Reconstructing from EOF using the modes: %s" % nmodes
    T = pcs.shape[-1]

    eof_amp = (eofs.real**2 + eofs.imag**2)**0.5
    eof_phase = numpy.arctan2(eofs.imag, eofs.real)
    pc_amp = (numpy.real(pcs)**2 + numpy.imag(pcs)**2)**0.5
    pc_phase = numpy.arctan2(numpy.imag(pcs), numpy.real(pcs))

    data = np.zeros([T] + list(eofs.shape[:-1]))

    assert (eofs.ndim == 3), "Sorry, I only able to handle a 3D EOF array."

    for t in range(T):
        for n in modes:
            data[t] = data[t] + eof_amp[:,:,n] * pc_amp[t,n] * \
                    numpy.cos(eof_phase[:,:,n] + pc_phase[t,n])

    return data


def gridto2D(self, var, ind=None):
    """
    """
    I, J, K = self.data[var].shape

    if ind == None:
        ind = numpy.ones((J,K))==1

    N = ((numpy.ones(ind.shape)[ind]).sum())

    self.data2D = {}
    for k in ['lat','lon']:
        self.data2D[k] = ma.masked_all(N, dtype=self.data[k].dtype)

    self.data2D['grid_index'] = ma.masked_all((N, 2))

    for k in [var]:
        self.data2D[k] = ma.masked_all((I,N), dtype=self.data[k].dtype)

    n=-1
    for j in range(J):
        for k in range(K):
            if ind[j, k]:
                n += 1
                self.data2D['grid_index'][n] = numpy.array([j,k])
                self.data2D['lat'][n] = self.data['Lat'][j,k]
                self.data2D['lon'][n] = self.data['Lon'][j,k]
                self.data2D[var][:,n] = self.data[var][:,j,k]
    return
