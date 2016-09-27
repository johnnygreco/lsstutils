from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
DATADIR = '/Volumes/HSC-20160523/'

__all__ = ['DATADIR', 'BetterButler']

class BetterButler(object):
    """
    Wrapper for interatcting with HSC pipeline outputs.

    Parameters
    ----------
    tract : int
        HSC tract number.
    patch : string
        HSC patch. e.g., '5,7'.
    band : string, optional
        The photometric band of the observation ('G', 'R', 'I', 'Z', or 'Y').
    butler : Butler object, optional
        If None, will create a butler at initialization.
    dataDIR : string, optional
        Data directory. 
    """

    def __init__(self, tract, patch, band='I', butler=None, dataDIR=DATADIR):

        if butler is None:
            import lsst.daf.persistence
            self._butler = lsst.daf.persistence.Butler(dataDIR)
        else:
            self._butler = butler

        band = band.upper()
        self.dataID = {'tract':tract, 'patch':patch, 'filter':'HSC-'+band}

        self._fn = self._butler.get('deepCoadd_calexp_filename', self.dataID)[0]
        self._cat = None
        self._calexp = None
        self._calib = None
        self._wcs = None
        self._maskedImg = None

    @property
    def butler(self):
        """
        The Butler.
        """
        return self._butler

    @property
    def cat(self):
        """
        The main coadd catalog for the given dataID.
        """
        if self._cat is None:
            self._cat = self.butler.get('deepCoadd_meas', self.dataID, immediate=True)
        return self._cat

    @property
    def calexp(self):
        """
        The main calibrated exposure object. 
        """
        if self._calexp is None:
            self._calexp = self.butler.get('deepCoadd_calexp', self.dataID, immediate=True)
        return self._calexp

    @property
    def calib(self):
        """
        Calibration object.
        """
        if self._calib is None:
            self._calib = self.calexp.getCalib()
        return self._calib
    
    @property
    def wcs(self):
        """
        World Coordinate System object for this dataID.
        """
        if self._wcs is None:
            self._wcs = self.calexp.getWcs()
        return self._wcs

    @property
    def maskedImg(self):
        """
        The masked image object.
        """
        return self.calexp.getMaskedImage()

    def get_fn(self):
        """
        Return the fits file name for this exposure.
        """
        return self._fn

    def get_pixscale(self):
        """
        Return the pixel scale.
        """
        return self.wcs.pixelScale().asArcseconds()

    def get_zptmag(self):
        """
        Return the zero point magnitude. 
        """
        return 2.5*np.log10(self.calib.getFluxMag0()[0])

    def get_psf(self):
        """
        Return the PSF as a 2D numpy array.
        """
        return self.calexp.getPsf().computeImage().getArray().copy()

    def get_img(self):
        """
        Return image as a 2D numpy array.
        """
        return self.maskedImg.getImage().getArray().copy()

    def get_mask(self):
        """
        Return complete pipeline mask as a 2D numpy array.
        """
        return self.maskedImg.getMask().getArray().copy()

    def get_badmask(self):
        """
        Return a bad pixel mask as 2D numpy array.
        """
        mask = self.maskedImg.getMask()
        detected = mask.getPlaneBitMask('DETECTED')
        bad = mask.getArray().copy()
        bad[bad==detected] = 0
        return bad

    def get_detmask(self):
        """
        Return a detected pixel mask as a 2D numpy array.
        """
        mask = self.maskedImg.getMask()
        detected = mask.getPlaneBitMask('DETECTED')
        det = mask.getArray().copy()
        det[det!=detected] = 0
        return det

    def get_var(self):
        """
        Return sigma image as a 2D numpy array.
        """
        return self.maskedImg.getVariance().getArray().copy()

    def write_fits(self, outfile):
        """
        Write the calibrated exposure to a fits file. 

        Parameters
        ---------- 
        outfile : sting
            The output file name. If you want the file to
            be saved in a directory other than the current
            one, give the full or relative path. 

        Notes
        -----
        The output fits file will be a multi-extension 
        file with frames image, mask, & variance. There
        will also be a header and 14 bin tables.
        """
        self.calexp.writeFits(outfile)
