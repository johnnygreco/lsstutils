"""
Image processing tools that use the lsst stack. 
"""
from __future__ import division, print_function

import lsst.afw.math as afwMath

__all__ = ['smooth_gauss']

def smooth_gauss(masked_image, sigma, nsigma=7.0):
    """
    Smooth image with a Gaussian kernel. 

    Parameters
    ----------
    masked_image : lsst.afw.image.imageLib.MaskedImageF
        Masked image object to be smoothed
    sigma : float
        Standard deviation of Gaussian
    nsigma : float, optional
        Number of sigma for kernel width

    Returns
    -------
    convolved_image : lsst.afw.image.imageLib.MaskedImageF
        The convolved masked image
    """
    width = (int(sigma*nsigma + 0.5) // 2)*2 + 1 # make sure it is odd
    gauss_func = afwMath.GaussianFunction1D(sigma)
    gauss_kern = afwMath.SeparableKernel(width, width, gauss_func, gauss_func)
    convolved_image = masked_image.Factory(masked_image.getBBox())
    afwMath.convolve(convolved_image, masked_image, gauss_kern, 
                     afwMath.ConvolutionControl())
    return convolved_image
