from __future__ import division, print_function

import os
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.display.rgb as afwRgb
import lsst.daf.base
import lsst.daf.persistence
from astropy import units as u
from .utils import get_psf, tracts_n_patches, sky_cone

ROOT = '/tigress/HSC/HSC/rerun/production-20160523'

__all__ = ['make_stamp', 'make_rgb_image']

def make_stamp(ra, dec, radius, band='i', skymap=None, butler=None, 
               pixscale=0.168, root=ROOT, return_psf=False, 
               coadd_label='deepCoadd_calexp'):
    """
    Generate HSC cutout image.
    
    Parameters
    ----------
    ra, dec: float
        Center of cutout.
    radius: astropy Quantity, float, or int
        Angular radius of cone. Must be in arcsec
        if not a Quantity object.
    band: str, optional
        HSC photometric band.
    skymap: lsst.skymap.ringsSkyMap.RingsSkyMap, optional
        Skymap object (pass one to function if making many cutouts).
    butler: lsst.daf.persistence.butler.Butler, optional
        HSC data butler (pass one to function if making many cutouts).
    pixscale: float, optional
        Image pixel scale in arcsec/pixel. 
    return_psf: bool, optional
        If True, return image PSF.

    Returns
    -------
    stamp: lsst.afw.image.ExposureF
        Cutout exposure object. 
    """
    
    if butler is None:
        butler = lsst.daf.persistence.Butler(root)
    if skymap is None:
        skymap = butler.get('deepCoadd_skyMap', immediate=True)
    if type(radius) != u.Quantity:
        radius *= u.arcsec

    size = int(radius.to('arcsec').value/pixscale)
    coord = afwCoord.IcrsCoord(ra*afwGeom.degrees, dec*afwGeom.degrees)
    stamp_shape = (size*2+1, size*2+1)
        
    ########################################
    # Generate ra/dec list & get patches
    ########################################
    
    radec_list = np.array(sky_cone(ra, dec, radius)).T    
    patches, _ = tracts_n_patches(radec_list, skymap)

    ########################################
    # Get exposures for all patches
    ########################################
    
    images = []
    for t, p in patches:
        data_id = {'tract':t, 'patch':p.decode(), 'filter':'HSC-'+band.upper()}
        if butler.datasetExists(coadd_label, data_id):
            img = butler.get(coadd_label, data_id, immediate=True)
            images.append(img)

    if len(images)==0:
        print('***** No data at {:.5f} {:.5f} *****'.format(ra, dec))
        return None

    ########################################
    # Get cutouts from each patch
    ########################################
    
    cutouts = []
    idx = []
    bbox_sizes = []
    bbox_origins = []
    
    for img in images:
        wcs = img.getWcs()
        pix = wcs.skyToPixel(coord)
        pix = afwGeom.Point2I(pix)
        bbox = afwGeom.Box2I(pix, pix)
        bbox.grow(size)
        x0, y0 = bbox.getBegin()
        bbox_origins.append([x0, y0])
        bbox.clip(img.getBBox(afwImage.PARENT))
        xnew, ynew = bbox.getBeginX()-x0, bbox.getBeginY()-y0
        idx.append([xnew, xnew+bbox.getWidth(), ynew, ynew+bbox.getHeight()])
        bbox_sizes.append(bbox.getWidth() * bbox.getHeight())
        cut = img.Factory(img, bbox, afwImage.PARENT)
        cutouts.append(cut)

    ########################################
    # Stitch cutouts together with the
    # largest bboxes inserted last
    ########################################

    stamp_bbox = afwGeom.BoxI(afwGeom.Point2I(0,0), afwGeom.Extent2I(*stamp_shape))
    stamp = afwImage.MaskedImageF(stamp_bbox)
    bbox_sorted_ind = np.argsort(bbox_sizes)
    for i in bbox_sorted_ind:
        mi = cutouts[i].getMaskedImage()
        stamp[idx[i][0]: idx[i][1], idx[i][2]: idx[i][3]] = mi

    ########################################
    # Build new WCS for cutout
    ########################################

    largest_cutout = cutouts[bbox_sorted_ind[-1]]
    subwcs = largest_cutout.getWcs()
    crpix_1, crpix_2 = subwcs.skyToPixel(coord)
    crpix_1 -= bbox_origins[bbox_sorted_ind[-1]][0] 
    crpix_2 -= bbox_origins[bbox_sorted_ind[-1]][1]
    cdmat = wcs.getCDMatrix()      
    
    md = lsst.daf.base.PropertyList()
    md.add('CRVAL1', ra)
    md.add('CRVAL2', dec)
    md.add('CRPIX1', crpix_1 + 1)
    md.add('CRPIX2', crpix_2 + 1)
    md.add('CTYPE1', 'RA---TAN')
    md.add('CTYPE2', 'DEC--TAN')
    md.add('CD1_1', cdmat[0, 0])
    md.add('CD2_1', cdmat[1, 0])
    md.add('CD1_2', cdmat[0, 1])
    md.add('CD2_2', cdmat[1, 1])
    md.add('RADESYS', 'ICRS')
    stamp_wcs = afwImage.makeWcs(md)

    stamp = afwImage.ExposureF(stamp, stamp_wcs)
    
    return (stamp, get_psf(largest_cutout, coord)) if return_psf else stamp


def make_rgb_image(ra, dec, radius, butler=None, skymap=None, 
                   rgb='irg', Q=8, dataRange=0.6, root=ROOT, img_size=None, 
                   return_wcs=False):

    if butler is None:
        butler = lsst.daf.persistence.Butler(root)
    if skymap is None:
        skymap = butler.get('deepCoadd_skyMap', immediate=True)
    if type(radius)==float or type(radius)==int:
        radius *= u.arcsec

    colors = []
    for band in rgb:
        stamp = make_stamp(ra, dec, radius, band=band, 
                           butler=butler, skymap=skymap)
        if stamp:
            colors.append(stamp.getMaskedImage())

        if return_wcs and band==rgb[0]:
            wcs = stamp.getWcs()

    if len(colors)==0:
        return None

    rgb_kws = {'Q': Q, 'dataRange': dataRange}
    if img_size is not None:
        rgb_kws['xSize'] = img_size
    img = afwRgb.makeRGB(*colors, **rgb_kws)

    return (img, wcs) if return_wcs else img
