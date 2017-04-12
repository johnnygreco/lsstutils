from __future__ import division, print_function

import os
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.daf.persistence
import lsst.afw.coord as afwCoord
import lsst.daf.base
from astropy import units as u
from .utils import tracts_n_patches, sky_cone

ROOT = '/tigress/HSC/HSC/rerun/production-20160523'

__all__ = ['make_stamp']

def make_stamp(ra, dec, radius, band='i', skymap=None, butler=None, 
               pixscale=0.168, root=ROOT):
    """
    Generate HSC cutout image.
    
    Parameters
    ----------
    ra, dec: float
        Center of cutout.
    radius: astropy Quantity, float, or int
        Angular radius of cone. Must be in degrees
        if not a Quantity object.
    band: str, optional
        HSC photometric band.
    skymap: lsst.skymap.ringsSkyMap.RingsSkyMap, optional
        Skymap object (pass one to function if making many cutouts).
    butler: lsst.daf.persistence.butler.Butler, optional
        HSC data butler (pass one to function if making many cutouts).
    pixscale: float, optional
        Image pixel scale in arcsec/pixel. 
    
    Returns
    -------
    stamp: lsst.afw.image.ExposureF
        Cutout exposure object. 
    """
    
    if butler is None:
        butler = lsst.daf.persistence.Butler(root)
    if skymap is None:
        skymap = butler.get('deepCoadd_skyMap', immediate=True)
        
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
        data_id = {'tract':t, 'patch':p, 'filter':'HSC-'+band.upper()}
        fn = butler.get('deepCoadd_calexp_filename', data_id)[0]
        if os.path.isfile(fn):
            img = butler.get('deepCoadd_calexp', data_id)
            images.append(img)

    ########################################
    # Get cutouts from each patch
    ########################################
    
    cutouts = []
    idx = []
    bbox_sizes = []
    size = int(radius.to('arcsec').value/pixscale)
    coord = afwCoord.Coord(ra*afwGeom.degrees, dec*afwGeom.degrees)
    
    for count, img in enumerate(images):
        wcs = img.getWcs()
        pix = wcs.skyToPixel(coord)
        pix = afwGeom.Point2I(pix)
        bbox = afwGeom.Box2I(pix, pix)
        bbox.grow(size)
        x0, y0 = bbox.getBegin()
        bbox.clip(img.getBBox(afwImage.PARENT))
        xnew, ynew = bbox.getBeginX()-x0, bbox.getBeginY()-y0
        idx.append([xnew, xnew+bbox.getWidth(), ynew, ynew+bbox.getHeight()])
        bbox_sizes.append(bbox.getWidth() * bbox.getHeight())
        cut = img.Factory(img, bbox, afwImage.PARENT)
        cutouts.append(cut)
        if count==0:
            subwcs = cut.getWcs()
            crpix_1, crpix_2 = subwcs.skyToPixel(coord)
            crpix_1 -= x0 - 1
            crpix_2 -= y0 - 1
            cdmat = wcs.getCDMatrix()      

    ########################################
    # Build new WCS for cutout
    ########################################
    
    md = lsst.daf.base.PropertyList()
    md.add('CRVAL1', ra)
    md.add('CRVAL2', dec)
    md.add('CRPIX1', crpix_1)
    md.add('CRPIX2', crpix_2)
    md.add('CTYPE1', 'RA---TAN')
    md.add('CTYPE2', 'DEC--TAN')
    md.add('CD1_1', cdmat[0, 0])
    md.add('CD2_1', cdmat[1, 0])
    md.add('CD1_2', cdmat[0, 1])
    md.add('CD2_2', cdmat[1, 1])
    md.add('RADESYS', 'ICRS')
    stamp_wcs = afwImage.makeWcs(md)
    
    ########################################
    # Stitch cutouts together with the
    # largest bboxes inserted last
    ########################################

    shape = (size*2+1, size*2+1)
    stamp_bbox = afwGeom.BoxI(afwGeom.Point2I(0,0), afwGeom.Extent2I(*shape))
    stamp = afwImage.MaskedImageF(stamp_bbox)
    for i in np.argsort(bbox_sizes):
        if cutouts[i].getDimensions() == afwGeom.Extent2I(*shape):
            return cutouts[i]
        mi = cutouts[i].getMaskedImage()
        stamp[idx[i][0]: idx[i][1], idx[i][2]: idx[i][3]] = mi
    stamp = afwImage.ExposureF(stamp, stamp_wcs)
    
    return stamp
