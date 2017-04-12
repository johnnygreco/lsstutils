"""
Generate HSC cutout images.
Code adapted from Song Huang:
https://github.com/dr-guangtou
"""

import os, fcntl
import warnings, copy
import numpy as np
from astropy import wcs as apWcs
from astropy.io import fits
import lsst.daf.persistence as dafPersist
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable

COM = '#' * 100
SEP = '-' * 100
WAR = '!' * 100

def getTractPatchList(matches):
    """Get the list of Tract, Patch."""
    tract = []
    patch = []

    for match in matches:
        tractInfo, patchInfo = match
        tractId = tractInfo.getId()
        for patchItem in patchInfo:
            tract.append(tractId)
            patch.append("%d,%d" % patchItem.getIndex())

    return tract, patch


def saveImageArr(arr, header, name, clobber=True):
    """
    Just save an array to a fits file.

    Parameters:
    """
    hduImg = fits.PrimaryHDU(arr, header=header)
    hduList = fits.HDUList([hduImg])
    hduList.writeto(name, clobber=clobber)
    hduList.close()


def getCircleRaDec(ra, dec, size):
    """
    Get the (RA, DEC) that describe a circle.
    Region around the central input coordinate
    """
    # Convert the size from pixel unit to degress
    sizeDegree = (size * 0.168) / 3600.0
    # representative set of polar angles
    angles = np.array([0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0])
    phi = np.array(angles * np.pi / 180.0)
    raList = ra + sizeDegree * np.cos(phi)
    decList = dec + sizeDegree * np.sin(phi)
    return raList, decList


def getCoaddPsfImage(calExp, coord):
    """Get the coadd PSF image."""
    # Get the WCS information
    wcs = calExp.getWcs()
    # The X,Y coordinate of the image center
    coordXY = wcs.skyToPixel(coord)
    # Get the PSF object for the exposure
    psf = calExp.getPsf()
    try:
        psfImg = psf.computeImage(coordXY)
        return psfImg
    except Exception:
        warnings.warn("### Can not compute PSF Image !!!")
        return None


def getCoaddMskPlane(calExp, bitmask):
    """Get the mask plane of the coadd."""
    # Get the mask image
    mskImg = calExp.getMaskedImage().getMask()
    newMsk = copy.deepcopy(mskImg)
    try:
        # Extract specific plane from it
        newMsk &= newMsk.getPlaneBitMask(bitmask)
    except Exception:
        mskImg.printMaskPlanes()
        newMsk = None

    return newMsk


def getCoaddBadMsk(calExp, pipeNew=False):
    """Get the BAD mask plane."""
    mskImg = calExp.getMaskedImage().getMask()

    badMsk = copy.deepcopy(mskImg)
    # Clear the "DETECTED" plane
    badMsk.removeAndClearMaskPlane('DETECTED', True)
    try:
        # Clear the "EDGE" plane
        badMsk.removeAndClearMaskPlane('EDGE', True)
    except Exception:
        pass
    try:
        # Clear the "DETECTED_NEGATIVE" plane
        badMsk.removeAndClearMaskPlane('DETECTED_NEGATIVE', True)
    except Exception:
        pass
    try:
        # Clear the "CLIPPED" plane
        badMsk.removeAndClearMaskPlane('CLIPPED', True)
    except Exception:
        pass
    try:
        # Clear the "CROSSTALK" plane
        badMsk.removeAndClearMaskPlane('CROSSTALK', True)
    except Exception:
        pass
    if pipeNew:
        try:
            # Clear the "NOT_DEBLENDED" plane
            badMsk.removeAndClearMaskPlane('NOT_DEBLENDED', True)
        except Exception:
            pass
        try:
            # Clear the "BRIGHT_OBJECT" plane
            badMsk.removeAndClearMaskPlane('BRIGHT_OBJECT', True)
        except Exception:
            pass

    return badMsk



def coaddImageCutFull(root, ra, dec, size, saveSrc=True, savePsf=True,
                      filt='HSC-I', prefix='hsc_coadd_cutout', verbose=True,
                      extraField1=None, extraValue1=None, butler=None,
                      visual=True, imgOnly=False):
    """Get the cutout around a location."""
    coaddData = "deepCoadd_calexp"
    pipeNew = True

    # Get the SkyMap of the database
    if butler is None:
        try:
            butler = dafPersist.Butler(root)
            if verbose:
                print SEP
                print "## Load in the Butler"
        except Exception:
            print WAR
            print '## Can not load the correct Butler!'
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    # (Ra, Dec) Pair for the center
    raDec = afwCoord.Coord(ra*afwGeom.degrees, dec*afwGeom.degrees)
    # [Ra, Dec] list
    raList, decList = getCircleRaDec(ra, dec, size)
    points = map(lambda x, y: afwGeom.Point2D(x, y), raList, decList)
    raDecList = map(lambda x: afwCoord.IcrsCoord(x), points)

    # Expected size and center position
    dimExpect = int(2 * size + 1)
    cenExpect = (dimExpect/2.0, dimExpect/2.0)
    sizeExpect = int(dimExpect ** 2)
    # Get the half size of the image in degree
    sizeDegree = size * 0.168 / 3600.0

    # Verbose
    if verbose:
        print SEP
        print " Input Ra, Dec: %10.5f, %10.5f" % (ra, dec)
        print " Cutout size is expected to be %d x %d" % (dimExpect, dimExpect)

    # Create empty arrays
    imgEmpty = np.empty((dimExpect, dimExpect), dtype="float")
    imgEmpty.fill(np.nan)
    if not imgOnly:
        mskEmpty = np.empty((dimExpect, dimExpect), dtype="uint8")
        varEmpty = np.empty((dimExpect, dimExpect), dtype="float")
        detEmpty = np.empty((dimExpect, dimExpect), dtype="float")
        mskEmpty.fill(np.nan)
        varEmpty.fill(np.nan)
        detEmpty.fill(np.nan)

    # Figure out the area we want, and read the data.
    # For coadds the WCS is the same in all bands,
    # but the code handles the general case
    # Start by finding the tract and patch
    matches = skyMap.findTractPatchList(raDecList)
    tractList, patchList = getTractPatchList(matches)
    nPatch = len(patchList)
    if verbose:
        print "### Will deal with %d patches" % nPatch
    # Prefix of the output file
    outPre = prefix + '_' + filt + '_full'

    newX = []
    newY = []
    boxX = []
    boxY = []
    boxSize = []
    #
    trList = []
    paList = []
    zpList = []
    #
    imgArr = []
    mskArr = []
    varArr = []
    detArr = []
    psfArr = []
    #
    srcArr = []
    refArr = []
    forceArr = []

    # Go through all these images
    for j in range(nPatch):
        # Tract, patch
        tract, patch = tractList[j], patchList[j]
        print SEP
        print "### Dealing with %d - %s" % (tract, patch)
        print SEP
        # Check if the coordinate is available in all three bands.
        try:
            # Get the coadded exposure
            coadd = butler.get(coaddData, tract=tract,
                               patch=patch, filter=filt,
                               immediate=True)
        except Exception, errMsg:
            print WAR
            print " No data is available in %d - %s" % (tract, patch)
            print "#########################################################"
            print WAR
        else:
            # Get the WCS information
            wcs = coadd.getWcs()
            # Check if cdMatrix has been assigned
            cdExist = 'cdMatrix' in locals()
            if not cdExist:
                # Get the CD Matrix of the WCS
                cdMatrix = wcs.getCDMatrix()
                # Get the pixel size in arcsec
                pixScale = wcs.pixelScale().asDegrees() * 3600.0
            # Convert the central coordinate from Ra,Dec to pixel unit
            pixel = wcs.skyToPixel(raDec)
            pixel = afwGeom.Point2I(pixel)
            # Define the bounding box for the central pixel
            bbox = afwGeom.Box2I(pixel, pixel)
            # Grow the bounding box to the desired size
            bbox.grow(int(size))
            xOri, yOri = bbox.getBegin()
            # Compare to the coadd image, and clip
            bbox.clip(coadd.getBBox(afwImage.PARENT))
            # Get the masked image
            try:
                subImage = afwImage.ExposureF(coadd, bbox,
                                              afwImage.PARENT)
            except Exception:
                print WAR
                print '### SOMETHING IS WRONG WITH THIS BOUNDING BOX !!'
                print "    %d -- %s -- %s " % (tract, patch, filt)
                print "    Bounding Box Size: %d" % (bbox.getWidth() *
                                                     bbox.getHeight())
            else:
                # Extract the image array
                imgArr.append(subImage.getMaskedImage().getImage().getArray())

                if not imgOnly:
                    # Extract the detect mask array
                    mskDet = getCoaddMskPlane(subImage, 'DETECTED')
                    detArr.append(mskDet.getArray())
                    # Extract the variance array
                    imgVar = subImage.getMaskedImage().getVariance().getArray()
                    varArr.append(imgVar)

                    # Extract the bad mask array
                    mskBad = getCoaddBadMsk(subImage, pipeNew=pipeNew)
                    mskArr.append(mskBad.getArray())

                # Save the width of the BBox
                boxX.append(bbox.getWidth())
                # Save the heigth of the BBox
                boxY.append(bbox.getHeight())
                # Save the size of the BBox in unit of pixels
                boxSize.append(bbox.getWidth() * bbox.getHeight())
                # New X, Y origin coordinates
                newX.append(bbox.getBeginX() - xOri)
                newY.append(bbox.getBeginY() - yOri)
                # Tract, Patch
                trList.append(tract)
                paList.append(patch)
                # Photometric zeropoint
                zpList.append(2.5 * np.log10(
                              coadd.getCalib().getFluxMag0()[0]))
                # If necessary, save the psf images
                if savePsf and (not imgOnly):
                    psfImg = getCoaddPsfImage(coadd, raDec)
                    psfArr.append(psfImg)
                # Get the new (X,Y) coordinate of the galaxy center
                newCenExist = 'newCenX' in locals() and 'newCenY' in locals()
                if not newCenExist:
                    subWcs = subImage.getWcs()
                    newCenX, newCenY = subWcs.skyToPixel(raDec)
                    newCenX = newCenX - xOri
                    newCenY = newCenY - yOri
    # Number of returned images
    nReturn = len(newX)
    if nReturn > 0:
        print "### Return %d Useful Images" % nReturn
        # Sort the returned images according to the size of their BBox
        indSize = np.argsort(boxSize)

        # Go through the returned images, put them in the cutout region
        for n in range(nReturn):
            ind = indSize[n]
            # Put in the image array
            imgUse = imgArr[ind]
            imgEmpty[newY[ind]:(newY[ind] + boxY[ind]),
                     newX[ind]:(newX[ind] + boxX[ind])] = imgUse[:, :]
            # Put in the mask array
            if not imgOnly:
                mskUse = mskArr[ind]
                mskEmpty[newY[ind]:(newY[ind] + boxY[ind]),
                         newX[ind]:(newX[ind] + boxX[ind])] = mskUse[:, :]
                # Put in the variance array
                varUse = varArr[ind]
                varEmpty[newY[ind]:(newY[ind] + boxY[ind]),
                         newX[ind]:(newX[ind] + boxX[ind])] = varUse[:, :]
                # Convert it into sigma array
                sigEmpty = np.sqrt(varEmpty)
                # Put in the detection mask array
                detUse = detArr[ind]
                detEmpty[newY[ind]:(newY[ind] + boxY[ind]),
                         newX[ind]:(newX[ind] + boxX[ind])] = detUse[:, :]
            if n is (nReturn - 1):
                # This is the largest available sub-image
                phoZp = zpList[ind]
                # Save the psf image if necessary
                if not imgOnly:
                    if savePsf:
                        psfOut = outPre + '_psf.fits'
                        if psfArr[ind] is not None:
                            psfUse = psfArr[ind]
                            psfUse.writeFits(psfOut)
                            noPsf = False
                        else:
                            warnings.warn("## Can not compute PSF image !!")
                            noPsf = True
                else:
                    noPsf = True
        # See if all the cutout region is covered by data
        nanPix = np.sum(np.isnan(imgEmpty))
        if nanPix < (sizeExpect * 0.1):
            cutFull = True
            if verbose:
                print "## > 90% of the cutout region is covered!"
        else:
            cutFull = False
            if verbose:
                print "## There are still %d NaN pixels!" % nanPix
        if not imgOnly:
            # For mask images, replace NaN with a large value: 999
            mskEmpty[np.isnan(mskEmpty)] = 999
            # For detections, replace NaN with 0
            detEmpty[np.isnan(detEmpty)] = 0
        # Create a WCS for the combined image
        outWcs = apWcs.WCS(naxis=2)
        outWcs.wcs.crpix = [newCenX + 1, newCenY + 1]
        outWcs.wcs.crval = [ra, dec]
        outWcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        outWcs.wcs.cdelt = np.array([cdMatrix[0][0],
                                     cdMatrix[1][1]])
        # Output to header
        outHead = outWcs.to_header()
        outHead.set("PIXEL", pixScale, "Pixel Scale [arcsec/pix]")
        outHead.set("PHOZP", phoZp, "Photometric Zeropoint")
        outHead.set("EXPTIME", 1.0, "Set exposure time to 1 sec")
        outHead.set("GAIN", 3.0, "Average GAIN for HSC CCDs")
        for m in range(nReturn):
            outHead.set("TRACT" + str(m), trList[m])
            outHead.set("PATCH" + str(m), paList[m])

        # Define the output file name
        if verbose:
            print SEP
            print "### Generate Outputs"
        # Save the image array

        saveImageArr(imgEmpty, outHead, outPre + '_img.fits')
        if not imgOnly:
            # Save the mask array
            saveImageArr(mskEmpty, outHead, outPre + '_bad.fits')
            """ 15/12/12 Stop saving the variance plane"""
            # Save the variance array
            # saveImageArr(varEmpty, outHead, outPre + '_var.fits')
            # Save the sigma array
            saveImageArr(sigEmpty, outHead, outPre + '_sig.fits')
            # Save the detection mask array
            saveImageArr(detEmpty, outHead, outPre + '_det.fits')

        if nReturn > 0:
            cutFound = True
            # Save a preview image
            if visual:
                pngOut = outPre + '_pre.png'
                if not imgOnly:
                    previewCoaddImage(imgEmpty, mskEmpty, varEmpty, detEmpty,
                                      oriX=newX, oriY=newY, boxW=boxX,
                                      boxH=boxY, outPNG=pngOut)
        else:
            cutFound = False
            print WAR
            print "### No data was collected for " + \
                  "this RA,DEC in %s band!" % filt
    else:
        print WAR
        print "### No data was collected for this RA,DEC in %s band!" % filt
        cutFound = False
        cutFull = False

    return cutFound, cutFull, nReturn


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("ra",   type=float, help="RA  to search")
    parser.add_argument("dec",  type=float, help="Dec to search")
    parser.add_argument("size", type=float, help="Half size of the cutout box")
    parser.add_argument('-f', '--filter', dest='filt', help="Filter",
                        default='HSC-I')
    parser.add_argument('-p', '--prefix', dest='outfile',
                        help='Prefix of the output file',
                        default='hsc_coadd_cutout')
    parser.add_argument('-i', '--imgOnly', action="store_true",
                        dest='imgOnly', default=False)
    args = parser.parse_args()
    coaddImageCutFull(args.root, args.ra, args.dec, args.size,
                      filt=args.filt, prefix=args.outfile,
                      imgOnly=args.imgOnly)
