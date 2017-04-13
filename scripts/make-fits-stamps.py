#!/tigress/HSC/LSST/stack_20160915/Linux64/miniconda2/3.19.0.lsst4/bin/python
from __future__ import print_function

from astropy import units as u

import lsst.log
Log = lsst.log.Log()
Log.setLevel(lsst.log.ERROR)

import lsst.daf.persistence
import lsst.afw.display
import lsstutils

def _get_skymap():
    root = '/tigress/HSC/HSC/rerun/production-20160523'
    butler = lsst.daf.persistence.Butler(root)
    skymap = butler.get('deepCoadd_skyMap', immediate=True)
    print()
    return butler, skymap


def single_stamp(ra, dec, radius, band='i', display=0, prefix=None, 
                 save_psf=False, butler=None, skymap=None):
    
    if butler is None:
        butler, skymap = _get_skymap()

    stamp = lsstutils.make_stamp(
        ra, dec, radius*u.arcsec, band=band, 
        butler=butler, skymap=skymap, return_psf=save_psf)
    if save_psf:
        stamp, psf = stamp

    if display:
        disp = lsst.afw.display.Display(display)
        disp.setMaskTransparency(70)
        disp.mtv(stamp)

    if prefix is not None:
        out_fn = prefix+'.fits'
        stamp.writeFits(out_fn)
        if save_psf:
            out_fn = prefix+'-psf.fits'
            psf.writeFits(out_fn)


def batch_stamps(cat_fn, radius, band, prefix, save_psf=False):
    from astropy.table import Table
    butler, skymap = _get_skymap()
    cat = Table.read(cat_fn, format='csv')

    print('making stamps...')
    for num, obj in enumerate(cat):
        print('cutting {}-band image for source {}'.format(band, num+1))
        stamp = lsstutils.make_stamp(
            obj['ra'], obj['dec'], radius*u.arcsec,
            band=band, butler=butler, skymap=skymap, return_psf=save_psf)
        if save_psf:
            stamp, psf = stamp
        out_fn = prefix+'-{}-{}.fits'.format(band, num+1)
        stamp.writeFits(out_fn)
        if save_psf:
            out_fn = prefix+'-{}-{}-psf.fits'.format(band, num+1)
            psf.writeFits(out_fn)


if __name__=='__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        '-s', '--single', type=float, nargs=2, default=None,
        help='single mode: ra dec (in deg)')
    parser.add_argument(
        '-b', '--batch', type=str, default=None,
        help='batch mode: csv catalog (with ra & dec) file name')
    parser.add_argument('-f', '--filter', type=str, default='i')
    parser.add_argument(
        '-r', '--radius', type=float, default=35,
        help='angular radius of cutout in arcsec')
    parser.add_argument(
        '-d', '--display', type=int, default=0, 
        help='display with ds9 (single mode)')
    parser.add_argument(
        '-o', '--out_prefix', type=str, default=None, 
        help='output file (single mode) or directory (batch mode) prefix.')
    parser.add_argument('--save_psf', action='store_true')
    args = parser.parse_args()

    if args.single:
        ra, dec, = args.single
        single_stamp(ra, dec, args.radius, args.filter, 
                     args.display, args.out_prefix, args.save_psf)
    elif args.batch:
        cat_fn = args.batch
        batch_stamps(cat_fn, args.radius, args.filter, 
                     args.out_prefix, args.save_psf)
    else:
        print('***** must select single or batch mode *****')
        parser.print_help()
