#!/tigress/HSC/LSST/stack_20160915/Linux64/miniconda2/3.19.0.lsst4/bin/python
from __future__ import print_function

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'

from astropy import units as u

try:
    import lsst.log
    Log = lsst.log.Log()
    Log.setLevel(lsst.log.ERROR)
except ImportError:
    pass

import lsst.daf.persistence
import lsst.afw.display.rgb as afwRgb
import lsstutils

ROOT = '/tigress/HSC/HSC/rerun/production-20160523'

def _get_skymap(root=ROOT):
    butler = lsst.daf.persistence.Butler(root)
    skymap = butler.get('deepCoadd_skyMap', immediate=True)
    print()
    return butler, skymap


def single_rgb_image(ra, dec, radius, prefix, Q=8., dataRange=0.6, scale=20, 
                     file_format='png', img_size=None, butler=None, 
                     skymap=None, root=ROOT):

    if butler is None:
        butler, skymap = _get_skymap(root)

    img = lsstutils.make_rgb_image(
        ra, dec, radius, Q=Q, dataRange=dataRange, 
        butler=butler, skymap=skymap, img_size=img_size)

    if img is not None:

        fig, ax = plt.subplots(
            subplot_kw={'yticks':[], 'xticks':[]})
        ax.imshow(img, origin='lower')

        if scale:
            shape = img.shape
            xmin = 15.0
            xmax = xmin + scale/0.168
            y=0.93*shape[0]
            ax.axhline(y=y, xmin=xmin/shape[1], xmax=xmax/shape[1], 
                       color='w', lw=3.0, zorder=1000)
            label = str(int(scale))
            ax.text((xmin+xmax)/2 - 0.042*shape[1], y - 0.072*shape[0], 
                    r'$'+label+'^{\prime\prime}$', color='w', fontsize=20)

        fig.savefig(prefix+'.'+file_format, bbox_inches='tight')
        plt.close('all')


def batch_rgb_images(cat_fn, radius, prefix, Q=8, dataRange=0.6, scale=20,
                     file_format='png', img_size=None, root=ROOT):
    from astropy.table import Table
    cat = Table.read(cat_fn)

    butler, skymap = _get_skymap(root)

    print('generating {} rgb images for...'.format(len(cat)))
    for num, obj in enumerate(cat):
        print('source:', num)
        new_prefix = prefix+'-'+str(num)
        single_rgb_image(
            obj['ra'], obj['dec'], radius, new_prefix, Q, dataRange, scale, 
            file_format, img_size, butler=butler, skymap=skymap)


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
        '-o', '--out_prefix', type=str, default=None,
        help='output file (single mode) or directory (batch mode) prefix.')
    parser.add_argument(
        '-Q', type=float, default=8, help='RGB function parameter')
    parser.add_argument(
        '--dataRange', type=float, default=0.6, help='RGB function parameter')
    parser.add_argument(
        '--format', type=str, default='png', help='file format')
    parser.add_argument(
        '--root', type=str, default=ROOT,
        help='Root data directory.')

    args = parser.parse_args()

    if args.out_prefix is None:
        print('***** must give output prefix *****')
        parser.print_help()
    elif args.single:
        ra, dec, = args.single
        single_rgb_image(
            ra, dec, args.radius, args.out_prefix, args.Q, 
            args.dataRange, file_format=args.format, root=args.root)
    elif args.batch:
        cat_fn = args.batch
        batch_rgb_images(
            cat_fn, args.radius, args.out_prefix, args.Q, 
            args.dataRange, file_format=args.format, root=args.root)
    else:
        print('***** must select single or batch mode *****')
        parser.print_help()
