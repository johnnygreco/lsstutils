#!/tigress/HSC/LSST/stack_20160915/Linux64/miniconda2/3.19.0.lsst4/bin/python
from __future__ import print_function

from astropy import units as u

import lsst.log
Log = lsst.log.Log()
Log.setLevel(lsst.log.ERROR)

import lsst.daf.persistence
import lsst.afw.display
import lsstutils

root = '/tigress/HSC/HSC/rerun/production-20160523'
butler = lsst.daf.persistence.Butler(root)
skymap = butler.get('deepCoadd_skyMap', immediate=True)
print()

def single_stamp(ra, dec, radius, band='i', display=0, out_fn=None):
    stamp = lsstutils.make_stamp(
        ra, dec, radius*u.arcsec, band=band, butler=butler, skymap=skymap)

    if display:
        disp = lsst.afw.display.Display(display)
        disp.setMaskTransparency(70)
        disp.mtv(stamp)

    if out_fn is not None:
        stamp.writeFits(out_fn)


def batch_stamps():
    pass


if __name__=='__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('ra', type=float)
    parser.add_argument('dec', type=float)
    parser.add_argument('radius', type=float, help='radius in arcsec')
    parser.add_argument('-b', '--band', type=str, default='i')
    parser.add_argument('--display', type=int, default=0)
    parser.add_argument('-o', '--out_fn', type=str, default=None)

    args = parser.parse_args()
    single_stamp(args.ra, args.dec, args.radius, args.band, args.display, args.out_fn)
