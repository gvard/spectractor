#!/usr/bin/env python

"""Spectractor script for fitting spectral orders with spline
"""

import os, shutil
from optparse import OptionParser

import pyfits as pf
import matplotlib.pyplot as plt

from spctr.plotter import Flatter
from spctr.spectractor import read_fits, read_ccmtxt, write_ccmtxt

def intclbck(option, opt, value, parser):
  setattr(parser.values, option.dest, map(int, value.split(',')))

desc = "Spectractor script for fitting spectral orders with spline"
usage = "Usage: %prog [options] [input_dir] [input_file]"
parser = OptionParser(description=desc, usage=usage, version="0.1.0")
parser.add_option("-o", "--orders", type="str", dest="ords", action="callback",
    callback=intclbck, default=(1, 2), help="Set orders to flat")
parser.add_option('-f', '--file', type="str", dest="outfile",
    default="flat_z", help="Set file name for fits with divided data")
parser.add_option('-c', '--ccm', type="str", dest="ccm",
    default=None, help="Set file name for ccm.txt file with spline points")
parser.add_option('-b', '--blaze', type="str", dest="blz",
    default=None, help="Set file name for fits with blaze curve")
parser.add_option("-a", "--add", type="float", dest="add", default=None,
    help="Set shift for orders with minimum data value below zero")
(opts, args) = parser.parse_args()
if not args:
    raise SystemExit

#dat = read_fits(args[0])[1]
dat = pf.getdata(args[0])
subdata = []
if opts.blz and os.path.isfile(opts.blz):
    blz = pf.getdata(opts.blz)
    subblz = []
else:
    subblz = None
if opts.ccm and os.path.isfile(opts.ccm):
    ccm = read_ccmtxt(opts.ccm)
    subccm = {}
else:
    subccm = None
for i, d in enumerate(dat):
    if i+1 in opts.ords:
        if opts.add is not None and min(d) < 0:
            add = abs(opts.add - min(d))
            d += add
            print i, "order shifted to", add, "counts"
        subdata.append(d)
        if opts.blz:
            subblz.append(blz[i])
        if opts.ccm and i+1 in ccm:
            subccm[len(subdata)] = ccm[i+1]

fig = plt.figure()
figsub = fig.add_subplot(111)
adjprops = dict(left=.05, right=.99, top=.99, bottom=.05)
fig.subplots_adjust(**adjprops)
plc = Flatter(figsub, subdata, blaze=subblz, ccm=subccm, grid=True, lw=1.5)
plt.connect('button_press_event', plc)
plt.connect('key_press_event', plc)
plt.show()

df = pf.PrimaryHDU(plc.out)
if os.path.isfile(opts.outfile+'.fits'):
    shutil.move(opts.outfile+'.fits', opts.outfile+'_old.fits')
df.writeto(opts.outfile+'.fits')

if plc.outccm:
    if os.path.isfile(opts.outfile+'.ccm.txt'):
        shutil.move(opts.outfile+'.ccm.txt', opts.outfile+'_old.ccm.txt')
    write_ccmtxt(plc.outccm, opts.outfile+'.ccm.txt')
