#!/usr/bin/env python

"""Spectractor script for ploting lines in spectra vectors.

Some optical line wavelengths:
"Ca": [3933.663, 3968.469],
"Na": [5889.951, 5895.923],
"H":  [4101.73, 4340.46, 4861.32, 6562.8],
"He": [4026.21, 4120.84, 4387.93, 4471.5, 4713.17, 4921.93,
       5015.68, 5875.66, 6678.15, 7065.25],
"[O1]": [5577.339, 6300.304, 6363.776],
"DIB": [5780.55, 5796.99, 6195.99, 6203.06, 6613.63],
"Si3": [4552.622, 4567.84, 4574.757, 5739.734],

Example usage:
liner.py ./ -l 4861.32,6562.8 -e 100,200
liner.py ./ -l 5889.951,5895.923 -i l,n -n 60
"""

import os
#import cPickle
from optparse import OptionParser
from itertools import cycle

#if X display is not avaliable, try this:
#import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from spectractor.spectractor import Spectractor
from spectractor.serve import Walker, read_journal, mod_journal
from spectractor.plotter import plt_opt, Clicker


def floclbck(option, opt, value, parser):
  setattr(parser.values, option.dest, map(float, value.split(',')))

def strclbck(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

desc = "Spectractor script for plotting data arrays"
usage = "Usage: %prog [options] [input_dir] [input_file]"
parser = OptionParser(description=desc, usage=usage, version="0.1.0")
parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
    default=False, help="Be more verbose")
parser.add_option("-e", "--ext", type="str", dest="ext", action="callback",
    callback=strclbck, default=["200"],
    help="Set files extension. Can be 100, 200 or fits")
parser.add_option("-f", "--flag", type="str", dest="flag", action="callback",
    callback=strclbck, default=None,
    help="Set possible flags. If not set, process all")
parser.add_option("-i", "--inst", type="str", dest="inst", action="callback",
    callback=strclbck, default=None,
    help="Set instrumental keys. If not set, process all")
parser.add_option('-l', '--lam', type="str", dest="lam", action="callback",
    callback=floclbck, default=[5889.951],
    help="Set center wavelengths for output chains, e.g. in angstroms")
parser.add_option('-F', '--figsz', type="str", action="callback", dest="figsz",
    default=(7, 6), callback=floclbck, help="Set size of the figure")
parser.add_option('-x', '--xlims', type="str", action="callback",
    dest="xlims", default=None, callback=floclbck,
    help="Set limits on X axis of the figure")
parser.add_option("-w", "--cut", type="float", dest="window", default=4,
    help="Set width of spectral window e.g. in angstroms")
parser.add_option("-o", "--shift", type="float", dest="shift", default=0,
    help="Set shift of the spectral window")
parser.add_option("-s", "--vscale", dest="vscl", action="store_false",
    default=True, help="Set wavelength scale, velocity scale is default")
parser.add_option("-n", "--sn", type="float", dest="sn", default=1,
    help="Set limit for S/N ratio")
parser.add_option("-r", "--rvcorr", dest="rvcorr", action="store_false",
    default=True, help="Use rvcorrect. Useful for telluric lines")
parser.add_option("-b", "--bw", dest="greyscale", action="store_true",
    default=False, help="Set color/greyscale mode for plotting")
parser.add_option("-g", "--legend", dest="legend", action="store_true",
    default=False, help="Plot a legend")
(opts, args) = parser.parse_args()
if not args:
    args.append('./')
if not opts.rvcorr:
    opts.shift = 0

walk = Walker(args, exts=opts.ext, verbose=opts.verbose)
walk.walk()

try:
    from data.obj import rvcorrect_dct
    from data.cuts import cuts_dct
except ImportError:
    rvcorrect_dct, cuts_dct = {}, {}
journ_pth = os.path.join(os.path.dirname(args[0]), 'data', 'journal.pkl')
if os.path.isfile(journ_pth):
    journal = read_journal(journ_pth)
else:
    print "No journal", journ_pth
    journal = {}

limlam = (min(opts.lam), max(opts.lam))
mod_journal(journal, opts.sn, rvc=opts.rvcorr, rvcorr_dct=rvcorrect_dct,
            instr=opts.inst, lams=limlam, verbose=opts.verbose)

for ext in opts.ext:
     walk.select_spectra(journal, ext=ext, iflag=opts.flag)
spectra_parms = walk.spath_set

if not len(spectra_parms):
    print "Sorry, no spectra!"
    raise SystemExit

if opts.verbose:
    print "Length of journal is", len(journal)

sp = Spectractor(spectra_parms, vscl=opts.vscl, cuts_dct=cuts_dct, dtype=list,
                edglim=4, verbose=opts.verbose)

for lam in opts.lam:
    if opts.xlims:
        sp.get_chains(lam, llim=opts.xlims[0]-10, hlim=opts.xlims[1]+10)
    else:
        sp.collect_chains(lam, width=opts.window, shift=opts.shift)
        opts.xlims = sp.get_lims(lam)
findata = sp.fdata

if not len(findata):
    print "No chains, sorry. Exit."
    raise SystemExit

print "From", len(sp.spectra), "spectra we select", len(findata), "chains"

#datafile = open(os.path.join('out', 'liner_findata.pkl'), 'wb')
#cPickle.dump(findata, datafile)
#datafile.close()

fig = plt.figure(figsize=opts.figsz, dpi=110)
figsub = fig.add_subplot(111)

#Some plotting properties are here
pltopt = dict(xlbl="Heliocentric radial velocity, km/s",
        ylbl="Relative intensity", lblsz=10, ticksz=10)
adjprops = dict(left=.12, right=.99, top=.98, bottom=.07)
if opts.xlims:
    pltopt.update({"xlim": opts.xlims})
if 5889.951 in opts.lam:
    pltopt.update(dict(ylim=(0, 1.2)))
clr = cycle(('r', 'g', 'b', 'k', 'c', 'm', 'olive', 'grey'))
if opts.greyscale:
    clr = cycle(('k', 'grey'))

fig.subplots_adjust(**adjprops)

#for i, (label, (lams, dats)) in enumerate(findata.items()):
for label, lams, dats in findata:
    figsub.plot(lams, dats, c=clr.next(), label=str(label[:2]), lw=0.8)

Clck = Clicker(figsub, xlims=opts.xlims, crdlim=30)
figsub.set_autoscale_on(False)
plt.connect('button_press_event', Clck)
plt.connect('key_press_event', Clck)

if pltopt:
    if opts.legend:
        pltopt['legend'] = .007 # spacing between labels in legend
        pltopt['loc'] = 'best' # location of the legend
    plt_opt(figsub, **pltopt)

plt.show()
print "clicked coordinates:", Clck.coords
