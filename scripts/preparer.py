#!/usr/bin/env python

"""Spectractor script for logging and preparing raw FITS data.
"""

import os
from optparse import OptionParser

from spectractor.serve import Walker
from spectractor.midworker import Preparer

def strclbck(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

desc = "Spectractor script for logging and preparing raw FITS data."
usage = "Usage: %prog [options] [input_dir] [input_file]"
parser = OptionParser(description=desc, usage=usage, version="0.1.0")
parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
    default=False, help="Be more verbose")
parser.add_option("-e", "--ext", type="str", dest="ext", action="callback",
    callback=strclbck, default=("fits", "mt"),
    help="Set files extension, e.g. mt, fits")
parser.add_option("-c", "--lims", type="int", dest="lims", default=None,
    nargs=2, help="Set limits on Y-axis cuts on NES images as two ints")
parser.add_option("-F", "--flats", type="int", dest="flats", default=None,
    nargs=2, help="Set flat field image numbers as two ints")
parser.add_option("-l", "--logonly", action="store_true", dest="logonly",
    default=False, help="Make a log and save it to file, without processing")
parser.add_option("-u", "--usedisp", action="store_false", dest="disp",
    default=True, help="Use display.")
parser.add_option("-f", "--filemove", action="store_false", dest="filmove",
    default=True, help="Move initial files.")
parser.add_option("-b", "--badrow", action="store_true", dest="badrow",
    default=False, help="Use badrow for removing bad rows in NES images.")
parser.add_option("-p", "--outpath", type="str", dest="outpath",
    default="arch", help="Set output path. Default is 'arch'")
parser.add_option("-B", "--biasname", type="str", dest="biasname",
    default="bias", help="Set BIAS file name without extension")
(opts, args) = parser.parse_args()
if not args:
    args.append('./')

walk = Walker(args, exts=opts.ext, verbose=opts.verbose)
walk.walk()
files = walk.select_images()

prg = Preparer(dest_dir=opts.outpath, usedisp=opts.disp, verbose=opts.verbose)
if opts.lims:
    prg.imgcut_dct[2052] =  opts.lims + list(prg.imgcut_dct[2052][2:])
    print "Set cuts", prg.imgcut_dct[2052]
prg.prepare(files, log="logonly")
prg.savelog(wdir=opts.outpath)
if not opts.logonly:
    if opts.flats:
        print "Get flats from", opts.flats[0], "to", opts.flats[1]
        prg.Lg.flats = prg.Lg.select_nums(opts.flats, prfx="f")
    if not os.path.isfile(opts.biasname+".fits"):
        prg.crebias(filmove=opts.filmove, biasname=opts.biasname)
    prg.prepare(files, log=False, badrow=opts.badrow, filmove=opts.filmove,
                substbias=opts.biasname+".fits")
    prg.creflat()
prg.mw.bye()
