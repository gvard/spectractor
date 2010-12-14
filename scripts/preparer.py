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
    callback=strclbck, default=["fits", "mt"],
    help="Set files extension, e.g. mt, fits")
parser.add_option("-l", "--log", action="store_false", dest="log", default=True,
    help="Set prepare mode instead of logonly.")
parser.add_option("-L", "--logonly", action="store_true", dest="logonly",
    default=False, help="Set prepare mode instead of logonly.")
parser.add_option("-u", "--usedisp", action="store_false", dest="disp",
    default=True, help="Use display.")
parser.add_option("-f", "--filemove", action="store_false", dest="filmove",
    default=True, help="Move initial files.")
parser.add_option("-b", "--badrow", action="store_true", dest="badrow",
    default=False, help="Use badrow for removing bad rows in NES images.")
parser.add_option("-p", "--outpath", type="str", dest="outpath",
    default="arch", help="Set output path. Default: arch")
parser.add_option("-B", "--biases", type="str", dest="biases", default=[""],
    action="callback", callback=strclbck, help="Set biases")
(opts, args) = parser.parse_args()
if not args:
    args.append('./')
if os.path.isdir(args[0]):
    wdir = args[0]
else:
    wdir = "./"

walk = Walker(args, exts=opts.ext, verbose=opts.verbose)
walk.walk()
files = walk.select_images()

prg = Preparer(dest_dir=opts.outpath, usedisp=opts.disp, verbose=opts.verbose)
biasname = "bias"
prg.prepare(files, log="logonly")
if not os.path.isfile(biasname+".fits"):
    prg.crebias(filmove=opts.filmove, biasname=biasname)
prg.savelog(wdir=wdir)
prg.prepare(files, log=False, badrow=opts.badrow, filmove=opts.filmove,
            substbias=biasname+".fits")
prg.creflat()
