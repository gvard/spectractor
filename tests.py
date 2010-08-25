#!/usr/bin/env python

'''Unit tests for spectractor python module.
After finishing tests remove tmp files.
Created on 13.03.2010
@author: Dmitry Nasonov
'''

from random import randint
import unittest

import numpy as np

try:
    from spectractor import *
except ImportError:
    print "Cannot import spbin module, check existence!"
    raise SystemExit


class TestFunctions(unittest.TestCase):
    def setUp(self):
        self.ordnum, self.ordlen = 25, 2040
        self.spec = [[randint(0, 32767) for _ in xrange(self.ordlen)]
                     for _ in xrange(self.ordnum)]
        self.fspec = [[randint(0, 650000) for _ in xrange(self.ordlen)]
                      for _ in xrange(self.ordnum)]
        self.spec[2][self.ordlen-20] = 32769
        self.fds = [list(np.linspace(4400+(70*i)-6, 4470+(70*i)+6, self.ordlen))
                    for i in xrange(self.ordnum)]
        self.spec1d = [randint(0, 32767) for _ in xrange(self.ordlen)]
        self.fds1d = list(np.linspace(4390, 4476, self.ordlen))
        self.polyp_dct = {}
        for ornm in xrange(self.ordnum):
            self.polyp_dct[ornm+1] = [
                 (randint(0, self.ordlen*100)/100.,
                  np.mean(self.spec[ornm])+randint(-100, 100)/50.)
                 for _ in xrange(randint(6, 24))]

    def test_rw_bin(self, hl=10, fmt='h', verbose=False):
        bin_pth, fds_pth = 'tmp.100', 'tmp.fds'
        for spec, fds in ((self.spec, self.fds), (self.spec1d, self.fds1d)):
            write_bin(spec, bin_pth, 'Test 2D spectrum',
                      headlen=hl, fmtdat=fmt, verbose=verbose)
            h, d = read_bin(bin_pth, headlen=hl, fmtdat=fmt, verbose=verbose)
            hh = get_bin_head(bin_pth, headlen=hl, verbose=verbose)
            if verbose:
                print 'Right length?', len(d) == hh[0],
                print len(d[0]) == self.ordlen
                onum = min(3, hh[0]-1)
                if onum == 3:
                    print 'Are orders equal?', list(d[onum]) == spec[onum]
                else:
                    print 'Are orders equal?', list(d[onum]) == spec
            write_fds(fds, fds_pth, verbose=verbose)
            onum, olen = hh[0], hh[1]
            fds = read_fds(fds_pth, onum, olen)
            wvrang = get_wv_range(fds_pth, onum, olen)
            if verbose:
                print 'Wavelength range:', wvrang

    def test_rw_fits(self, hc=6.4, crdnm="DECH-HELIOCORRECTION", verbose=False):
        fts_pth = 'tmp.200'
        write_fits(self.fspec, fts_pth, 'Test 2D FITS spectrum',
                   helicorr=hc, hist="A lot of history")
        head, dat = read_fits(fts_pth, verbose=verbose)
        ftscrd = get_fits_card(fts_pth, cardnam=crdnm)
        if verbose:
            print "Are orders equal?", list(dat[3]) == list(self.fspec[3])
            print "FITS card:", ftscrd, head.get(crdnm)

    def test_spectr_func(self, lam=4320.67, verbose=False):
        """Test some functions for manipulating data arrays"""
        write_ccmtxt(self.polyp_dct, 'tmp.ccm.txt')
        ccm_dct = read_ccmtxt('tmp.ccm.txt')
        if verbose:
            print "Are ccm.txt equal?", ccm_dct == sorted(self.polyp_dct)
        take_orders(self.spec, self.fds, lam, ccm_dct, olim=3)
        get_line(np.array(self.fds[0]), lam, Va=-12.3, verbose=verbose)
        interp_chain(self.fds1d, self.spec1d)
        gen_split_array(3500, 200)
        int_split_array(np.arange(1000, 15000, .5))
        #import pyfits
        #dat = pyfits.getdata('tmp.tfits')
        #prt = PrepareTfits(dat, verbose=verbose)
        #prt.get_split()

if __name__ == '__main__':
    unittest.main()
