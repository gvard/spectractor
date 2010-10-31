
'''This module is main part of Spectractor. It contains functions for
read/write binary data (ex. Dech data formats, FITS) and for manipulating
data arrays.

@author: Dmitry Nasonov
'''

import os, sys
from struct import unpack, pack, calcsize, error as StructError
import numpy as np

try:
    import pyfits
except ImportError, err:
    print >> sys.stderr, str(err), "PyFITS is required for read and write FITS",
    print >> sys.stderr, "files only in functions read_fits and write_fits."
    # Get it here: http://www.stsci.edu/resources/software_hardware/pyfits

try:
    from scipy.interpolate import splrep, splev, InterpolatedUnivariateSpline
except ImportError, err:
    print >> sys.stderr, str(err), "SciPy is required for spline fitting only",
    print >> sys.stderr, "in two functions: flat_order and interp_chain."


def read_bin(bin_pth, headlen=10, fmtdat='h', npuse=False, verbose=False):
    """Read binary file, for example *.100 files.
    100 format is a simple representation of 2D array with short integer values.
    The structure of given binary files: First 'headlen' bytes: Object name.
    Next 2 + 2 bytes: number of orders and number of pixels in one order.
    The rest is an array of integer values with given shape.
    @param bin_pth: path for input binary file.
    @param headlen: Length of text header, used for saving object name.
    @param fmtdat: struct format character for L{pack}/L{unpack}, usually
            one of the "bBhHi". Default is signed short.
    @param npuse: use NumPy and return data in np.array, lists otherwise.
    @return: object name, data.
    """
    bindat = open(bin_pth, 'rb')
    objname = bindat.read(headlen)
    ord_num, ord_len = unpack('HH', bindat.read(4))
    fmt_str = fmtdat * ord_len
    bsize = calcsize(fmtdat)
    if verbose:
        print "binary file header:", objname, ord_num, ord_len
    if npuse:
        data = np.zeros((ord_num, ord_len), int)
    else:
        data = [[] for _ in xrange(ord_num)]
    for i in xrange(ord_num):
        data[i] = unpack(fmt_str, bindat.read(bsize*ord_len))
    return (objname, ord_num, ord_len), data


def get_bin_head(bin_pth, headlen=10, verbose=False):
    """Get shape of the data and object name from binary (ex. *.100) file.
    @param bin_pth: path for input binary file.
    @param headlen: Length of text header, used for saving object name.
    @return: number of orders, length of one order, header 'object name' field.
    """
    bindat = open(bin_pth, 'rb')
    objname = bindat.read(headlen)
    ord_num, ord_len = unpack('HH', bindat.read(4))
    #bindat.close()
    if verbose:
        print "Header for binary file", bin_pth+":", objname, ord_num, ord_len
    return ord_num, ord_len, objname


def write_bin(data, bin_pth, objname="", headlen=10, fmtdat='h', verbose=False):
    """Write binary file, for example *.100 files.
    100 format is simple representation of 2D array with short integer values.
    @param data: 2D array of int values.
    @param bin_pth: path for output file.
    @param objname: Object name (or some other text). objname strip (or fit)
            to 'headlen' chars.
    @param headlen: length of objname to write. We strip or fit objname to given
            length. Default value is 10 (for "100" format).
    @param fmtdat: struct format character for L{pack}/L{unpack}, usually
            one of the "bBhHi". Default is signed short.
    """
    # Get data shape:
    ord_num = len(data)
    try:
        ord_len = len(data[0])
    except TypeError:
        ord_num, ord_len = 1, len(data)
        data = [data, ]
        if verbose:
            print 'write_bin: Input data is 1D! Set ord_num = 1.'
    # Start writing to file
    binf = open(bin_pth, 'wb')
    binf.write(objname[:headlen].ljust(headlen) + pack('HH', ord_num, ord_len))
    fmt_str = fmtdat * ord_len
    # Write data:
    for i in xrange(ord_num):
        try:
            binf.write(pack(fmt_str, *data[i]))
        except StructError, err:
            print >> sys.stderr, "StructError in write_bin:", err
            _coef = 2
            if fmtdat.isupper():
                _coef = 1
            maxval = 2**(8*calcsize(fmtdat))/_coef - 1
            for j in xrange(ord_len):
                try:
                    binf.write(pack(fmtdat, data[i][j]))
                except StructError:
                    if verbose:
                        print >> sys.stderr, "Overflow in", i, "order:",
                        print >> sys.stderr, j, "value", data[i][j]
                    binf.write(pack(fmtdat, maxval))
    binf.close()


def read_fits(file_pth, byteswap=False, delempty=True, verbose=False):
    """Read FITS, but with byteswap=False by default for *.200 files.
    @param file_pth: Path to input FITS file.
    @param byteswap: False for "200" format, True for normal FITS.
    @param delempty: delete empty cards from header.
    @param verbose: print header length, OBJECT card, data.shape
            and data.dtype.
    @return: header as PyFITS header and data as numpy.array.
    """
    #ftsdat = pyfits.open(file_pth)
    #head, data = ftsdat[0].header, ftsdat[0].data.byteswap(byteswap)
    head = pyfits.getheader(file_pth)
    if delempty:
        # Clean header: remove empty lines.
        del head[""]
    data = pyfits.getdata(file_pth).byteswap(byteswap)
    if byteswap:
        data = data.newbyteorder()
    if verbose:
        print "Header length:", len(head.ascard), "Object:", head['OBJECT'],
        print "Shape:", data.shape, "Data type:", data.dtype
    return head, data


def write_fits(data, fits_pth, objname="", byteswap=True, helicorr=0, hist=""):
    """Write FITS, but with byteswapped by default data for *.200 files.
    Please keep in mind, that FITS standard only allows up to 8 characters for
    the keyword name, but C{len('DECH-HELIOCORRECTION') == 20}. All questions
    to U{author of "200" format<http://gazinur.com>}. To avoid this, we use
    trick described in the PyFITS handbook (2009, p.35).
    Also, "200" format supports only integer values, not float.
    @param data: 2D array. may be list of lists or numpy.array with NAXIS=2.
    @param byteswap: True for "200" format, False for writing normal FITS.
    @param objname: object name. If not empty, write it to card OBJECT.
    @param hist: history. If not empty, write it to card HISTORY.
    @param helicorr: Write heliocentric correction in header, if byteswap.
    """
    data = np.array(data).byteswap(byteswap)
    if not byteswap:
        data = data.newbyteorder()
    # Create fits object
    fits = pyfits.PrimaryHDU(data)
    if objname:
        fits.header.update('OBJECT', objname)
    if byteswap and helicorr and type(helicorr) in (int, float):
        # We need it only for .200 files, which always have byteswap == True.
        hcard = pyfits.Card().fromstring("DECH-HELIOCORRECTION= "+str(helicorr))
        fits.header.ascardlist().append(hcard)
        # This code don't work because of incorrect keyword name:
        #fits.header.update('DECH-HELIOCORRECTION', helicorr)
    if hist:
        fits.header.add_history("This file created by Spectractor.")
        fits.header.add_history(hist)
    fits.writeto(fits_pth, output_verify='ignore')


_CRDSPL = lambda crd: crd.split('=')[1].split('/')[0].strip().strip("'").strip()

def get_fits_card(fts_pth, cardnam="OBJ"):
    """Read FITS header, while FITS card not equal END or cardnam.
    Pure Python, no PyFITS ;)
    @param cardnam: first chars of card name. Case insensitive.
    @type cardnam: string
    @return: value of last matching card or None if card not found.
    """
    fts, card = open(fts_pth, 'rb'), ''
    while card[:len(cardnam)] not in ('END'.ljust(len(cardnam)),
                                      cardnam.upper()):
        card = fts.read(80)
    if card[:len(cardnam)] == cardnam:
        return _CRDSPL(card)
    else:
        return


def read_fds(fds_pth, ord_num, ord_len):
    """Read .fds file with array of wavelengths.
    For reading fds we must know about data shape of reference spectrum:
    number of orders and length of one order. Output data is an array of
    wavelengths (floating point numbers) for each point of reference spectrum.
    @param fds_pth: path to fds file.
    @param ord_num: number of orders.
    @type ord_num: unsigned int
    @param ord_len: length of order.
    @type ord_len: unsigned int
    @return: numpy array of float wavelengths.
    """
    fds = open(fds_pth, 'rb')
    fmt_str = 'i' * ord_len * ord_num
    data = unpack(fmt_str, fds.read(4*ord_len*ord_num))
    fds.close()
    data = np.reshape(np.array(data) / 10000., (ord_num, ord_len))
    return data


def get_wv_range(fds_pth, shift=None):
    """Read only first and last wavelengths from .fds file, so get wl range.
    @param fds_pth: path to fds file.
    @param shift: seek to given position for reading last value.
        Default is os.path.getsize(fds_pth) - 4
    @return: two floats: first and last wavelength.
    """
    fds = open(fds_pth, 'rb')
    if not shift:
        shift = os.path.getsize(fds_pth) - 4
    beg = unpack('i', fds.read(4))[0] / 10000.
    fds.seek(shift)
    end = unpack('i', fds.read(4))[0] / 10000.
    fds.close()
    return beg, end


def write_fds(data, fds_pth, verbose=False):
    """Write .fds binary file with array of wavelengths.
    @param data: 2D array data (np.array or list).
    @param fds_pth: path for saving fds file.
    """
    # Get data shape:
    ord_num = len(data)
    try:
        ord_len = len(data[0])
    except TypeError:
        ord_len = ord_num
        ord_num = 1
        data = [data, ]
    if verbose:
        print "Write fds. data.shape:", ord_num, ord_len
    arr = np.around(np.array(data)*10000).astype(int).reshape(-1)
    fmtstr = 'i' * ord_len * ord_num
    fds = open(fds_pth, 'wb')
    fds.write(pack(fmtstr, *arr))
    fds.close()


def read_ccmtxt(ccm_pth):
    """Read *.ccm.txt files. This files creates by Dech 7.4.x.
    ccm.txt files is simple way to represent a continuum polynome for each
    order in echelle spectrum.
    @param ccm_pth: path to .ccm.txt file.
    @return: dictionary with points for spline fitting. keys are number
             of orders, values are lists of tuples with points.
    """
    dots = {}
    for dot in open(ccm_pth, 'r').read().splitlines():
        dot = dot.split(",")
        try:
            dots[int(dot[0])].append((float(dot[1]), float(dot[2])))
        except KeyError:
            dots[int(dot[0])] = [(float(dot[1]), float(dot[2]))]
        except (IndexError, ValueError):
            continue
    return dots


def write_ccmtxt(ccm_dct, ccm_pth):
    """Write *.ccm.txt files. This files creates by Dech 7.4.x.
    ccm.txt files is simple way to represent a continuum polynome for each
    order in echelle spectrum.
    @bug: difference for resulting .ccm.txt in rounding of floats:
    8.0 -> '8.0' (Dech: 8.0 -> '8').
    @param ccm_dct: dictionary with order numbers as keys and lists of
            polynome points as values.
    @param ccm_pth: path to .ccm.txt file.
    """
    ccm = open(ccm_pth, 'wb')
    for ordnum, pnt in sorted(ccm_dct.items()):
        for pnx, pny in sorted(pnt):
            ccm.write(",".join(map(str, (ordnum, pnx, pny)))+'\r\n')
    ccm.close()


def take_orders(data, fds, lam, polypt, edglim=5, ish=0, wcut=.5, olim=2):
    """Take orders from data which contains given wavelength.
    Here we assume that orders are sorted by wavelength.
    @param data: data as 2D array.
    @param fds: dispersion curve as 2D array of wavelengths.
    @param lam: wavelength, probably containing in fds' wavelength range.
    @param polypt: dictionary with points for flatten polynome plotting.
    @param edglim: index for breaking points at the edges.
    @param ish: intensity shift for flatted order.
    @param wcut: wavelength shift in angstroms: break order with
            'incomplete' line with wcut<0 or take it otherwise).
    @param olim: limit for number of extracted orders.
    @return: dictionary with tuples of fds[i], data[i] as values
            and order numbers as keys.
    """
    ords_dct = {}
    for i in xrange(len(data)):
        if fds[i][edglim]-wcut < lam < fds[i][-edglim]+wcut:
            if i+1 in polypt:
                ords_dct[i] = (flat_order(data[i], polypt[i+1],
                                          intlev=ish), fds[i])
            else:
                print "polynome points dict has no", i+1, "order, pass"
        if len(ords_dct) == olim:
            break
    return ords_dct


def flat_order(order, flatdots, s=None, k=3, intlev=0):
    """Make interpolation of given order by spline based on given points.
    This is wrapper function for SciPy spline interpolation. Description
    of s and k parameters are from docstrings of scipy.interpolate.splrep.
    @param order: 1D array of int data.
    @param flatdots: points for plot flat spline.
    @type flatdots: tuple or list of [dotsx], [dotsy]
    @param k: The order of the spline fit. It is recommended to use cubic
            splines; 1 <= k <= 5.
    @param s: A smoothing condition. The amount of smoothness is determined by
    satisfying the conditions: sum((w * (y - g))**2,axis=0) <= s where g(x) is
    the smoothed interpolation of (x, y). Larger s means more smoothing while
    smaller values of s indicate less smoothing. Recommended values of s depend
    on the weights, w. default: s=m-sqrt(2*m) if weights are supplied.
    @param intlev: Intensity shift for flatted order. By default
            0 <= output_data <= 1.
    @return: normalized order.
    """
    pixnums = np.arange(1, len(order)+1, 1)
    xflat, yflat = zip(*flatdots)
    try:
        tckp = splrep(xflat, yflat, s=s, k=k)
    except ValueError, err:
        print >> sys.stderr, "ValueError in flat_order:", err,
        print >> sys.stderr, "Maybe you should sort input arrays."
        raise SystemExit
    return order / splev(pixnums, tckp) + intlev


def interp_chain(lams, dats, mf=5, k=1):
    """Interpolate chain with multiply factor mf.
    @param f: Multiply factor. length of output arrays is calculated as
            length of input arrays multiplied by mf.
    @param k: spline degree.
    @return interpolated input arrays.
    """
    ius = InterpolatedUnivariateSpline(lams, dats, k=k)
    outlams = np.linspace(lams[0], lams[-1], len(dats)*mf)
    return outlams, ius(outlams)


def get_line(lams, lam, Va=0, width=.9, off=0, cutedg=4, vscl=False,
             verbose=False):
    """Get the spectrum fds and wavelength, find limits and return them.
    First we calculate dlam = (Va*lam)/C, where C is a speed of light in km/s.
    Then we get range: lam1, lam2 = (lam+dlam-width-off, lam+dlam+width-off).
    Finally, find wavelengths nearest to lam1, lam2 in lams and return
    slice of lams.
    @param lams: array of wavelengths. Must be np.array!
    @param lam: reference wavelength.
    @param Va: spectrum shift in km/s, for example, heliocentric correction.
    @param width: get piece of lams with length equal to 2*width angstroms.
    @param off: spectrum shift in angstroms, for example, systemic velocity
            shift.
    @param cutedg: drop given number of pixels on edges of the order.
    @param vscl: scale mode. If True, transform wavelength scale to radial
                 velocities scale.
    @return: array of wavelengths or velocities, start and end indices
            of input array' slice.
    """
    _C = 299792.5
    vadlam = (Va * lam) / _C
    if vadlam:
        lams += vadlam
    lam1, lam2 = lam - width - off, lam + width - off
    # Get the slice of lams!
    try:
        beg = max(np.where(lams>lam1)[0][0], cutedg)
        end = min(np.where(lams<lam2)[0][-1], len(lams)-cutedg)
    except IndexError, err:
        print >> sys.stderr, "IndexError in get_line:", err
        return (), None, None
    lams = lams[beg:end]
    if verbose:
        print "Get line. delta lam:", round(vadlam, 3),
        print "l1, l2:", round(lam1, 3), round(lam2, 3),
        print "beg, end indices:", beg, end, "length:", len(lams)
    if vscl:
        # Simply v = _C * deltalam / lam
        lams = (lams - lam) * _C / lam
    return lams, beg, end

class PrepareTfits:
    """Split (UVES) spectrum in tfits format into orders for Dech software.
    First we reshape structured numpy array, split spectrum to non-zero chains,
    than split each chain into 'orders' with overlap and construct separate
    arrays for data and wavelengths as list of orders.
    @param data: structured numpy array with values containing wavelength, data
            value and data error with dtype '>f8'. UVES tfits has this format
    """
    def __init__(self, data, dtype='>f8', dlen=3, verbose=False):
        self.verbose = verbose
        # Reshape structured numpy array:
        data = data.view(dtype).reshape(-1, dlen)
        lams, cha = data[:, 0], data[:, 1]
        self.maxarr = np.max(cha)
        # Split to non-zero chains
        self.split_tfits(cha, lams)

    def split_tfits(self, cha, lams):
        """Split spectrum with its fds to non-zero chains"""
        self.arr, self.lam = [], []
        # Find all zeros:
        azers, = np.where(cha == 0)
        if azers[0]:
            self.arr.append(cha[:azers[0]])
            self.lam.append(lams[:azers[0]])
            if self.verbose:
                print "split_tfits chain:", azers[0], cha[azers[0]]
        for i, pix in enumerate(azers):
            if i and pix - azers[i-1] > 1:
                self.arr.append(cha[azers[i-1]:pix])
                self.lam.append(lams[azers[i-1]:pix])
                if self.verbose:
                    print "split_tfits chain:", azers[i-1], pix, cha[i-1]
        if azers[-1] != len(cha)-1:
            self.arr.append(cha[azers[-1]:])
            self.lam.append(lams[azers[-1]:])
            if self.verbose:
                print "split_tfits chain:", azers[-1], cha[azers[-1]]
        self.lens = [len(chain) for chain in self.arr]
        if self.verbose:
            print "Num of chains", len(self.arr), "Min length", min(self.lens),
            print "max val", self.maxarr, "lengths array:", self.lens

    def get_split(self, olen=6350, ovr=350, bigval=False, maskbad=True):
        """Split each chain in array to 'orders' with equal length.
        @param data: raw data from UVES tfits as structured array of
                wavelengths, counts and errors. Errors are ignored
        @param maskbad: treat values <0 as 'bad' and replace them with previous
                    value
        @param bigval: Multiply spectrum values by big constant factor
                    900000/max(data) instead of 32767/max(data)
        @param ordlen: 'order' length
        @param ovr: order's overlap length
        @return data and fds as lists of orders with equal length
        """
        if bigval:
            divr = int(900000./self.maxarr)
        else:
            divr = int(32767./self.maxarr)
        if self.verbose:
            print "Set division coefficient", divr
        self.ndat, self.nlam = [], []
        for j, chain in enumerate(self.arr):
            if divr:
                chain = np.around(chain*divr).astype(int)
            self.split_to_orders(chain, self.lam[j], olen, ovr, maskbad)
        return self.ndat, self.nlam

    def split_to_orders(self, chain, lam, ordlen=6350, ovr=350, maskbad=False):
        """Split spectrum chain into 'orders' for Dech software.
        @param chain: raw data as array of floats/ints.
        @param lam: wavelengths as array of floats/ints.
        @param ordlen: order length
        @param ovr: overlap length
        @param maskbad: treat values <0 as 'bad' and replace them with previous
                value
        """
        if maskbad and np.where(chain < 0)[0].all():
            for i in np.where(chain < 0)[0]:
                if self.verbose:
                    print "Replace point < 0 in chain", i, chain[i], chain[i-1]
                chain[i] = chain[i-1]
        for i in xrange(len(chain) / ordlen):
            self.ndat.append(chain[i*ordlen:(i+1)*ordlen+ovr])
            self.nlam.append(lam[i*ordlen:(i+1)*ordlen+ovr])
        ostk = (i + 1) * ordlen
        lastord = np.concatenate((chain[ostk-ovr:],
                                  np.zeros(ordlen-len(chain)+ostk, dtype=int)))
        lastlam = np.concatenate((lam[ostk-ovr:], np.linspace(lam[-1],
                                  lam[-1]+40, ordlen-len(chain)+ostk)))
        self.ndat.append(lastord)
        self.nlam.append(lastlam)


def split_array(data, ordlen=4000, ovr=100):
    """Split single big array on equal chains.
    @param ordlen: length of one chain ("order")
    @param ovr: length of chains overlapping
    @return list of chains.
    """
    ndata = []
    for i in xrange(len(data)/ordlen-1):
        ndata.append(data[i*ordlen:(i+1)*ordlen+ovr])
    ndata.append(data[(i+1)*ordlen-ovr:])
    return ndata


def int_split_array(data, ordlen=4000, ovr=50, divr=30., maskbad=True):
    """Split single big array on equal chains.
    @param data: data as list or np.array
    @param ordlen: length of one chain ("order")
    @param ovr: length of chains overlapping
    @param divr: division coefficient. If zero, skip dividing
    @param maskbad: rewrite bad pixels
    @return data as list
    """
    if divr:
        data = np.around(data/divr).astype(int)
    if maskbad:
        for i, pnt in enumerate(data):
            if pnt < 0:
                data[i] = data[i-1]
    return split_array(data, ordlen, ovr)


def gen_split_array(ordlen=4000, ovr=100, beg=4000, end=6800, step=.05):
    """Generate array with np.arange and split them on equal chains.
    Default values are for making ELODIE spectra FDS.
    @param ordlen: length of one chain ("order")
    @param ovr: length of chains overlapping
    @param beg: start of interval
    @param end: end of interval
    @param step: spacing between values
    @return list of chains.
    """
    return split_array(np.arange(beg, end, step), ordlen, ovr)
