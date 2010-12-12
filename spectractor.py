
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
    print >> sys.stderr, str(err), "SciPy is required for spline fitting"

_C = 299792.458

def read_bin(bin_pth, headlen=10, fmtdat='h', npuse=True, verbose=False):
    """Read binary file, for example *.100 files.
    100 format is a simple representation of 2D array with short integer values.
    The structure of given binary files: First 'headlen' bytes: Object name.
    Next 2 + 2 bytes: number of orders and number of pixels in one order.
    The rest is an array of integer values with given shape.
    @param bin_pth: path for input binary file
    @param headlen: Length of text header, used for saving object name
    @param fmtdat: struct format character for L{pack}/L{unpack}, usually
        one of the "bBhHi". Default is signed short
    @param npuse: use NumPy and return data in np.array, lists otherwise
    @return: object name, data
    """
    bindat = open(bin_pth, 'rb')
    objname = bindat.read(headlen)
    ord_num, ord_len = unpack('HH', bindat.read(4))
    fmt_str = fmtdat * ord_len
    bsize = calcsize(fmtdat)
    if verbose:
        print "binary file header:", objname, ord_num, ord_len
    if npuse:
        data = np.zeros((ord_num, ord_len), np.int32)
    else:
        data = [[] for _ in xrange(ord_num)]
    for i in xrange(ord_num):
        data[i] = unpack(fmt_str, bindat.read(bsize*ord_len))
    return (objname, ord_num, ord_len), data


def get_bin_head(bin_pth, headlen=10, verbose=False):
    """Get shape of the data and object name from binary file.
    @param bin_pth: path for input binary file
    @param headlen: Length of text header, used for saving object name
    @return: number of orders, length of one order, header 'object name' field
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
    @param data: 2D array of int values
    @param bin_pth: path for output file
    @param objname: Object name (or some other text). objname strip (or fit)
        to 'headlen' chars
    @param headlen: length of objname to write. We strip or fit objname to given
        length. Default value is 10 (for "100" format)
    @param fmtdat: struct format character for L{pack}/L{unpack}, usually
        one of the "bBhHi". Default is signed short
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
    @param file_pth: Path to input FITS file
    @param byteswap: False for "200" format, True for normal FITS
    @param delempty: delete empty cards from header
    @param verbose: print header length, OBJECT card, data.shape
        and data.dtype
    @return: header as PyFITS header and data as numpy.array
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

    Also, "200" format supports only integer values, not float. The best way
    for converting array of floats is use following commands:

        >>> import numpy as np
        >>> data = np.around(data).astype(np.int32)

    @param data: 2D array. may be list of lists or numpy.array with NAXIS=2
    @param byteswap: True for "200" format, False for writing normal FITS
    @param objname: object name. If not empty, write it to card OBJECT
    @param hist: history. If not empty, write it to card HISTORY
    @param helicorr: Write heliocentric correction in header, if byteswap
    @raise IOError: file fits_pth already exist.
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
    @param cardnam: first chars of card name. Case insensitive
    @type cardnam: string
    @return: value of last matching card or None if card not found
    """
    fts, card = open(fts_pth, 'rb'), ''
    while card[:len(cardnam)] not in ('END'.ljust(len(cardnam)),
                                      cardnam.upper()):
        card = fts.read(80)
    if card[:len(cardnam)] == cardnam:
        return _CRDSPL(card)
    else:
        return


def read_dis(dis_pth, ord_num, maxpnt=30):
    """Read .dis file with dispersion curves.
    @param maxpnt: maximum number of points in one curve
    @return: list with dispersion curve for each order
    """
    fsz = 4 #calcsize('f')
    dis = open(dis_pth, 'rb')
    data = []
    for j in xrange(ord_num):
        #number of points in current curve
        curlen = int(unpack('f', dis.read(fsz))[0])
        data.append([])
        for i in xrange(curlen):
            xlm = unpack('ff', dis.read(fsz*2))
            data[j].append(xlm)
        shift = (maxpnt - curlen) * fsz * 2
        dis.seek(dis.tell()+shift)
    dis.close()
    return data


def read_fds(fds_pth, ord_num, ord_len):
    """Read .fds file with array of wavelengths.
    For reading fds we must know about data shape of reference spectrum:
    number of orders and length of one order. Output data is an array of
    wavelengths (floating point numbers) for each point of reference spectrum.
    @param fds_pth: path to fds file
    @param ord_num: number of orders
    @type ord_num: unsigned int
    @param ord_len: length of order
    @type ord_len: unsigned int
    @return: numpy array of float wavelengths
    """
    fds = open(fds_pth, 'rb')
    fmt_str = 'i' * ord_len * ord_num
    data = unpack(fmt_str, fds.read(4*ord_len*ord_num))
    fds.close()
    data = np.reshape(np.array(data) / 10000., (ord_num, ord_len))
    return data


def get_wv_range(fds_pth, shift=None):
    """Read only first and last wavelengths from .fds file, so get wl range.
    @param fds_pth: path to fds file
    @param shift: seek to given position for reading last value.
        Default is os.path.getsize(fds_pth) - 4
    @return: two floats: first and last wavelength
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
    @param data: 2D array data (np.array or list)
    @param fds_pth: path for saving fds file
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
    @param ccm_pth: path to .ccm.txt file
    @return: dictionary with points for spline fitting. keys are number
    of orders, values are lists of tuples with points
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
        polynome points as values
    @param ccm_pth: path to .ccm.txt file
    """
    ccm = open(ccm_pth, 'wb')
    for ordnum, pnt in sorted(ccm_dct.items()):
        for pnx, pny in sorted(pnt):
            ccm.write(",".join(map(str, (ordnum, pnx, pny)))+'\r\n')
    ccm.close()


def take_orders(data, fds, lam, polypt, edglim=5, ish=0, wcut=.5, olim=2):
    """Take orders from data which contains given wavelength.
    Here we assume that orders are sorted by wavelength.
    @param data: data as 2D array
    @param fds: dispersion curve as 2D array of wavelengths
    @param lam: wavelength, probably containing in fds' wavelength range
    @param polypt: dictionary with points for flatten polynome plotting
    @param edglim: index for breaking points at the edges
    @param ish: intensity shift for flatted order
    @param wcut: wavelength shift in angstroms: break order with
        'incomplete' line with wcut<0 or take it otherwise)
    @param olim: limit for number of extracted orders
    @return: dictionary with tuples of fds[i], data[i] as values
        and order numbers as keys
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


def interp_chain(lams, dats, mf=5, k=1):
    """Interpolate chain with multiply factor mf.
    @param mf: Multiply factor: length of output arrays is calculated as
        length of input arrays multiplied by mf
    @param k: spline degree
    @return: interpolated input arrays
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
    @param lam: reference wavelength
    @param Va: spectrum shift in km/s, for example, heliocentric correction
    @param width: get piece of lams with length equal to 2*width angstroms
    @param off: spectrum shift in angstroms, for example, systemic velocity
        shift
    @param cutedg: drop given number of pixels on edges of the order
    @param vscl: scale mode. If True, transform wavelength scale to radial
        velocities scale
    @return: array of wavelengths or velocities, start and end indices
        of input array' slice
    """
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
        @param maskbad: treat values <0 as 'bad' and replace them with previous
            value
        @param bigval: Multiply spectrum values by big constant factor
            900000/max(data) instead of 32767/max(data)
        @param olen: 'order' length
        @param ovr: order's overlap length
        @return: data and fds as lists of orders with equal length
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
        @param chain: raw data as array of floats/ints
        @param lam: wavelengths as array of floats/ints
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
    @return: list of chains
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
    @return: data as list
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
    @return: list of chains
    """
    return split_array(np.arange(beg, end, step), ordlen, ovr)


class Spliner:
    """Plot polynome using a set of reference points in a given range
    @param spk: degree of spline
    """
    def __init__(self, beg=None, end=None, smooth=None, spk=3):
        self.smooth = smooth
        self.spk = spk
        self.beg, self.end = beg, end

    def splpl(self, xflat, yflat):
        self.tckp = splrep(xflat, yflat, s=self.smooth, k=self.spk)
        return self.mk_spline()

    def mk_spline(self):
        return splev(np.arange(self.beg, self.end), self.tckp)

    def flat_order(self, order, flatdots, beg=0, end=None, intlev=0):
        """Make interpolation of given order by spline based on given points.
        This is wrapper function for SciPy spline interpolation.
        See help of splrep and splev functions from scipy.interpolate.
        Spline parameters are set to default (cubic spline). If you want
        to change them, use class variables self.spk and self.sps.
        beg and end indices are used for slicing (cutting) data.
        @param order: 1D array of data
        @param flatdots: points for plotting flat spline. Must be sorted!
        @type flatdots: tuple or list of X and Y array: ([dotsx], [dotsy])
        @param intlev: intensity shift for flatten order
        @return: normalized (flatten) order, shifted on intlev
        """
        xflat, yflat = zip(*flatdots)
        if end is None:
            end = len(order)
        self.beg, self.end = beg + 1, end + 1
        try:
            splfit = self.splpl(xflat, yflat)
        except ValueError, err:
            print >> sys.stderr, "ValueError in flat_order:", err,
            print >> sys.stderr, "Maybe you should sort input arrays."
            raise SystemExit
        return order / splfit + intlev


class Spectractor(Spliner):
    """Class for manipulating set of 1D spectra: collect and unify data.
    First we read raw data and fds. Than select spectral orders, shift, cut
    and flat each of them, finally, write to resulting array.
    @param spectra: iterable array (named tuple) with spectra parameters:
        file paths, observation dates etc. These parameters must be present:
        id, pth, fdspth (path to fds), plpth (path to ccm.txt), Vr (spectrum
        shift, e.g. heliocentric correction). Others, e.g. jd, are optional
        and can be used as part of keys in resulting data.
        This array can be easily constructed from database table
    @param vscl: get spectrum chains in velocity scale if True (zero velocity
        wavelength must be present), use wavelengths scale otherwise
    @param cuts_dct: dictionary with number of pixels to cut for spectral
        orders
    @param dtype: data type of resulting array, could be list or dict
    @param edglim: number of pixels to cut at the edges of each spectral order
    @type edglim: integer
    """
    def __init__(self, spectra, vscl=True, cuts_dct=None, dtype=list,
                edglim=0, verbose=False):
        self.dtype = dtype
        self.fdata = self.dtype()
        self.spectra = spectra
        self.cuts_dct = cuts_dct
        self.vscl = bool(vscl)
        self.edglim = edglim
        self.verbose = verbose
        self.llam, self.hlam = None, None
        self.vshift = True
        self.wcut = 1
        self.mkey = lambda sp, x: ((sp.id, sp.flag, sp.pth, sp.Vr, sp.jd, x))
        self.Spl = Spliner()

    def get_raw(self, data_pth, fds_pth, pp_pth):
        """Read binary files containing data, wavelength scale (fds) and ccm
        (file with flatten polynome points).
        @param data_pth: path to spectrum data file
        @param fds_pth: path to spectrum dispersion curve file
        @param pp_pth: path to polynome points file
        @return: data array, fds array, polynome points dictionary
        """
        if data_pth.endswith(".200"):
            head, data = read_fits(data_pth)
            self.nx, self.ny = head.get('NAXIS2'), head.get('NAXIS1')
            #self.nx, self.ny = len(data), len(data[0])
        elif data_pth.endswith(".fits"):
            head, data = read_fits(data_pth, byteswap=True)
            self.nx, self.ny = head.get('NAXIS2'), head.get('NAXIS1')
        elif data_pth.endswith(".100"):
            (objname, self.nx, self.ny), data = read_bin(data_pth)
        else:
            print >> sys.stderr, "Sorry, unknown data file format in", data_pth
            return None, None, None
        fds = read_fds(fds_pth, self.nx, self.ny)
        polys_dct = read_ccmtxt(pp_pth)
        return data, fds, polys_dct

    def filler(self, key, lams, dats):
        """Fill self.fdata depending of its dtype (dictionary or list)."""
        if self.dtype is dict:
            self.fdata[key] = (lams, dats)
        elif self.dtype is list:
            self.fdata.append((key, lams, dats))

    def take_orders(self, data, fds, polypt, opts, lam=None, ish=0, cuts=None):
        """Take some orders from given spectrum depending of some conditions.
        The main conditions are lower and higher wavelength limits and presence
        of points for flatten polynome (polypt).
        Here we assume that each order has a fds sorted by wavelength.
        After selection each selected order becomes flatten and cut using
        class settings. Resulted chains are collect in self.fdata.

        Here we use an opts.Vr property - spectrum shift, for example,
        heliocentric correction. If self.vshift is True, Vr must be in
        velocity scale (e.g. km/s), wavelength scale otherwise.
        @param data: data as iterable 2D array
        @param fds: dispersion curve as 2D array of wavelengths
        @param polypt: dictionary with points for flatten polynome plotting
        @param lam: wavelength for converting into velocity scale,
            if self.vscl is True
        @param ish: intensity shift for flatted order
        """
        #if not self.nx:
            #self.nx = len(data)
        for i in xrange(self.nx):
            if self.llam and fds[i][-self.edglim] + self.wcut < self.llam:
                continue
            if self.hlam and fds[i][self.edglim] - self.wcut > self.hlam:
                continue
            if i+1 in polypt:
                pols = polypt[i+1]
            else:
                print >> sys.stderr, "polynome points dict in", opts.id, \
                    opts.flag, "has no", i+1, "order, continue"
                continue
            if cuts and i in cuts:
                cut = cuts[i]
            else:
                cut = None
            lams, beg, end = self.shcut_lams(fds[i], shift=opts.Vr, cuts=cut)
            if not beg and not end:
                print >> sys.stderr, "Zero length of chain", i, "continue"
                continue
            if self.verbose:
                print "Get", opts.id, opts.flag, "order", str(i+1), "pix num",
                print beg, end, "its wavelengths", round(lams[0], 2),
                print round(lams[-1], 2), "chain length:", end - beg
            dats = self.Spl.flat_order(data[i][beg:end], pols, beg, end,
                        intlev=ish)
            # The same:
            #dats = self.Spl.flat_order(data[i], pols, intlev=ish)[beg:end]
            if lam and self.vscl:
                # Convert to radial velocity scale: v = C * deltalam / lam
                lams = (lams - lam) * _C / lam
            key = self.mkey(opts, i)
            self.filler(key, lams, dats)

    def shcut_lams(self, lams, shift=0, cuts=None):
        """Get fds, shift and cut them depending on given shift and limits.
        If self.vshift is True, we calculate shift array dlam = (shift*lams)/C,
        where C is the speed of light.
        @param lams: dispersion curve as array of wavelengths
        @param shift: spectrum shift, e.g. heliocentric correction
        @param cuts: tuple with number of start and end pixels to cut for
            given order
        @return: fds chain of wavelengths or velocities, pixel numbers of
            limiting wavelengths
        """
        if shift and self.vshift:
            # Compute array of shifts for each pixel
            wshift = (shift * lams) / _C
            # shift dispersion curve array
            lams += wshift
        elif shift:
            lams += shift
        # Set cuts:
        if cuts:
            beg, end = cuts
        else:
            beg, end = self.edglim, self.ny-self.edglim
        # Change cuts depending on limits if it was defined:
        if self.llam:
            lind = np.where(lams>=self.llam)[0]
            try:
                beg = max(lind[0], beg)
            except IndexError:
                return None, None, None
        if self.hlam:
            hind = np.where(lams<=self.hlam)[0]
            try:
                end = min(hind[-1], end)
            except IndexError:
                return None, None, None
        if beg >= end:
            return None, None, None
        lams = lams[beg:end]
        return lams, beg, end

    def runner(self, lam=None, ish=None):
        """Iterate over the spectra."""
        for spec in sorted(self.spectra):
            data, fds, polys = self.get_raw(spec.pth, spec.fdpth, spec.plpth)
            if self.cuts_dct and spec.id in self.cuts_dct:
                cuts = self.cuts_dct[spec.id]
            else:
                cuts = None
            self.take_orders(data, fds, polys, spec, lam, cuts=cuts, ish=ish)

    def get_lims(self, lam=None):
        """Return boundary wavelengths or velocities depending on self.vscl.
        """
        if self.vscl and lam:
            llim = (self.llam - lam) * _C / lam
            hlim = (self.hlam - lam) * _C / lam
            return llim, hlim
        else:
            return self.llam, self.hlam

    def collect_chains(self, lam, width=4, shift=0, ish=0, wcut=None):
        """Collect spectra chains containing given wavelength in a given
        spectral window.
        @param lam: wavelength of interest
        @param width: half-width of spectral window in wavelength scale
        @param shift: shift of the spectral window for placing lam at the center
            of the window
        @param ish: relative intensity shift for flatted order
        @param wcut: tolerance parameter for rejecting partly matching orders.
            The meaning of the parameter is a shift of chain limiting
            wavelengths before comparing with window limits.
            So, some chains can be break with wcut<0 and can be processed
            with wcut>0. Default is 1
        @return: dictionary or list with spectra chains, values are tuples and
        keys are parameters. Type of resulting data are given in dtype
        """
        if wcut:
            self.wcut = wcut
        self.llam = lam - width - shift
        self.hlam = lam + width - shift
        self.runner(lam=lam, ish=ish)
        return self.fdata

    def get_chains(self, lam, llim=None, hlim=None, ish=0, wcut=None):
        """Collect spectra chains containing given wavelength in a given
        velocity range. At least one of the limits in velocity scale must be set
        @param lam: wavelength of interest
        @param ish: relative intensity shift for flatted order
        @param wcut: tolerance parameter for rejecting partly matching orders.
            The meaning of the parameter is a shift of chain limiting
            wavelengths before comparing with window limits.
            So, some chains can be break with wcut<0 and can be processed
            with wcut>0. Default is 1
        @return: dictionary or list with spectra chains, values are tuples and
        keys are parameters. Type of resulting data are given in dtype
        """
        if wcut:
            self.wcut = wcut
        if hlim and not llim:
            llim = - hlim
        elif llim and not hlim:
            hlim = - llim
        if llim and hlim:
            self.llam = (llim * lam)/_C + lam
            self.hlam = (hlim * lam)/_C + lam
        else:
            return
        self.runner(lam, ish=ish)
        return self.fdata

    def get_atlas(self, llam=None, hlam=None, ish=0):
        """Collect spectra chains in a given wavelength range.
        If limit is not set, get all orders in this direction.
        @param ish: relative intensity shift for flatted order
        @return: dictionary or list with spectra chains, values are tuples and
        keys are parameters. Type of resulting data are given in dtype
        """
        self.llam = llam
        self.hlam = hlam
        self.runner(ish=ish)
        return self.fdata
