
'''This is Spectractor module for some specific calculations.
With pyMidas we calculate heliocentric correction and
with pyephem we calculate precession.
'''

import os, sys, shutil
from datetime import datetime
from math import pi

import numpy as np

try:
    from pymidas import midas
except ImportError, err:
    print >> sys.stderr, str(err) + ": We use ESO MIDAS and PyMidas for, eg. ",
    print "heliocentric RV correction. Get it here: http://www.eso.org/esomidas"
    sys.exit(True)

try:
    from coords.astrodate import utc2jd, jd2jyear
except ImportError, err:
    print >> sys.stderr, str(err) + ": coords module required, ",
    print >> sys.stderr, "download and install it from ",
    print >> sys.stderr, "http://www.scipy.org/AstroLib/ ",
    print >> sys.stderr, "or https://www.stsci.edu/trac/ssb/astrolib"
    # Another link: http://stsdas.stsci.edu/astrolib/coords_api/
    sys.exit(True)

try:
    from ephem import Equatorial
    #import ephem._libastro as libastro
except ImportError, err:
    print >> sys.stderr, str(err) + ": pyephem module required ",
    print >> sys.stderr, "for coordinates precessing."

try:
    import pyfits
except ImportError, err:
    print >> sys.stderr, str(err), "PyFITS is required for read and write FITS"


spacer = lambda *args: " %s " % ' '.join(map(str, args))
commer = lambda arg: ",".join(map(str, arg))
convcoord = lambda cr, rest: (int(cr), int(divmod(rest*60, 1)[0]),
                               round(divmod(rest*60, 1)[1]*60, 3))

def clocker(val, divc=1.):
    """Convert float value to tuple of (int h, int m, float s).
    Useful for coordinates and times.
    @param divc: division coefficient for val
    """
    h, rest = divmod(val/divc, 1)
    m, rest = divmod(rest * 60, 1)
    s = round(rest * 60, 2)
    return int(h), int(m), s


class HelVelCor:
    """Calculate heliocentric RV correction for number of observation times.
    @param params: list of observation parameters: date, time, coordinates.
        See docstrings of methods of this class for details
    @param dtype: data type of resulting array, must be dict or list
    @param obscoord: Observatory coordinates. Default is for 6-m telescope
        (BTA), SAO RAS
    """
    def __init__(self, params, dtype, obscoord=None, verbose=False):
        self.params = params
        self.verbose = verbose
        if not obscoord:
            obscoord = ((41, 26, 30), (43, 39, 12)) # 41.4414 43.6467
        self.obscoord = spacer(commer(obscoord[0]), commer(obscoord[1]))
        self.dtype = dtype
        self.fdata = dtype()
        self.intcoord = lambda h, m, s: (int(h), int(m), float(s))
        self.cmpbary = lambda date, time: ' '.join(map(','.join, (
                                   map(str, date), map(str, time))))

    def precess_midas(self, date, time, objcoords):
        """Wrapper for comp/prec MIDAS routine.
        @param date: tuple of year, month and day
        @param time: tuple of hours and minutes
        @return: ra and dec as tuples
        """
        minuhour = (time[0]/24.) + (time[1]/1440.)
        compspec_str = ",".join(map(str, date)) + '.' + \
                       str(int(minuhour*10000))
        ra, dec, epoch = objcoords
        coords_str = spacer(commer(ra), commer(dec), str(epoch))
        cmpr = midas.computePrec(coords_str, compspec_str)
        ra = self.intcoord(cmpr[6], cmpr[8], cmpr[10])
        dec = self.intcoord(cmpr[13], cmpr[15], cmpr[17])
        #jyear_epoch = float(cmpr[3])
        return ra, dec

    def compbary(self, date, time, ra, dec):
        """Wrapper for comp/bary MIDAS routine.
        @param date: tuple of year, month and day
        @param time: tuple of hours and minutes
        @param ra: tuple of coords
        @param dec: tuple of coords
        @return: heliocentric and barycentric radial velocity corrections
        """
        compbary_str = self.cmpbary(date, time)
        coords_str = spacer(','.join(map(str, ra)), ','.join(map(str, dec)))
        cmbr = midas.computeBary(compbary_str, coords_str, self.obscoord)
        hcor = round(float(cmbr[20]), 4)
        bcor = round(float(cmbr[14]), 4)
        return hcor, bcor

    def filler(self, ID, jd, jyr, ra, dec, brcor, hlcor):
        """Dummy function for filling resulting array"""
        if self.dtype is dict:
            self.fdata[ID] = (jyr, round(jd-2400000.5, 5), (ra, dec),
                brcor, hlcor)
        elif self.dtype is list:
            self.fdata.append((ID, jyr, round(jd-2400000.5, 5), (ra, dec),
                brcor, hlcor))

    def getdatime(self, spec):
        """Convert datetime of observetion to JD. Pay attention to indices!
        """
        if type(spec[0]) is not type(datetime(1970, 1, 1)):
            date, time = spec[:3], spec[3:5]
            jd = utc2jd(datetime(*spec[:5]))
        else:
            date, time = spec[0].timetuple()[:5]
            jd = utc2jd(spec[0])
        return date, time, jd

    def veltimes(self, objcoords, midprecess=False):
        """Compute heliocentric correction for the set of observation times.
        @param objcoords: Object coordinates with its epoch as tuple:
            (ra, dec, epoch). ra and dec are also tuples.
        @param midprecess: Precess coordinates using ESO MIDAS comp/bary
            routine, use pyephem instead
        @return: list of (spectra) parameters, including heliocentric and
            barycentric radial velocity correction
        """
        for ID, spec in sorted(self.params.items()):
            date, time, jd = self.getdatime(spec)
            jyr = round(jd2jyear(jd), 5)
            if midprecess or abs(objcoords[2]-2000) > .01:
                # Precess coords to jyear_epoch using MIDAS COMP/PREC
                ra, dec = self.precess_midas(date, time, objcoords)
            elif not midprecess:
                # Precess coords using pyephem
                ra, dec = precess(jd-2400000.5, objcoords, verbose=self.verbose)
                if self.verbose:
                    print "Precessed coordinates:", ra, dec
            helrvcor, barrvcor = self.compbary(date, time, ra, dec)
            self.filler(ID, jd, jyr, ra, dec, barrvcor, helrvcor)
        return self.fdata

    def velcoords(self, midprecess=False):
        """Compute heliocentric correction for the set of coordinates.
        """
        for ID, spec in sorted(self.params.items()):
            date, time, jd = self.getdatime(spec)
            jyr = round(jd2jyear(jd), 5)
            objcoords = spec[5:8]
            if not objcoords[2]:
                objcoords[2] = jyr
            if midprecess or abs(objcoords[2]-2000) > .01:
                # Precess coords to jyear_epoch using MIDAS COMP/PREC
                ra, dec = self.precess_midas(date, time, objcoords)
            elif not midprecess:
                # Precess coords using pyephem
                ra, dec = precess(jd-2400000.5, objcoords, verbose=self.verbose)
                if self.verbose:
                    print "Precessed coordinates:", ra, dec
            helrvcor, barrvcor = self.compbary(date, time, ra, dec)
            self.filler(ID, jd, jyr, ra, dec, barrvcor, helrvcor)
        return self.fdata


def precess(mjd, objcoords, verbose=False):
    """Precess coords using pyephem from J2000 to mjd.
    @param mjd: Modified julian date: JD - 2400000.5
    @param objcoords: Object coordinates with its epoch as tuple
    @return: ra, dec as tuples
    """
    # Dublin Julian Day = JD - 2415020.0
    #J2000 = 2451545.0 - 2415020.0  <-- it is from libastro sources
    # Maybe we should try direct call to libastro:
    #libastro.precess(2000, 2007.794, 1.386,  0.282)
    #2415020.0 - 2400000.5 = 15019.5
    mjd = mjd - 15019.5
    if verbose:
        print "We get the coords", objcoords
    # coords as strings, separated by ":"
    dv = lambda arg: ":".join(map(str, arg))
    coord = Equatorial(dv(objcoords[0]), dv(objcoords[1]))
    precessed_coords = Equatorial(coord, epoch=mjd)
    ra, dec = precessed_coords.to_radec()
    # Get tuple for each coordinate
    ra, dec = convcoord(*divmod(ra*12/pi, 1)), convcoord(*divmod(dec*180/pi, 1))
    if verbose:
        print "Precess to", mjd, ra, dec, "with J2000 coords", objcoords
    return ra, dec


class MidWorker:
    """Class with simple wrappers for some ESO MIDAS routines
    """
    def __init__(self, usedisp=False, verbose=False):
        self.verbose = bool(verbose)
        self.usedisp = usedisp
        self.ext = 'bdf'
        self.med = 35
        self.avermed = lambda m: ",".join(('median', str(m), str(m), 'data'))

    def plot(self, bnam, plt="row", cut=2010, over=False):
        """Wrapper for Midas functions plot/row, plot/colu"""
        if over:
            over = "over"
        else:
            over = ""
        if self.verbose:
            print "Plot", over, plt, bnam, "cut on pixel", cut
        cmd = " ".join((over+"plot/"+plt, bnam, "@"+str(cut/2)))
        midas.do(cmd)

    def readesc(self, fname, desc="npix"):
        """Wrapper for read/desc Midas procedure"""
        nx, ny = midas.readDesc(fname, desc)[-2:]
        return int(nx), int(ny)

    def loima(self, img, chan=0):
        """Load image with autocuts"""
        if not (os.path.isfile(img) or os.path.isfile(img+'.'+self.ext)):
            return None, None
        nx, ny = self.readesc(img)
        dx, dy = map(int, midas.writeOut('{ididev(2)} {ididev(3)}'))
        if nx > dx:
            i = round(- (nx/(dx - 10.)) - 1, 2)
        if ny > dy:
            j = round(- (ny/(dy - 10.)) - 1, 2)
        scale = ",".join((str(i), str(j), 'a'))
        if self.verbose:
            print "Scale to", scale
        midas.loadImag(img, str(chan), 'scale='+scale, 'center=c,c')
        midas.writeKeyw("log/i/4/1", "2")
        #midas.putKeyword("log", 2)
        midas.statistImag(img)
        midas.do('@a autocuts ' +img+ ' SIGMA')
        midas.writeKeyw("log/i/4/1", "0")
        return i, j

    def savebdf(self, name, bdf_pth, wridesc=True):
        """Save a MIDAS bdf image, using existing FITS file
        Also write descriptors, load it on display
        """
        midas.indiskFits(name, bdf_pth)
        if wridesc:
            midas.writeDesc(bdf_pth, 'start/d/1/2 1.,1.')
            midas.writeDesc(bdf_pth, 'step/d/1/2 1.,1.')
        if self.usedisp:
            self.loima(bdf_pth)

    def crebias(self, files, biasname="bias", med=None, delfits=False):
        if not med:
            med = self.med
        catname = biasname+".cat"
        midas.createIcat(catname, "NULL")
        for file_pth in files:
            midas.addIcat(catname, file_pth)
        midas.averageImag(biasname, "=", catname, "? +", self.avermed(med))
        if self.usedisp:
            self.loima(biasname)
        midas.outdiskFits(biasname, os.path.join('bias.fits'))

    def creflat(self, files, flatname="flat", delcat=False, averopt=""):
        """Create flat field average image
        @param delcat: delete created cat, if icat is present
        @param pfix: postfix of image filename, the "d" value means
                dark substracted
        """
        catname = flatname+".cat"
        midas.createIcat(catname, "NULL")
        for file_pth in files:
            midas.addIcat(catname, file_pth)
        midas.averageImag(flatname, "=", catname, "?", averopt)
        if self.usedisp:
            self.loima(flatname)

    def filtcosm(self, specname, ns=2, cosm="cosm"):
        """Wrapper for filter/cosm Midas routine, which removes cosmic
        ray events from CCD image and replace them by a local median value.
        @param specname: input image name. Resulting image will be with the
                same name and the 'f' postfix
        @param ns: threshold for the detection of cosmic rays
        param cosm: frame containing mask for the giving image
        """
        #if specname.endswith(".bdf"):
        resnam = specname.split(".")[0]+"f"
        params = ",".join(map(str, (0, 2, 8, ns, 2)))
        if self.usedisp:
            i, j = self.loima(specname)
        midas.filterCosm(specname, resnam, params, cosm)
        if self.usedisp:
            scale = ",".join((str(i), str(j), 'max'))
            midas.loadImag(cosm, "scale="+scale, "cuts=0,1")

    def cycle(self, beg, end, plt, fcat=False, pfix=""):
        k = True
        if fcat:
            midas.createIcat(fcat, "NULL")
        else:
            namlst = []
        for bnum in xrange(beg, end+1):
            bnam = 'l'+str(bnum).zfill(3) + pfix
            if os.path.isfile(bnam+'.'+self.ext):
                if k:
                    if self.verbose:
                        print "Add first", bnam, "image"
                    nx, ny = self.readesc(bnam)
                    cut = ny if plt == "row" else nx
                    self.plot(bnam, plt, cut)
                    firstnam, k = bnam, False
                else:
                    if self.verbose:
                        print "Add", bnam, "image"
                    self.plot(bnam, plt, cut, over=True)
                if fcat:
                    midas.addIcat(fcat, bnam)
                else:
                    namlst.append(bnam)
        if not fcat:
            namlst = ",".join(namlst)
        else:
            namlst = fcat
        return firstnam, namlst, cut

    def pltgra(self, resnam, plt, cut, over=True, midlo=False):
        midas.setGrap("color=2")
        self.plot(resnam, plt, cut, over=over)
        midas.setGrap("color=1")
        if midlo:
            midas.do("@@ lo_ima " + resnam)
        elif self.usedisp:
            self.loima(resnam)

    def crebiasref(self, beg, end, resnam='bias', icat=False, plt="row",
                med=None):
        """Create bias frame from set of biases.
        @param resnam: resulting frame name
        @param med: value of low and high interval in median averaging
        @param plt: plotting mode, could be row or col.
        """
        firstnam, namlst, cut = self.cycle(beg, end, plt, fcat=icat)
        if not med:
            med = self.med
        midas.averageImag(resnam, "=", namlst, "? +",  self.avermed(med))
        self.pltgra(resnam, plt, cut)

    def creflatref(self, beg, end, resnam=None, icat='flat.cat', plt="col",
                pfix="d", averopt="", delcat=False):
        """Create flat field average image
        @param delcat: delete created cat, if icat is present
        @param pfix: postfix of image filename, the "d" value means
                dark substracted
        """
        firstnam, namlst, cut = self.cycle(beg, end, plt, fcat=icat, pfix=pfix)
        if not resnam:
            resnam = firstnam.replace('d', 'ad')
        midas.averageImag(resnam, "=", namlst, "?", averopt)
        self.pltgra(resnam, plt, cut)
        if delcat and icat:
            os.remove(icat)

    def subdark(self, beg, end, biasname="bias", pfix="d"):
        """Substract bias from images"""
        for bnum in xrange(beg, end+1):
            bnam = 'l'+str(bnum).zfill(3)
            if os.path.isfile(bnam+'.'+self.ext):
                nameo = bnam + pfix
                midas.computeImag(nameo, "=", bnam, "-", biasname)
                midas.statistImag(nameo)
                if self.usedisp:
                    self.loima(nameo)
                midas.copyDd(bnam, "*,3", nameo)
                os.remove(bnam+'.'+self.ext)


class Preparer(MidWorker):
    """Class for image processing with MIDAS, against MIDAS sripting language =)
    We use MIDAS routines through pyMidas, but replace it with numpy and pyfits,
    if it is possible.

    This code also heavily depends of some specific characteristics of SAO RAS
    spectral archive, mostly on spectra taken on Nesmyth Echelle Spectrograph of
    6-meter telescope.
    It is not pipeline package, but a set of non-interactive methods for packet
    processing. It contains methods for averaging images, substracting bias
    frame etc.
    @param dest_dir: destination dir for moving processed fits files. Could be
        omitted.
    @param usedisp: use Midas display for loading images using Midas routines
    """
    def __init__(self, dest_dir="arch", usedisp=False, verbose=False):
        self.dest_dir = dest_dir
        self.verbose = bool(verbose)
        self.usedisp = usedisp
        self.fixfitslog = {True: 'fix', False: 'silentfix'}[self.verbose]
        self.print_str = "[01;31m%s not found in header! Set it to zero.[0m"
        # Define keys, presented in FITS headers:
        self.stdkeys = ['AUTHOR', 'REFERENC', 'DATE', 'ORIGIN', 'DATE-OBS',
            'TELESCOP', 'INSTRUME', 'OBSERVER', 'OBJECT', 'EQUINOX', 'EPOCH']
        self.nonstdkeys = ['OBSERVAT', 'CHIPID', 'DETNAME', 'CCDTEMP',
                           'P_DEWAR', 'SEEING', 'PROG-ID']
        self.otherkeys = ['WIND', 'T_OUT', 'T_IN', 'T_MIRR', 'PRES',
                           'ST', 'MT', 'FOCUS']
        self.journal_dct = {}
        self.obs_dct = {}
        self.sep, self.postfix, self.ext = '_', '.', 'fits'
        self.key_dct = {"bias": 'b', "thar": 't', "flatfield": 'f', "ff": 'f',
                        "ll": 'f', "sky": 's'}
        #self.imgcut_dct = {2052: (4, 1992, 4, 2039), # NES 2052x2052 chip
        self.imgcut_dct = {2052: (3, 2048, 4, 2039), # NES 2052x2052 chip
        # Resulting shape after cut for NES images with cuts (4, 1992, 4, 2039)
        # will be (2035, 1988). NES with full use: (3, 2048, 0, 2046),
        # shape will be (2045, 2046)
        2068: (0, 2048, 2, 2048), # LPR 2068x2072 chip
        1050: (5, 1040, 2, 1157) # PFES/Lynx 1050x1170 chip
        }
        self.keyset = lambda key: self.key_dct.get(key.lower()) \
                if self.key_dct.get(key.lower()) else 'o'
        self.trygetkey = lambda head, *args: filter(head.has_key, args)
        self.seper = lambda fstr, h, m, s: \
                "%s" % ":".join(["%02d" % h, "%02d" % m, fstr % s])
        self.nesbadrows = [(338, 337, 339), (568, 567, 569),
        (336, 336, 335), (335, 336), (1044, 1044, 1043, 1042),
        (1043, 1044), (1042, 1044), (433, 433, 434, 435),
        (434, 433), (435, 433), (185, 185, 184), (184, 185), (1903, 1903, 1902),
        (1902, 1903), (1229, 1229, 1228), (1228, 1229), (47, 47, 48), (48, 47)]
        self.rowdoc = lambda d, xs: sum([d[:, 2052-r-60] for r in xs]) / len(xs)
        self.mw = MidWorker(usedisp=usedisp, verbose=verbose)
        self.biases, self.flats = {}, {}, {}

    def prepare(self, files, logonly=False, crop=True, rot=True, badrow=False,
            substbias=False, cleanhead=True, filmove=True, process=True):
        """Method for packet preparing and/or logging of raw images,
        mostly taken on spectrographs of 6-meter telescope.
        @param files: list of raw image parameters with its paths
        @param logonly: create only log, without image processing.
        @param crop: crop each image depending on its shape, which is
            determine used CCD chip
        @param rot: rotate image clockwise
        @param badrow: mask bad rows on NES CCD chip
        @param cleanhead: remove empty strings in fits headers
        @param filmove: move or not source files (raw images) to
                self.arch_dir, after they are have been processed
        """
        if substbias and type(substbias) is str:
            self.bias = pyfits.getdata(substbias)
            self.substbias = True
        else:
            self.substbias = False
        self.crop = crop
        self.rot = rot
        self.badrow = badrow
        if files and filmove and not os.path.isdir(self.dest_dir):
            os.mkdir(self.dest_dir)
        for (file_pth, num, date, fflag, keymode) in sorted(files):
            if self.verbose:
                print pyfits.info(file_pth)
                print "Open file", file_pth
            curfts = pyfits.open(file_pth) #mode="update"
            #head = curfts[0].header
            self.logger(curfts, num, date, file_pth, keymode=keymode)
            if logonly:
                continue
            if cleanhead:
                # Clean header - delete empty strings:
                del curfts[0].header['']
                #head.add_blank(after='SHSTAT')
            # Write header to FITS object (not to disk!)
            #curfts[0].header = head
            if file_pth in self.biases:
                if filmove:
                    shutil.move(file_pth, self.dest_dir)
                continue
            elif file_pth in self.flats:
                newname = self.flats[file_pth]
                hilim = 9000
            else:
                date = "".join((str(date[0]), str(date[1]).zfill(2),
                                              str(date[2]).zfill(2)))
                newname = os.path.join(self.dest_dir, "".join((date, self.sep,
                            num, self.flag, self.postfix, self.ext)))
                hilim = 2600
            if process:
                self.processdata(curfts, hilim=hilim)
            self.savefits(curfts, newname)
            if filmove:
                shutil.move(file_pth, self.dest_dir)
            if file_pth in self.flats:
                continue
            self.mw.savebdf(newname, "l"+num+"d.bdf")
            if self.verbose:
                print pyfits.info(newname)

    def crebias(self, biasname="bias", files=None, filmove=True, log=True):
        if files:
            self.prepare(files, process=False, cleanhead=False,
                        log=log, filmove=filmove)
        elif not self.biases:
            print "No biases, return"
            return
        for file_pth, newname in self.biases.items():
            if not os.path.isfile(newname):
                shutil.copy2(file_pth, newname)
        self.mw.crebias(self.biases.values(), biasname=biasname)

    def creflat(self, flatname="flat"):
        if self.flats:
            self.mw.creflat(self.flats.values(), flatname=flatname)
            self.mw.filtcosm(flatname)

    def processdata(self, curfts, hilim=600):
        """Process image data: substract bias, crop, rotate, mask bad rows"""
        data = curfts[0].data
        dshape = data.shape
        if self.substbias:
            data -= self.bias
        if self.badrow:
            self.nesbadrow(data, hilim=hilim)
        if self.crop and dshape[0] in self.imgcut_dct:
            x1, x2, y1, y2 = self.imgcut_dct[dshape[0]]
            #midas.extractImag({name} = {namef}[@{x1},@{y1}:@{x2},@{y2}])
            data = data[y1:y2, x1:x2]
        if self.rot and dshape[0] >= 2052:
            data = np.rot90(data)
        #if self.flip:
            ## flip for Lynx/Pfes
            #data = np.flipud(data)
        # Write result
        curfts[0].data = data
        curfts[0].scale('float32', 'old') # int16

    def savefits(self, curfts, newname):
        """Save current fits"""
        if os.path.isfile(newname):
            bcknam = newname.replace("."+self.ext, "-backup."+self.ext)
            shutil.move(newname, bcknam)
        try:
            curfts.writeto(newname)
        except pyfits.core.VerifyError:
            curfts.verify(self.fixfitslog)
            curfts.writeto(newname)
        curfts.close()

    def nesbadrow(self, d, lowcut=-25, hilim=600):
        """Mask bad pixels in raw images taken with Uppsala NES chip.
        """
        for row in self.nesbadrows:
            d[:, 2052-row[0]-60] = self.rowdoc(d, row[1:])
        d[:, 1992] = (d[:, 1991] + d[:, 1993]) / 2
        d[:, 1997] = (d[:, 1996] + d[:, 1998]) / 2
        d[:, 2005] = (d[:, 2004] + d[:, 2006]) / 2
        d[:, 2008] = (d[:, 2007] + d[:, 2009]) / 2
        d[:, 2010] = (d[:, 2009] + d[:, 2011]) / 2
        d[:, 2014] = (d[:, 2013] + d[:, 2015]) / 2
        d[:, 2027] = (d[:, 2026] + d[:, 2028]) / 2
        for i in xrange(2036, 2038):
            d[i, :22] = (d[i-1, :22] + d[i-2, :22])/2
            d[i, 1337:1569] = (d[i-1, 1337:1569] + d[i-2, 1337:1569])/2
        for i in xrange(2038, 2040):
            d[i, :24] = (d[i-1, :24] + d[i-2, :24])/2
            d[i, 1337:1569] = (d[i-1, 1337:1569] + d[i-2, 1337:1569])/2
        for i in xrange(2040, 2042):
            d[i, :26] = d[i-1, :26]
            d[i, 1125:1692] = (d[i-1, 1125:1692] + d[i-2, 1125:1692])/2
        for i in xrange(2042, 2044):
            d[i, :34] = d[i-1, :34]
            #d[i, 1087:27] = d[i-1, 1087:27]
        #for i in xrange(2044, 2046):
            #d[i, :2052-2005] = d[i-1, :2052-2005]
            #d[i, 2052-1130:2052-15] = d[i-1, 2052-1130:2052-15]
        #zonlim = (0, 250, 1979) #[31:96, 1988:2040]
        zonlim = (1, 118, 1988)
        zone = d[zonlim[0]:zonlim[1], zonlim[2]:]
        med = min(np.median(zone), hilim)
        d[zonlim[0]:zonlim[1], zonlim[2]:][np.where(zone > hilim)] = med
        #d[:280, 2004:2006][np.where(d[:280, 2004:2006] > hilim)] = med
        d[:280, 2004:2011][np.where(d[:280, 2004:2011] > hilim)] = med
        d[np.where(d < lowcut)] = lowcut

    def modkey(self, val, keymode=False):
        """Modify input value, depending on keymode.
        Because of different standards of writing FITS headers in SAO RAS,
        we need to extract the right value. For example, val maybe string,
        float(string) or int(float * const).
        @param keymode: False (None), 'flstr' or 'str'. Define, how we can
            extract correct information from given value
        """
        if val and not keymode and 3000 < abs(val) < 1250000:
            val = val / 3600.
        elif keymode == "flstr": # or abs(val) > 1250000:
            val = str(val)
            val = float(val[-4:])/3600 + int(val[-6:-4])/60. + int(val[:-6])
        elif keymode == "str":
            if type(val) is str and " " in val:
                val = val.split()
            elif type(val) is str and ":" in val:
                val = val.split(":")
            val = float(val[2])/3600 + int(val[1])/60. + int(val[0])
        return val

    def pickhead(self, head, keymode, *args):
        """Look into FITS header for existance of keys, return head[key].
        First, we look on existence of keys, containing in args.
        Second, we get existing value from header or return zero.
        """
        key = self.trygetkey(head, *args)
        if key:
            val = head.get(key[0])
        else:
            print self.print_str % " ".join(args)
            return 0
        if keymode == "raw":
            return val
        else:
            return self.modkey(val, keymode=keymode)

    def getcoords(self, head, keymode=False):
        ra = self.pickhead(head, keymode, 'RA')
        if not keymode:
            ra = ra / 15.
        dec = self.pickhead(head, keymode, 'DEC')
        az  = self.pickhead(head, keymode, 'A_PNT', 'AZIMUTH')
        zen = self.pickhead(head, keymode, 'Z_PNT', 'ZDIST', 'ZENDIST')
        if type(az) in (int, float) and az < 0:
            az += 360
        #sidertime_start = self.pickhead(head, keymode, 'ST')
        ra, dec = clocker(ra), clocker(dec)
        az, zen = clocker(az), clocker(zen)
        return ra, dec, az, zen

    def logger(self, curfts, num, date, file_pth, keymode=False):
        """Create an observation log for input images.
        Search for some specific keys in given fits header and write it
        to 'journal': journal_dct and obs_dct dictionaries.
        """
        # Read coords
        try:
            ra, dec, az, zen = self.getcoords(curfts[0].header, keymode=keymode)
        except (ValueError, TypeError, AttributeError):
            print "Trying to fix", file_pth, "FITS header."
            curfts.verify(self.fixfitslog)
            keymode = "str"
            ra, dec, az, zen = self.getcoords(curfts[0].header, keymode=keymode)
        head = curfts[0].header
        #except (TypeError, AttributeError):
            #ra, dec, az, zen = (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)
        # 199xxxxx(DVD01) 20030810(DVD07)--20080331(DVD31), since 20081103_002
        wind = self.pickhead(head, 'raw', 'WIND')
        # This block in 20030810(DVD07)--20080331(DVD31), since 20081103_002
        moscowtime_start = self.pickhead(head, 'raw', 'MT') / 3600.
        focus = self.pickhead(head, 'raw', 'FOCUS')
        t_out = self.pickhead(head, 'raw', 'T_OUT')
        t_in = self.pickhead(head, 'raw', 'T_IN')
        t_mirr = self.pickhead(head, 'raw', 'T_MIRR')
        pres = self.pickhead(head, 'raw', 'PRES')
        # HUMD in 20070902(DVD28)--20080331(DVD31), since 20081103_002
        hum = self.pickhead(head, 'raw', 'HUMD')
        # JD since 20080422(DVD31)
        jd = head.get('JD')
        # keymode is here!
        self.journal_dct['KEYMODE'] = [keymode]
        for key in (self.nonstdkeys + self.stdkeys):
            if key not in self.journal_dct and head.get(key):
                self.journal_dct[key] = [head.get(key),]
            elif head.get(key) and head.get(key) not in self.journal_dct[key]:
                self.journal_dct[key].append(head.get(key))
            elif key not in self.journal_dct and not head.get(key):
                self.journal_dct[key] = []
        # Block of standard "observation" keys
        time_end = self.pickhead(head, 'raw', 'TM_END') / 3600.
        try:
            time_start = self.pickhead(head, 'raw', 'TM-START',
                                                    'TM_START') / 3600.
            time = clocker(time_start)
        except TypeError:
            time = self.pickhead(head, 'raw', 'TM-START', 'TM_START')
        self.exp = int(self.pickhead(head, 'raw', 'EXPTIME'))
        self.obj = self.pickhead(head, 'raw', 'OBJECT')
        self.flag = self.keyset(self.obj)
        if self.flag == "b" and self.exp != 0:
            self.flag = "o"
        elif self.exp == 0 and self.flag != "b":
            self.flag = "b"

        if self.flag == "b":
            newname = "".join(('b', num, self.postfix, self.ext))
            self.biases[file_pth] = newname
        elif self.flag == "f":
            newname = "".join(('f', num, self.postfix, self.ext))
            self.flats[file_pth] = newname

        dateobs = head.get('DATE-OBS')
        if self.verbose:
            self.warner(self.exp, self.obj, num)
        self.obs_dct[int(num)] = [date, self.flag, self.obj, time, self.exp, ra, dec,
                az, zen, dateobs]

    def warner(self, exp, obj, num):
        """Warn about some problems."""
        if exp == 0 and str(obj).lower() != "bias":
            print "Wrong object name:", obj, "with exp=0 in spectrum", num
        #if head.get('CCDTEMP') != None and -103 > head.get('CCDTEMP') > -97:
            #print "incorrect CCDTEMP", ccdtemp, "in spectrum number", num

    def savelog(self, filenam="night", wdir=""):
        """Save observation log to file."""
        if os.path.isfile(os.path.join(wdir, filenam + ".log")):
            shutil.move(os.path.join(wdir, filenam + ".log"),
                        os.path.join(wdir, filenam + ".old"))
        obslog = open(os.path.join(wdir, filenam+".log"), "w")
        for num, val in self.obs_dct.items():
            obsstr = self.logstringer(num, val)
            print >> obslog, obsstr
        obslog.close()

    def savejournal(self, filenam="journal", wdir=""):
        """Save observation journal to file."""
        journpth = os.path.join(wdir, filenam + ".log")
        if os.path.isfile(journpth):
            shutil.move(journpth, os.path.join(wdir, filenam + ".old"))
        journ = open(journpth, 'w')
        for key, val in sorted(self.journal_dct.items()):
            #print "{0:8s}: {1:s}".format(key, ", ".join(map(str, sorted(val))))
            if not val or val == [" "]:
                continue
            elif key in ("CCDTEMP", "P_DEWAR") and val:
                val = [str(np.mean(val))]
            elif key in ("DATE-OBS", "ORIGIN", "AUTHOR", "OBSERVER",
                         "SEEING") and val:
                val = [", ".join(map(str, val))]
            print >> journ, "%9s %s" % (key, " ".join(map(str, val)))
            if self.verbose:
                print "%8s: %s" % (key, ", ".join(map(str, sorted(val))))
        journ.close()

    def logstringer(self, num, val):
        """Make a string for given parameters of observation log."""
        date, flag, obj, time, exp, ra, dec, az, zen, dateobs = val
        date = "".join((str(date[0]), str(date[1]).zfill(2),
                                     str(date[2]).zfill(2)))+" "+flag
        #date = "".join(map(str, date))+" "+flag
        time = self.seper("%02d", *time)
        ra = self.seper("%04.1f", *ra)
        dec = self.seper("%04.1f", *dec)
        az = self.seper("%04.1f", *az)
        zen = self.seper("%04.1f", *zen)
        return "%3d %13s %15s %8s %4d %11s %11s %11s %10s %s" % (int(num),
                date, obj, time, exp, ra, dec, az, zen, dateobs)
