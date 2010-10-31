#!/usr/bin/env python2

'''This is Spectractor module for some specific calculations.
With pyMidas we calculate heliocentric correction and
with pyephem we calculate precession.
'''

import sys
from datetime import datetime
from math import pi

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


spacer = lambda *args: " %s " % ' '.join(map(str, args))
commer = lambda arg: ",".join(map(str, arg))
convcoord = lambda cr, rest: (int(cr), int(divmod(rest*60, 1)[0]),
                               round(divmod(rest*60, 1)[1]*60, 3))


class HelVelCor:
    """Calculate heliocentric RV correction for number of observation times.
    @param params: list of observation parameters: date, time, coordinates.
            See docstrings of methods of this class for details
    @param dtype: data type of resulting array, must be list or tuple
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
        @return ra and dec as tuples
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
        @return list of (spectra) parameters, including heliocentric and
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
    @return ra, dec as tuples
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
