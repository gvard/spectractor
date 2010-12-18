
'''Spectractor submodule containing functions for managing various files,
journals, creating logs etc.

@author: Dmitry Nasonov
'''

import re
import os, sys, shutil, glob
import cPickle
import collections

from spectractor import get_wv_range
data_dir, out_dir = 'data', 'out'


def filler(indct, key, val):
    """Fill dictionary which values are lists."""
    try:
        if val not in indct[key]:
            indct[key].append(val)
    except KeyError:
        indct[key] = [val, ]
    except TypeError:
        print >> sys.stderr, "TypeError in filler! Key, value:", key, val
        sys.exit(True)


def fill_dctset(indct, key, val):
    """Fill dictionary which values are sets."""
    try:
        indct[key].add(val)
    except KeyError:
        indct[key] = set([val])


def backuper(file_pth, verbose=False):
    """Backup file if it is exist: move it with '-backup' postfix."""
    if os.path.isfile(file_pth):
        shutil.move(file_pth, file_pth + '-backup')
        if verbose:
            print "Move", file_pth, "to", file_pth + '-backup'


def rexper(regexp, path, reflag=None, splitnig=False):
    """Get parameters from file name using regexp.
    Parameters are: int ID (for example, night number), flag.
    For now, regexp.findall must return 3 values! If not, return (None, None).
    @param regexp: regular expression describing file name
    @param path: path to file
    @param reflag: tuple of two values. If not empty, flag containing in first
        value will be replaced to second one
    @param splitnig: split ID to tuple of 'night' and spectrum number
    @return: ID (night number): int or tuple of ints with splitnig, flag
    """
    try:
        bname = os.path.basename(path)
        night, img, flag = regexp.findall(bname)[0]
    except (IndexError, ValueError):
        print >> sys.stderr, "Filename", bname, "does not match regexp!"
        return None, None
    if reflag and flag in reflag[0]:
        flag = reflag[1]
    if not splitnig:
        night = int(night + img.zfill(3))
    else:
        night = (int(night), int(img))
    return night, flag


class Walker:
    """Directory walker: search files and dirs in args matching regexp.
    @param exts: set of file extensions for matching
    @param namrxp: filename regular expression, but without extension
    @param args: list of paths to files and dirs
    """
    def __init__(self, args, exts=("200", "fds"), namrxp=".*", verbose=False):
        self.files_dct = {}
        self.workdirs = []
        self.path_set, self.spath_set = set(), set()
        self.verbose = verbose
        # rexp_spec contains 3 groups: night number, spectrum number and flag.
        self.rexp_spec = "^[a-z_]{0,3}(\d{3})(\d{1,3})[-_]?([a-z_\-]*).*\."
        self.fdsrexp = re.compile(self.rexp_spec+'fds'+'$', re.I | re.U)
        # rexp_mt and rexp_fits contain groups for date (year, month, day),
        # spectrum number and flag.
        self.rexp_mt = \
                "^[lns](\d{1})([A-C0-9])(\d{2})[_-]?(\d{2,3})([a-z_]?).*"
        self.rexp_fits = \
                "^(?:Bn)?(\d{4})(\d{2})(\d{2})[_-]?(\d{3})([spbfdot_]?).*"
        exts = "|".join(exts)
        self.rxp = re.compile(namrxp+'\.('+exts+')$', re.I | re.U)
        # Pathfinder: find dirs and files in args, append them.
        for arg in args:
            if os.path.isdir(arg):
                self.workdirs.append(arg)
                if verbose:
                    print "Appending", arg, "to the list of working dirs"
            elif os.path.isdir(os.path.join(os.getcwd(), arg)):
                self.workdirs.append(os.path.join(os.getcwd(), arg))
                if verbose:
                    print "Appending", os.path.join(os.getcwd(), arg), \
                          "to the list of working dirs"
            elif os.path.isfile(arg) and len(self.rxp.findall(arg)):
                filler(self.files_dct, self.rxp.findall(arg)[0], arg)
                if verbose:
                    print "Appending", arg, "to the list of working files"
            elif verbose:
                print arg, "is not a directory or needed file, continue."

    def dirwalker(self):
        """Walk on workdirs, collect subdirectories.
        @return: list of subdirectories
        """
        if not self.workdirs:
            print "Sorry, no dirs for walk."
            return
        subdirs = []
        for dirw in self.workdirs:
            # dir path, inner dirs, inner filenames
            for root, dirs, files in os.walk(dirw):
                for wdir in sorted(dirs):
                    subdirs.append(os.path.join(root, wdir))
        return subdirs

    def walk(self, exclude="excluded"):
        """Walk into working dirs, search some kind of files.
        @param exclude: if it in file path, exclude file. Useful when
            we want to exclude some subdir
        """
        if not self.workdirs:
            if not self.files_dct:
                print "Sorry, no files found."
            return
        for dirw in self.workdirs:
            # dir path, inner dirs, inner filenames
            for root, dirs, files in os.walk(dirw):
                #sorted(dirs + files) == sorted(os.listdir(root))
                for wfile in sorted(files):
                    if len(self.rxp.findall(wfile)) and exclude not in root:
                        filler(self.files_dct, self.rxp.findall(wfile)[0],
                               os.path.join(root, wfile))

    def resolver(self, ext='rv', filtr=None, reflag=None, splitnig=False):
        """Select data, get ID and 'flag' from filename, collect paths.
        Here we assume, that file name contains information about ID
        (night number, number of spectrum) and flag. Others will be broken
        (or you must change default regular expression).
        @param ext: select file type, must be file extension
        @param filtr: selecting only given flags, exclude others. Select all
            when filtr in ("all", None)
        @param reflag: tuple of two values. If not empty, flag containing in
            first value will be replaced to second one
        @param splitnig: split ID to tuple of 'night' and spectrum number
        @return: list of tuples with 'night', flag and path to file
        """
        try:
            files = self.files_dct[ext]
        except KeyError:
            print >> sys.stderr, "Sorry,", ext, "files seem not to exist."
            return self.path_set
        rexp = re.compile(self.rexp_spec + ext + '$', re.I | re.U)
        for path in files:
            night, flag = rexper(rexp, path, reflag=reflag, splitnig=splitnig)
            # Collect "kosher" paths:
            if filtr in ("all", None) or flag in filtr:
                self.path_set.add((night, flag, path))
            elif self.verbose:
                print "Exclude", night, path
        if not len(self.path_set):
            print >> sys.stderr, "No", ext, "files in path!"
        return self.path_set

    def select_spectra(self, journal, ext='fits', iflag=None, splitnig=False):
        """Select data, get ID ('night') and flag from filename, collect paths.
        Here we assume, that file name contains information about ID
        (night number, number of spectrum) and flag. Others will be broken
        (or you must change default regular expression). Then we check presence
        of ID in journal keys, get spectrum parameters from journal and append
        this data to resulting list.
        @param journal: journal dictionary with spectra IDs as keys and spectra
            parameters (e.g. JYear, MJD, spectrograph key, limiting
            wavelengths, S/N etc) as values
        @param ext: select file type, must be file extension
        @param iflag: Flags to select. Can be a string, if flag is one symbol
        @param splitnig: split 'night' to tuple of night and spectrum number
        @return: list of tuples with ID ('night'), flag, path to file and
            other parameters. Currently parameters are fds path, dictionary
            of points, defining the flatten polynome, heliocentric
            correction and S/N value
        """
        try:
            files = self.files_dct[ext]
        except KeyError:
            print >> sys.stderr, "Sorry,", ext, "files seem not to exist"
            return self.spath_set
        rexp = re.compile(self.rexp_spec + ext + '$', re.I | re.U)
        if not journal:
            self.journal = self.mk_dummy_journal(rexp, files)
        else:
            self.journal = journal
        Sp = collections.namedtuple("Sp", "id flag pth fdpth plpth Vr sn jy jd")
        for path in files:
            night, flag = rexper(rexp, path, reflag=False, splitnig=splitnig)
            if flag and iflag and flag not in iflag:
                if self.verbose:
                    print "Exclude", night, "with flag", flag
                continue
            if night not in self.journal:
                continue
            jy, mjd, spc, lwv, hwv, Vr, sn, fdsnam = self.journal[night][:8]
            #jy, mjd, spc, lwv, hwv, Vr, sn, fdsnum = self.journal[night]
            pls_pth = os.path.splitext(path)[0]+'.ccm.txt'
            if not os.path.isfile(pls_pth):
                print >> sys.stderr, "Exclude", path+": no ccm.txt file"
                continue
            fds_pth = os.path.join(os.path.dirname(path), str(fdsnam)+'.fds')
            if not os.path.isfile(fds_pth):
                fds_pth = self.find_fds(path, str(night).zfill(6))
                if not fds_pth:
                    print >> sys.stderr, "Exclude", path+": no fds file"
                    continue
            if self.verbose:
                print "Append", night, flag, "with shift", Vr
            params = Sp(night, flag, path, fds_pth, pls_pth, Vr, sn, jy, mjd)
            self.spath_set.add(params)
        return self.spath_set

    def mk_dummy_journal(self, rexp, files):
        """Make dummy journal with paths to fds and wavelength limits
        """
        journal = {}
        for path in files:
            night, flag = rexper(rexp, path, reflag=False, splitnig=False)
            if night not in journal:
                fds_pth = self.find_fds(path, str(night).zfill(6))
                if not fds_pth:
                    print "Exclude", path, "no fds"
                    continue
                lwv, hwv = get_wv_range(fds_pth)
                fds_pth = os.path.basename(fds_pth)
                #jy, mjd, spc, lwv, hwv, Vr, sn, fdsnam
                journal[night] = [0, 0, '', lwv, hwv, 0, 1, fds_pth]
        return journal

    def find_fds(self, path, night):
        """Try to find fds, using path to file with spectrum.
        We assume that fds and spectrum name are similar and search fds
        given at the same night, that of spectrum.
        """
        fdses = {}
        nig, num = night[:3], int(night[3:])
        rxp = '*'+nig+'*'+'.fds'
        if os.path.dirname(path):
            rxp = os.path.dirname(path) + os.sep + rxp
        fds_pths = glob.glob(rxp)
        if not fds_pths:
            return
        for fdpth in fds_pths:
            fdnam = os.path.basename(fdpth)
            night, img, flag = self.fdsrexp.findall(fdnam)[0]
            if img:
                fdses[abs(num-int(img))] = fdpth
        return fdses[min(fdses)]

    def select_images(self):
        """Select image files with filenames matching rexp_fits or rexp_mt.
        Currently only fits, fts and mt file extensions are supported.
        """
        rexpfts = re.compile(self.rexp_fits, re.I | re.U)
        rexpmt = re.compile(self.rexp_mt, re.I | re.U)
        Sp = collections.namedtuple("Sp", "pth num date flag keymode")
        for ext in self.files_dct:
            for pth in self.files_dct[ext]:
                fnam = os.path.basename(pth)
                if ext.lower() in ("fits", "fts") \
                                        and len(rexpfts.findall(fnam)):
                    year, mon, day, num, flag = rexpfts.findall(fnam)[0]
                elif ext.lower() == "mt" and len(rexpmt.findall(fnam)):
                    year, mon, day, num, flag = rexpmt.findall(fnam)[0]
                    mon = str(int(mon, 16)).zfill(2)
                    year = str(2000 + int(year)) if int(year) else "2010"
                else:
                    continue
                keymode = False
                # Set flag for incorrect FITS header cards in NES spectra
                # taken in 2008:
                if 20080422 < int("".join((year, mon, day))) < 20081103:
                    keymode = "flstr"
                date = tuple(map(int, (year, mon, day)))
                params = Sp(pth, num, date, flag, keymode)
                self.path_set.add(params)
        return self.path_set


def read_journal(j_pth): return cPickle.load(open(j_pth, 'rb'))


def mod_journal(journal, sn=10, rvc=True, rvcorr_dct=None, instr=None,
                lams=None, jd=None, verbose=False):
    """Exclude or correct items in journal depending on input parameters.
    @param journal: dictionary with spectra main parameters
    @param sn: limiting signal-to-noise ratio value: exclude all spectra with
        lower S/N
    @param rvc: Set Va shift to zero. Useful for telluric line measurements
    @param rvcorr_dct: Dictionary of heliocentric velocity corrections,
        defined, for example, by telluric lines measurements
    @param instr: Select spectra obtained with given spectrograph
    @param lams: Exclude all spectra, which no contain given vawelength
    @param jd: Exclude all spectra, which obtained earlier, than given JYear
        value
    """
    if not rvc:
        print "Set Va to zero"
    for nig, dat in journal.items():
        jy, mjd, spc, lwv, hwv, Vr, jsn, fn = dat[:8]
        if sn and jsn < sn:
            if verbose:
                print "Exclude", nig, "S/N", sn
            del journal[nig]
            continue
        if instr and spc not in instr:
            if verbose:
                print "Exclude", nig, "instr", instr
            del journal[nig]
            continue
        if lams and (lwv > lams[1] or hwv < lams[0]):
            if verbose:
                print "Exclude", nig, "wv limits", lwv, hwv
            del journal[nig]
            continue
        if jd and jy < jd:
            if verbose:
                print "Exclude", nig, "JD", jd
            del journal[nig]
            continue
        if not rvc:
            if type(journal[nig]) is not list:
                journal[nig] = list(journal[nig])
            journal[nig][5] = 0
        if rvcorr_dct and nig in rvcorr_dct:
            if type(journal[nig]) is not list:
                journal[nig] = list(journal[nig])
            journal[nig][5] += rvcorr_dct[nig]


def clocker(val, divc=1.):
    """Convert float value to tuple of (int h, int m, float s).
    Useful for coordinates and times.
    @param divc: division coefficient for val
    """
    h, rest = divmod(val/divc, 1)
    m, rest = divmod(rest * 60, 1)
    s = round(rest * 60, 2)
    return int(h), int(m), s


class Logger:
    """Extract information from FITS headers of given FITS files, create log
    @param fitsfixmode: when FIX header is incorrect, try to fix it with given
        option. Can be 'fix' and 'silentfix'
    """
    def __init__(self, fitsfixmode='silentfix', verbose=False):
        self.fitsfixmode = fitsfixmode
        self.verbose = verbose
        # Define keys, presented in FITS headers:
        self.stdkeys = ['AUTHOR', 'REFERENC', 'DATE', 'ORIGIN', 'DATE-OBS',
            'TELESCOP', 'INSTRUME', 'OBSERVER', 'OBJECT', 'EQUINOX', 'EPOCH']
        self.nonstdkeys = ['OBSERVAT', 'CHIPID', 'DETNAME', 'CCDTEMP',
                           'P_DEWAR', 'SEEING', 'PROG-ID']
        #self.otherkeys = ['WIND', 'T_OUT', 'T_IN', 'T_MIRR', 'PRES',
                           #'ST', 'MT', 'FOCUS']
        self.journal_dct = {}
        self.obs_dct = {}
        self.key_dct = {"bias": 'b', "thar": 't', "flatfield": 'f', "ff": 'f',
                        "ll": 'f', "sky": 's', "ref.ima.": 'r'}
        self.postfix, self.ext = '.', 'fits'
        self.print_str = "[01;31m%s not found in header! Set it to zero.[0m"
        self.keyset = lambda key: self.key_dct.get(key.lower()) \
                if self.key_dct.get(key.lower()) else 'o'
        self.trygetkey = lambda head, *args: filter(head.has_key, args)
        self.seper = lambda fstr, h, m, s: \
                "%s" % ":".join(["%02d" % h, "%02d" % m, fstr % s])
        self.Sp = collections.namedtuple("Sp",
                "date dateobs flag obj time exp ra dec az zen pth")
        self.biases, self.flats, self.calib = {}, {}, {}

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
            curfts.verify(self.fitsfixmode)
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
        exp = int(self.pickhead(head, 'raw', 'EXPTIME'))
        obj = self.pickhead(head, 'raw', 'OBJECT')
        flag = self.keyset(obj)
        if flag == "b" and exp != 0:
            flag = "o"
        elif exp == 0 and flag != "b":
            flag, obj = "b", "BIAS"

        if flag == "b":
            newname = "".join(('b', num, self.postfix, self.ext))
            self.biases[file_pth] = newname
        elif flag == "f":
            newname = "".join(('f', num, self.postfix, self.ext))
            self.flats[file_pth] = newname
        elif flag == "t":
            newname = "".join(('l', num, "t", self.postfix, self.ext))
            self.calib[file_pth] = newname

        dateobs = head.get('DATE-OBS')
        if self.verbose:
            self.warner(exp, obj, num)
        self.obs_dct[int(num)] = self.Sp(date, dateobs, flag, obj,
                time, exp, ra, dec, az, zen, file_pth)

    def warner(self, exp, obj, num):
        """Warn about some problems."""
        if exp == 0 and str(obj).lower() != "bias":
            print "Wrong object name:", obj, "with exp=0 in spectrum", num
        #if head.get('CCDTEMP') != None and -103 > head.get('CCDTEMP') > -97:
            #print "incorrect CCDTEMP", ccdtemp, "in spectrum number", num

    def select_nums(self, nums, prfx="f"):
        """Select numbers in given range from obs_dct
        @param nums: two boundary spectrum numbers as tuple or list
        @param prfx: prefix in file name
        """
        out = {}
        for i in xrange(nums[0], nums[1]+1):
            if i in self.obs_dct:
                num = str(i).zfill(3)
                newname = "".join((prfx, num, self.postfix, self.ext))
                out[self.obs_dct[i].pth] = newname
        return out

    def logtostr(self):
        log = []
        for num, val in sorted(self.obs_dct.items()):
            log.append(self.logstringer(num, val))
        return log

    def journaltostr(self):
        jrn = []
        for key, val in sorted(self.journal_dct.items()):
            if not val or val == [" "]:
                continue
            elif key in ("CCDTEMP", "P_DEWAR") and val:
                val = [str(sum(val)/len(val))]
            elif key in ("DATE-OBS", "ORIGIN", "AUTHOR", "OBSERVER",
                         "SEEING") and val:
                val = [", ".join(map(str, val))]
            jrn.append("%8s %s" % (key, " ".join(map(str, val))))
        return jrn

    def logstringer(self, num, val):
        """Make a string for given parameters of observation log."""
        #date, flag, obj, time, exp, ra, dec, az, zen, dateobs = val
        date = "".join((str(val.date[0]), str(val.date[1]).zfill(2),
                        str(val.date[2]).zfill(2))) +" "+ val.flag
        #date = "".join(map(str, date))+" "+flag
        time = self.seper("%02d", *val.time)
        ra = self.seper("%04.1f", *val.ra)
        dec = self.seper("%04.1f", *val.dec)
        az = self.seper("%04.1f", *val.az)
        zen = self.seper("%04.1f", *val.zen)
        return "%3d %10s %15s %8s %4d %11s %11s %11s %10s %s" % (int(num),
                date, val.obj, time, val.exp, ra, dec, az, zen, val.dateobs)
