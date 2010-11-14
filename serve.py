#!/usr/bin/env python

'''This module is part of Spectractor. It contains functions for
managing various files and journals.

@author: Dmitry Nasonov
'''

import re
import os, sys, shutil
import collections

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
    @param regexp: regular expression describing file name.
    @param path: path to file.
    @param reflag: tuple of two values. If not empty, flag containing in first
            value will be replaced to second one.
    @param splitnig: split ID to tuple of 'night' and spectrum number.
    @return ID (night number): int or tuple of ints with splitnig, flag.
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
    @param namrxp: filename regular expression, but without extension.
    @param args: list of paths to files and dirs.
    """
    def __init__(self, args, exts=("200", "fds"), namrxp=".*", verbose=False):
        self.files_dct = {}
        self.workdirs = []
        self.path_set, self.spath_set = set(), set()
        self.verbose = verbose
        # rexp_spec contains 3 groups: night number, spectrum number and flag.
        self.rexp_spec = "^[a-z_]{0,3}(\d{3})(\d{1,3})[-_]?([a-z_\-]*).*\."
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
        """Walk on workdirs, return list of all subdirs.
        @return list of subdirectories.
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
                we want exclude some subdir.
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
        @param ext: select file type, must be file extension.
        @param filtr: selecting only given flags, exclude others. Select all
                when filtr in ("all", None).
        @param reflag: tuple of two values. If not empty, flag containing in
                first value will be replaced to second one.
        @param splitnig: split ID to tuple of 'night' and spectrum number.
        @return: list of tuples with 'night', flag and path to file.
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

    def select_spectra(self, journal, ext='200', iflag=None, splitnig=False):
        """Select data, get ID ('night') and flag from filename, collect paths.
        Here we assume, that file name contains information about ID
        (night number, number of spectrum) and flag. Others will be broken
        (or you must change default regular expression). Then we check presence
        of ID in journal keys, get spectrum parameters from journal and append
        this data to resulting list.
        @param journal: journal dictionary with spectra IDs as keys and spectra
        parameters (for example, JYear, MJD, Spectrograph key, limiting
                wavelengths, S/N etc) as values.
        @param ext: select file type, must be file extension.
        @param iflag: Flags to select. Can be a string, if flag is one symbol.
        @param splitnig: split 'night' to tuple of night and spectrum number.
        @return: list of tuples with ID ('night'), flag, path to file and
                other parameters. Currently parameters are fds path, dictionary
                of points, defining the flatten polynome, heliocentric
                correction and S/N value.
        """
        Sp = collections.namedtuple("Sp", "id flag pth fdpth plpth Vr sn")
        try:
            files = self.files_dct[ext]
        except KeyError:
            print >> sys.stderr, "Sorry,", ext, "files seem not to exist."
            return self.spath_set
        rexp = re.compile(self.rexp_spec + ext + '$', re.I | re.U)
        for path in files:
            night, flag = rexper(rexp, path, reflag=False, splitnig=splitnig)
            if (flag and iflag) and flag not in iflag:
                if self.verbose:
                    print "Exclude", night, "with flag", flag
                continue
            if night not in journal:
                continue
            jy, mjd, spc, lwv, hwv, Vr, sn, fdsnum = journal[night]
            pls_pth = os.path.splitext(path)[0]+'.ccm.txt'
            if not os.path.isfile(pls_pth):
                if self.verbose:
                    print "Exclude", path, "with no ccm.txt"
                continue
            #night = int(str(spe[0][0]) + str(spe[0][1]).zfill(3))
            fds_pth = os.path.join(os.path.dirname(path),
                        str(night)[:3]+str(fdsnum).zfill(3)+'.fds')
            if self.verbose:
                print "Append", night, "with corr", Vr
            params = Sp(night, flag, path, fds_pth, pls_pth, Vr, sn)
            self.spath_set.add(params)
        return self.spath_set

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
