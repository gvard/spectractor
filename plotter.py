'''This module is part of Spectractor. It contains functions for
plotting.
'''


import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np

from spectractor import Spliner


class Flatter(Spliner):
    """Class for interactive flatting set of 1D data.
    @param splt: matplotlib subplot object
    @param data: iterable array of 1D data
    @param ccm: optional array of supplement reference points
    @param blaze: optional array of supplement curves
    @param zw: half width of zooming window
    @param spk: degree of spline
    """
    def __init__(self, splt, data, ccm=None, blaze=None,
                clr="c", lw=1, grid=False, zw=100, spk=3):
        self.splt = splt
        self.data = data
        self.ccm = ccm
        self.blz = blaze
        # Plot params
        self.clr = clr
        self.lw = lw
        self.grid = grid
        self.zw = zw
        self.kpe, self.bpe = "key_press_event", "button_press_event"
        self.tb = plt.get_current_fig_manager().toolbar
        self.nord, self.end = len(self.data), len(self.data[0])
        self.xlims = (0, self.end)
        # Initialization
        self.Sl = Spliner(0, self.end, smooth=None, spk=spk)
        self.out = np.ones((self.nord, self.end))
        self.init_prm()
        self.i = 0
        self.dat = self.data[0]
        self.plotdat()

    def init_prm(self):
        self.plcrv, self.plspl, self.cmpl = [], [], []
        self.coords = []
        self.isdiv, self.zoomed = False, False

    def plotdat(self):
        _tmp, = self.splt.plot(self.dat, '-k', lw=self.lw)
        self.plcrv.append(_tmp)
        if self.ccm is not None and self.i+1 in self.ccm:
            for x, y in self.ccm[self.i+1]:
                self.cmpl.append(self.splt.plot(x, y, 'o', c='g', ms=5)[0])
            self.splfit = self.mksplfit(self.ccm[self.i+1])
            self.plspl.append(self.splt.plot(self.splfit, '-g')[0])
        elif self.blz is not None:
            _tmp, = self.splt.plot(self.blz[self.i], '-g', lw=self.lw)
            self.plcrv.append(_tmp)
        self.set_lims(self.dat)
        self.splt.grid(self.grid)
        plt.draw()

    def plotspl(self, delprev=True):
        self.splfit = self.mksplfit(self.coords)
        if delprev:
            [plelm.remove() for plelm in self.plspl]
            self.plspl = []
        self.plspl.append(self.splt.plot(self.splfit, '-r')[0])
        plt.draw()

    def mksplfit(self, points):
        xfl, yfl = zip(*sorted(points))
        return self.Sl.splpl(xfl, yfl)

    def plotdiv(self):
        self.clear()
        self.dat = self.dat / self.splfit
        self.isdiv = True
        self.plcrv.append(self.splt.plot(self.dat, '-k', lw=self.lw)[0])
        for x, y in self.coords:
            self.cmpl.append(self.splt.plot(x, 1, 'o', c=self.clr, ms=5)[0])
        self.set_lims(self.dat)
        plt.draw()

    def clear(self):
        [plelm.remove() for plelm in self.plcrv]
        [plelm.remove() for plelm in self.plspl]
        [plelm.remove() for plelm in self.cmpl]
        self.init_prm()

    def delplt(self, ind):
        self.cmpl[ind].remove()
        del self.cmpl[ind]
        del self.coords[ind]
        plt.draw()

    def set_lims(self, data, ysh=None):
        self.splt.set_xlim(self.xlims)
        dmin, dmax = np.min(data), np.max(data)
        if not ysh:
            ysh = (dmax - dmin) / 50.
        self.splt.set_ylim(dmin-ysh, dmax+ysh)

    def getcoords(self, event, rndl=3):
        return round(event.xdata, rndl), round(event.ydata, rndl)

    def addpnt(self, event):
        x, y = self.getcoords(event)
        if event.button == 3:
            y = self.dat[int(round(x))]
        self.coords.append((x, y))
        self.cmpl.append(self.splt.plot(x, y, 'o', c=self.clr, ms=5)[0])
        plt.draw()

    def next(self):
        if self.isdiv:
            self.out[self.i] = self.dat
        self.clear()
        if self.i < self.nord - 1:
            self.i += 1
            self.dat = self.data[self.i]
            self.plotdat()
        else:
            plt.close()

    def zoom(self, event):
        if not self.zoomed:
            x, y = self.getcoords(event)
            x0, x1 = max(x-self.zw, 0), min(self.end, x+self.zw)
            self.zoomed = True
        else:
            x0, x1 = 0, self.end
            self.zoomed = False
        self.xlims = (x0, x1)
        self.set_lims(self.dat[x0:x1])
        plt.draw()

    def __call__(self, evnt):
        """@param evnt: event"""
        if not self.tb.mode and evnt.inaxes and evnt.name == self.bpe \
                and evnt.button in (1, 3):
            self.addpnt(evnt)
        elif evnt.name == self.kpe and evnt.inaxes and evnt.key == ' ':
            self.zoom(evnt)
        elif evnt.name == self.kpe and evnt.key == 'f10' \
                and len(self.coords) > 3:
            self.plotspl()
        elif evnt.name == self.kpe and evnt.key == 'f9':
            self.plotdiv()
        elif evnt.name == self.kpe and evnt.key == 'd' and len(self.coords):
            self.delplt(-1)
        elif evnt.name == self.kpe and evnt.key == 'D' and len(self.coords):
            self.delplt(0)
        elif evnt.name == self.kpe and evnt.key == 'n':
            self.next()
        elif evnt.name == self.kpe and evnt.key == 'q':
            plt.close()


def plt_opt(figsub, **kwargs):
    """Set default parameters for kosher plotting"""
    lblsz = 'x-small'
    if 'lblsz' in kwargs:
        lblsz = kwargs['lblsz']
    if 'xlbl' in kwargs:
        figsub.set_xlabel(kwargs['xlbl'], fontsize=lblsz)
    if 'ylbl' in kwargs:
        figsub.set_ylabel(kwargs['ylbl'], fontsize=lblsz)
    if 'grid' in kwargs:
        figsub.grid(True)
    if 'xlim' in kwargs:
        figsub.set_xlim(kwargs['xlim'])
    if 'ylim' in kwargs:
        figsub.set_ylim(kwargs['ylim'])
    if 'xlc' in kwargs:
        figsub.xaxis.set_major_locator(MultipleLocator(kwargs['xlc'][0]))
        figsub.xaxis.set_minor_locator(MultipleLocator(kwargs['xlc'][1]))
    if 'ylc' in kwargs:
        figsub.yaxis.set_major_locator(MultipleLocator(kwargs['ylc'][0]))
        figsub.yaxis.set_minor_locator(MultipleLocator(kwargs['ylc'][1]))
    if 'xmjf' in kwargs:
        figsub.xaxis.set_major_formatter(FormatStrFormatter(kwargs['xmjf']))
    if 'ticksz' in kwargs:
        #figsub.set_xticklabels(fontsize=kwargs['ticksz'])
        figsub.ticklabel_format(fontsize=kwargs['ticksz'])
    if 'legend' in kwargs:
        leglblsz = 'x-small'
        if 'leglblsz' in kwargs:
            leglblsz = kwargs['leglblsz']
        leg = figsub.legend(numpoints=1, loc=kwargs['loc'], #'best'
                    labelspacing=kwargs['legend'])
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=leglblsz)
