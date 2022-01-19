# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Grid importance sampling using a k-d tree.

:copyright:
    2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import itertools
import numpy as np
from scipy.interpolate import griddata


class KDTCell():
    def __init__(self, extent, calc_pdf, min_cell_prob=0,
                 ndiv=None, maxdiv=None):
        self.extent = extent
        self.coords = np.mean(extent, axis=1)
        self.delta = np.diff(extent)
        self.calc_pdf = calc_pdf
        self.pdf = calc_pdf(self.coords)
        self.volume = np.prod(self.delta)
        self.prob = self.pdf * self.volume
        self.min_cell_prob = min_cell_prob
        if self.prob <= min_cell_prob:
            self.is_divisible = False
            self.prob_divisible = 0.
        else:
            self.is_divisible = True
            self.prob_divisible = self.prob
        self.ndiv = ndiv
        if self.ndiv is None:
            self.ndiv = np.zeros(self.extent.shape[0])
        self.maxdiv = maxdiv

    def divide(self, parts):
        if not self.is_divisible:
            return [self, ]
        # dim is the dimension of the parameter space
        dim = self.extent.shape[0]
        # define a set of increments for all the spatial directions
        # the number of increments is "parts"
        # The code "np.mgrid...[1]" produces a 2D array of the type:
        #  [[0 1 ... parts-1]
        #   ... (dim) ...
        #   [0 1 ... parts-1]]
        # increments = np.mgrid[0:dim, 0:parts][1] * self.delta/(parts-1)
        increments = np.mgrid[0:dim, 0:parts][1]
        if self.maxdiv is not None:
            # dimensions that will not be divided
            increments[self.ndiv >= self.maxdiv] *= 0
            # dimensions that will be divided
            self.ndiv[self.ndiv < self.maxdiv] += parts
        increments = increments * self.delta/(parts-1)
        # take the minimum coordinate and transform it to a column vector
        mincoord = self.extent[:, 0].reshape(-1, 1)
        # add minimum coordinates to the increments
        coords = mincoord + increments
        cells = []
        # loop over all the possible n-uplets of coordinates
        # we use set() to avoid repetitions for dimensions that
        # will not be dvided
        for c in set(itertools.product(*coords)):
            # c is a coordinate n-uplet. Let's transform it to column vector
            c = np.array(c).reshape(-1, 1)
            delta = self.delta/parts
            # FIXME: I forgot my own code!
            #  should this `2` be `parts`? or is it related to delta/2?
            inc = np.mgrid[0:dim, 0:2][1] * delta - delta/2
            extent = c + inc
            cells.append(
                KDTCell(
                    extent, self.calc_pdf, self.min_cell_prob,
                    self.ndiv, self.maxdiv))
        return cells


class KDTree():
    def __init__(self, extent, init_parts, calc_pdf, min_cell_prob=0.,
                 maxdiv=None):
        # extent defines the size of search hypervolume
        # reshape extent to (dim, 2), where dim is the
        # arbitrary dimension of the parameter space
        self.extent = np.array(extent).reshape(-1, 2)
        # create the first cell, with the same size
        # of the search hypervolume
        cell0 = KDTCell(self.extent, calc_pdf, maxdiv=maxdiv)
        self.init_prob = cell0.prob
        cell0.min_cell_prob = cell0.prob*min_cell_prob
        self.cells = cell0.divide(init_parts)
        self.ncells = len(self.cells)

    def divide(self):
        # find the cell with highest probability and
        # divide it in 2 parts along every dimension
        # self.cells.sort(key=lambda c: c.prob)
        self.cells.sort(key=lambda c: c.prob_divisible)
        cell0 = self.cells.pop()
        # print(self.init_prob, cell0.prob, cell0.prob/self.init_prob)
        self.cells += cell0.divide(2)
        self.ncells = len(self.cells)

    def get_pdf(self, deltas):
        deltas = np.array(deltas).reshape(-1, 1)
        extent = self.extent.copy()
        ranges = []
        extent_new = []
        for v in np.hstack((extent, deltas)):
            start, stop, step = v
            # add a small number to make sure end value is inclued
            stop += stop/1e5
            rng = np.arange(start, stop, step)
            ranges.append(rng)
            extent_new += [rng[0], rng[-1]]
        xi = np.meshgrid(*ranges, indexing='ij')
        coords = np.array([cell.coords for cell in self.cells])
        pdf_values = np.array([cell.pdf for cell in self.cells])
        pdf = griddata(coords, pdf_values, tuple(xi), method='linear')
        return pdf, extent_new
