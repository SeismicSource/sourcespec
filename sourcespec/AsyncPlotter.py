# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
AsyncPlotter class for sourcespec.

:copyright:
    2016-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)

Modified from: https://gist.github.com/astrofrog/1453933
"""
import multiprocessing as mp
import time
from sourcespec.savefig import savefig


def _async_plotter(nc, processes, fig, filename, fmt, **kwargs):
    while nc.value >= processes:
        time.sleep(0.1)
    nc.value += 1
    savefig(fig, filename, fmt, **kwargs)
    nc.value -= 1


class AsyncPlotter():
    def __init__(self, processes=mp.cpu_count()):
        self.manager = mp.Manager()
        self.nc = self.manager.Value('i', 0)
        self.pids = []
        self.processes = processes

    def save(self, fig, filename, fmt, **kwargs):
        p = mp.Process(
            target=_async_plotter,
            args=(self.nc, self.processes, fig, fmt, filename),
            kwargs=kwargs
        )
        p.start()
        self.pids.append(p)

    def terminate(self):
        for p in self.pids:
            p.terminate()
            p.join()

    def join(self):
        for p in self.pids:
            p.join()
