# -*- coding: utf8 -*-
# Modified from: https://gist.github.com/astrofrog/1453933
import multiprocessing as mp
import time


def _async_plotter(nc, processes, fig, filename, **kwargs):
    while nc.value >= processes:
        time.sleep(0.1)
    nc.value += 1
    fig.savefig(filename, **kwargs)
    nc.value -= 1


class AsyncPlotter():
    def __init__(self, processes=mp.cpu_count()):
        self.manager = mp.Manager()
        self.nc = self.manager.Value('i', 0)
        self.pids = []
        self.processes = processes

    def save(self, fig, filename, **kwargs):
        p = mp.Process(
            target=_async_plotter,
            args=(self.nc, self.processes, fig, filename),
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
