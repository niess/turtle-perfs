#!/usr/bin/env python3

import sys

import numpy
import matplotlib.pyplot as plot
from matplotlib.colors import LogNorm
plot.style.use("share/mplstyle/l3.mplstyle")


def load(path):
    data = numpy.loadtxt(path, comments="#")
    if len(data) == 601 * 301:
        shape = (601, 301)
        azimuth = numpy.linspace(-6., 54., 601)
        elevation = numpy.linspace(0., 30., 301)
        xticks = (0, 10, 20, 30, 40, 50)
    else:
        shape = (1801, 241)
        azimuth = numpy.linspace(0., 360, 1801)
        elevation = numpy.linspace(0., 12., 241)
        xticks = (0, 90, 180, 270, 360)
    depth = data[:,2].reshape(shape).T
    return azimuth, elevation, depth, xticks


def show(path1, path2):
    _, _, depth1, _ = load(path1)
    azimuth, elevation, depth2, xticks = load(path2)

    delta = numpy.absolute(depth1 - depth2)
    K = (depth1 > 0) | (depth2 > 0)

    print("mean: {:.3E}".format(delta[K].mean()))
    print("median: {:.3E}".format(numpy.median(delta[K])))
    print("max: {:.3E}".format(numpy.median(delta[K].max())))

    cmap = plot.get_cmap("jet")
    cmap.set_under("white")

    plot.figure()
    vmin, vmax = 1E-06, 1E-02
    delta[K] = numpy.maximum(delta[K], vmin)
    plot.pcolor(azimuth, elevation, delta, cmap=cmap,
                vmin=vmin, vmax=vmax, norm=LogNorm())
    plot.colorbar()
    plot.xlabel("azimuth (deg)")
    plot.ylabel("elevation (deg)")
    plot.xticks(xticks)

    plot.show()


if __name__ == "__main__":
    show(sys.argv[1], sys.argv[2])
