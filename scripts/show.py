#!/usr/bin/env python3

import sys

import numpy
import matplotlib.pyplot as plot
plot.style.use("share/mplstyle/l3.mplstyle")


def show(path):
    data = numpy.loadtxt(path, comments="#")
    if len(data) == 601 * 301:
        shape = (601, 301)
        azimuth = numpy.linspace(-6., 54., 601)
        elevation = numpy.linspace(0., 30., 301)
        vmax = 4.
        xticks = (0, 10, 20, 30, 40, 50)
    else:
        shape = (1801, 241)
        azimuth = numpy.linspace(0., 360, 1801)
        elevation = numpy.linspace(0., 12., 241)
        vmax = 100.
        xticks = (0, 90, 180, 270, 360)
    depth = data[:,2].reshape(shape).T

    cmap = plot.get_cmap("jet")
    cmap.set_under("white")

    plot.figure()
    depth *= 1E-03
    depth[depth <= 0] = -1.
    plot.pcolor(azimuth, elevation, depth, cmap=cmap,
                vmin=0, vmax=vmax)
    plot.colorbar()
    plot.xlabel("azimuth (deg)")
    plot.ylabel("elevation (deg)")
    plot.xticks(xticks)
    plot.show()


if __name__ == "__main__":
    show(sys.argv[1])
