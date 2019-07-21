#!/usr/bin/env luajit
require "strict"
local argparse = require "argparse"

local geometry = require "geometry"
local simulator = require "simulator"
local view = require "view"


-- Parse settings from the command line
local parser = argparse()
    :name "scan-simulate"
    :description "Run a Monte Carlo scan"

parser:argument("site", "Site name.")
parser:argument("algorithm", "Raytracing algorithm.")
parser:argument("period", "Downsampling period.")
    :convert(tonumber)
parser:flag("--embree", "Use Embree instead of CGAL.")
parser:flag("--precompute", "Precompute the geometry.")
parser:option("-n --events", "Number of Monte Carlo events")
    :default(1)
    :convert(tonumber)
parser:option("--resolution", "Resolution parameter (optimistic)")
    :default(1E-02)
    :convert(tonumber)
parser:option("--range", "Range parameter (optimistic)")
    :default(1E+00)
    :convert(tonumber)
parser:option("--slope", "Slope parameter (optimistic)")
    :default(0.4)
    :convert(tonumber)

local args = parser:parse()
if args.embree then args.algorithm = "embree" end


-- Instanciate the geometry
local topography = function()
    -- The topography data and their downsampling period
    local map_name = { CDC = "pdd.png", ULS = "tianshan.png" }
    return "share/topography/" .. map_name[args.site], args.period
end

local geo = geometry.Geometry(args.algorithm, topography())
if args.algorithm == "opti" then
    geo:configure{resolution = args.resolution, range = args.range,
                  slope = args.slope}
elseif (args.algorithm == "BVH") or (args.algorithm == "embree") then
    geo:configure{mode = 0}
end


-- Run the scan
local site = view.View(args.site)
simulator.configure{events = args.events, precompute = args.precompute}

local cpu, steps, depth, flux = 0., 0, 0., 0.
for los in site:lines_of_sight()
do
    local n, d, f, t = simulator.simulate(geo, site, los.azimuth, los.elevation)
    cpu, steps, depth, flux = cpu + t, steps + n, depth + d, flux + f
    print(string.format("%4d %4d  %.5E %.5E  %d  %.5E",
        los.i, los.j, d, t, n, f))
    io.flush()
end


-- Print summary statistics
print(string.format("# mean depth    : %.5E m", depth / site.LOS))
print(string.format("# mean flux     : %.5E m", flux / site.LOS))
print(string.format("# total steps   : %d", steps))
print(string.format("# steps per LOS : %.3f", steps / site.LOS))
print(string.format("# total cpu     : %.5E s", cpu))
print(string.format("# cpu per LOS   : %.3f ms", cpu / site.LOS * 1E+03))
print(string.format("# cpu per step  : %.3f mus", cpu / steps * 1E+06))
