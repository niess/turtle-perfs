#!/usr/bin/env luajit
require "strict"
local argparse = require "argparse"

local geometry = require "geometry"
local tracer = require "tracer"
local view = require "view"


-- Parse settings from the command line
local parser = argparse()
    :name "scan-trace"
    :description "Perform a raytracing scan"

parser:argument("site", "Site name (CDC or ULS).")
parser:argument("algorithm", "Raytracing algorithm (BVH, mesh or opti).")
parser:argument("period", "Downsampling period.")
    :convert(tonumber)
parser:flag("--embree", "Use Embree instead of CGAL.")
parser:option("--resolution", "Resolution parameter (optimistic)")
    :default(1E-02)
    :convert(tonumber)
parser:option("--range", "Range parameter (optimistic)")
    :default(1E+00)
    :convert(tonumber)
parser:option("--slope", "Slope parameter (optimistic)")
    :default(0.4)
    :convert(tonumber)
parser:flag("--all", "All mode (BVH)")
parser:flag("--first", "First mode (BVH/CGAL)")

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
elseif (args.algorithm == "BVH") then
    local mode
    if args.first then
        mode = 2
    elseif args.all then
        mode = 0
    else
        mode  = 1
    end
    geo:configure{mode = mode}
elseif (args.algorithm == "embree") then
    local mode
    if args.all then
        mode = 0
    else
        mode  = 1
    end
    geo:configure{mode = mode}
end


-- Run the scan
local site = view.View(args.site)

local cpu, steps, depth = 0., 0, 0.
for los in site:lines_of_sight()
do
    local n, d, t = tracer.trace(geo, site, los.azimuth, los.elevation)
    cpu, steps, depth = cpu + t, steps + n, depth + d
    print(string.format("%4d %4d  %.5E %.5E  %d", los.i, los.j, d, t, n))
    io.flush()
end


-- Print summary statistics
print(string.format("# mean depth    : %.5E m", depth / site.LOS))
print(string.format("# total steps   : %d", steps))
print(string.format("# steps per LOS : %.3f", steps / site.LOS))
print(string.format("# total cpu     : %.5E s", cpu))
print(string.format("# cpu per LOS   : %.3f ms", cpu / site.LOS * 1E+03))
print(string.format("# cpu per step  : %.3f mus", cpu / steps * 1E+06))
