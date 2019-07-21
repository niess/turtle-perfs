#!/usr/bin/env luajit
require "strict"
local argparse = require "argparse"

local geometry = require "geometry"
local tracer = require "tracer"
local view = require "view"


-- Parse the command line arguments
local parser = argparse()
    :name "trace"
    :description "Perform a ray tracing along a line of sight."

parser:argument("site", "Site name (CDC or ULS).")
parser:argument("azimuth", "Azimuth angle (deg).")
    :convert(tonumber)
parser:argument("elevation", "Elevation angle (deg).")
    :convert(tonumber)
parser:flag("--embree", "Use Embree instead of CGAL.")
parser:option("-p --period", "Downsampling period.")
    :default(100)
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
parser:flag("--all", "All mode (BVH)")
parser:flag("--first", "First mode (BVH/CGAL)")

local args = parser:parse()


-- Instanciate the geometries
local topography = function()
    -- The topography data and their downsampling period
    local map_name = { CDC = "pdd.png", ULS = "tianshan.png" }
    return "share/topography/" .. map_name[args.site], args.period
end

local bvh, mode
if args.embree then
    bvh = geometry.Geometry("embree", topography())
    mode = (not args.all) and 1 or 0
else
    bvh = geometry.Geometry("BVH", topography())
    if args.first then
        mode = 2
    elseif args.all then
        mode = 0
    else
        mode = 1
    end
end
bvh:configure{mode = mode}

local mesh = geometry.Geometry("mesh", topography())

local opti = geometry.Geometry("opti", topography())
opti:configure{resolution = args.resolution, range = args.range,
               slope = args.slope}


-- Bind a recorder to the ray tracer
local function Recorder()
    local last_medium
    return {
        record = function (step, medium, distance, r)
            if medium ~= last_medium then
                print(string.format("%4d     %2d    %9.3f",
                                    step, medium, distance))
            end
            last_medium = medium
        end,
        reset = function() last_medium = -1 end
    }
end
recorder = Recorder()
tracer.record(recorder.record)


-- Do the ray tracing
local site = view.View(args.site)

for tag, geo in pairs({BVH = bvh, mesh = mesh, opti = opti})
do
    recorder.reset()
    print("=== " .. tag .. " ===")
    print(" step  medium  length (m)")
    n, d, t = tracer.trace(geo, site, args.azimuth, args.elevation)
    print(" steps  depth (m)     time (s)")
    print(string.format("%4d    %9.6f  %10.7f", n, d, t))
end
