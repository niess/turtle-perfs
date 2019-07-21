#!/usr/bin/env luajit
require "strict"
local argparse = require "argparse"

local geometry = require "geometry"
local simulator = require "simulator"
local view = require "view"


-- Parse the command line arguments
local parser = argparse()
    :name "simulate"
    :description "Run a Monte Carlo simulation along a line of sight."

parser:argument("site", "Site name (CDC or ULS).")
parser:argument("azimuth", "Azimuth angle (deg).")
    :convert(tonumber)
parser:argument("elevation", "Elevation angle (deg).")
    :convert(tonumber)
parser:flag("--embree", "Use Embree instead of CGAL.")
parser:flag("--precompute", "Pre-compute the geometry.")
parser:option("-p --period", "Downsampling period.")
    :default(100)
    :convert(tonumber)
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


-- Instanciate the geometries
local topography = function()
    -- The topography data and their downsampling period
    local map_name = { CDC = "pdd.png", ULS = "tianshan.png" }
    return "share/topography/" .. map_name[args.site], args.period
end

local bvh
if args.embree then
    bvh = geometry.Geometry("embree", topography())
else
    bvh = geometry.Geometry("BVH", topography())
end
bvh:configure{mode = 0}

local mesh = geometry.Geometry("mesh", topography())

local opti = geometry.Geometry("opti", topography())
opti:configure{resolution = args.resolution, range = args.range,
               slope = args.slope}


-- Bind a recorder to the simulator
local function Recorder()
    local last_medium
    return {
        record = function(event, step, flag, medium, energy, distance, position)
            if (step == 0) or (medium ~= last_medium) or
               ((flag > 0) and ((flag<= 2016) or (flag >= 8192))) then
                print(string.format("%3d   %4d  %5d   %2d     %.3E   %9.3f",
                                    event, step, flag, medium, energy, distance))
            end
            last_medium = medium
        end,
        reset = function() last_medium = -1 end
    }
end
recorder = Recorder()
simulator.record(recorder.record)


-- Do the simulation
local site = view.View(args.site)
simulator.configure{events = args.events, precompute = args.precompute}

for tag, geo in pairs({BVH = bvh, mesh = mesh, opti = opti})
do
    recorder.reset()
    print("=== " .. tag .. " ===")
    print(" event step  flag medium  energy (GeV) length (m)")
    n, d, f, t = simulator.simulate(geo, site, args.azimuth, args.elevation)
    print(" steps  depth (m)     flux (a.u.)    time (s)")
    print(string.format("%4d  %12.6f    %.5E  %10.7f", n, d, f, t))
end
