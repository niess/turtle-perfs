#!/usr/bin/env luajit
require "strict"
local argparse = require "argparse"

local geant4 = require "geant4"
local view = require "view"


-- Parse the command line arguments
local parser = argparse()
    :name "geant4"
    :description "Run a Geant4 simulation along a line of sight."

parser:argument("site", "Site name (CDC or ULS).")
parser:argument("azimuth", "Azimuth angle (deg).")
    :convert(tonumber)
parser:argument("elevation", "Elevation angle (deg).")
    :convert(tonumber)
parser:option("-p --period", "Down-sampling period.")
    :default(100)
    :convert(tonumber)
parser:option("-e --events", "Number of Monte Carlo events.")
    :default(1)
    :convert(tonumber)
parser:option("--particle", "Primary particle.")
    :default("Geantino")
parser:option("--energy", "Kinetic energy (GeV).")
    :default(1E+04)
    :convert(tonumber)
parser:flag("-t --tessellate", "Use the native tessellation.")
parser:flag("--stl", "Tessellate and dump an STL file.")
parser:flag("--secondaries", "Track secondaries.")
parser:option("-c --cut", "Production cut (m).")
    :default(1.)
    :convert(tonumber)
parser:flag("-v --verbose", "Extra verbosity.")

local args = parser:parse()


-- Instantiate the Geant4 engine
local topography = function()
    -- The topography data and their down-sampling period
    local map_name = { CDC = "pdd.png", ULS = "tianshan.png" }
    return {area = map_name[args.site], period = args.period}
end

local mode
if args.stl then
    mode = 2
elseif args.tessellate then
    mode = 1
else
    mode = 0
end

geant4.initialise{
    area = ({ CDC = "pdd", ULS = "tianshan" })[args.site],
    period = args.period,
    mode = mode,
    secondaries = args.secondaries,
    cut = args.cut
}


-- Do the simulation
local site = view.View(args.site)
geant4.configure{events = args.events, verbosity = args.verbose and 2 or 0}

print("--- Geant4")
local n, d, t = geant4.run(site, args.particle, args.energy, args.azimuth,
                           args.elevation)
print(" steps  depth (m)     time (s)")
print(string.format("%4d    %9.6f  %10.7f", n, d, t))

