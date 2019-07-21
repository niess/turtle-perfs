#!/usr/bin/env luajit
require "strict"
local argparse = require "argparse"

local geant4 = require "geant4"
local view = require "view"


-- Parse settings from the command line
local parser = argparse()
    :name "scan-geant4"
    :description "Perform a Geant4 scan"

parser:argument("site", "Site name.")
parser:argument("period", "Downsampling period.")
    :default(1)
    :convert(tonumber)
parser:option("-e --events", "Number of Monte Carlo events.")
    :default(1)
    :convert(tonumber)
parser:option("--particle", "Primary particle.")
    :default("Geantino")
parser:option("--energy", "Kinetic energy (GeV).")
    :default(1E+04)
    :convert(tonumber)
parser:flag("-t --tessellate", "Use native tessellation.")
parser:flag("--stl", "Tessellate and dump an STL file.")
parser:flag("--secondaries", "Track secondaries.")
parser:option("-c --cut", "Production cut (m).")
    :default(1.)
    :convert(tonumber)

local args = parser:parse()


-- Instanciate the geant4 engine
local topography = function()
    -- The topography data and their downsampling period
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


-- Run the scan
local site = view.View(args.site)
geant4.configure{events = args.events, verbosity = 0}

local cpu, steps, depth = 0., 0, 0.
for los in site:lines_of_sight()
do
    local n, d, t = geant4.run(site, args.particle, args.energy, los.azimuth,
                               los.elevation)
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
