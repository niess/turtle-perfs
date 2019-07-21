local ffi = require "ffi"

local M = {}


-- Load the C library
local include = function(tag)
    local f = io.open("include/" .. tag .. ".h")
    ffi.cdef(f:read("*a"))
    f:close()

    ffi.load("lib/lib" .. tag .. ".so", true)
end
include "geant4"


-- Wrap the geant4 functions
local Geant4 = ffi.metatype("struct geant4", {
    __new = function(ct, ...)
        local o = ffi.new(ct)
        ffi.C.geant4_initialise(o, ...)
        return o
    end,
    __gc = function(o)
        o.clear()
    end
})
local g4 = {}


function M.initialise(opts)
    local function args()
        local period = opts.period ~= nil and opts.period or 1
        local secondaries = not opts.secondaries
        local cut = opts.cut ~= nil and opts.cut or 1.
        local mode = opts.mode ~= nil and opts.mode or
                     ffi.C.GEANT4_MODE_G4TURTLE
        return opts.area, period, mode, secondaries, cut
    end
    g4 = Geant4(args())
end


function M.run(view, particle, energy, azimuth, elevation)
    local pid = ({
        Geantino = ffi.C.GEANT4_PARTICLE_GEANTINO,
        MuonMinus = ffi.C.GEANT4_PARTICLE_MUON_MINUS,
        MuonPlus = ffi.C.GEANT4_PARTICLE_MUON_PLUS
    })[particle]

    local depth = ffi.new "double[1]"
    local time = ffi.new "double[1]"
    local direction = ffi.new "double[3]"
    view:direction(azimuth, elevation, direction)

    local n = g4.run(pid, energy, view.position, direction, depth, time)
    return n, depth[0], time[0]
end


function M.configure(opts)
    local function args()
        local events = opts.events ~= nil and opts.events or -1
        local verbosity = opts.verbosity ~= nil and opts.verbosity or -1
        return events, verbosity
    end
    ffi.C.geant4_configure(g4, args())
end

return M
