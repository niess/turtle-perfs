local ffi = require "ffi"

local M = {}


-- Load the C library
local include = function(tag)
    local f = io.open("include/" .. tag .. ".h")
    local header = f:read("*a")
    header = header:gsub("#[%w%p ]*", "")
    ffi.cdef(header)
    f:close()

    ffi.load("lib/lib" .. tag .. ".so", true)
end
include "simulator"


-- Wrap the simulator
local Simulator = ffi.metatype("struct simulator", {
    __new = function(ct)
        local o = ffi.new(ct)
        ffi.C.simulator_initialise(o)
        return o
    end
})
local s = Simulator()


function M.record(cb) s.record = cb end


function M.simulate(geometry, view, azimuth, elevation)
    local depth = ffi.new "double[1]"
    local flux = ffi.new "double[1]"
    local time = ffi.new "double[1]"
    local direction = ffi.new "double[3]"
    view:direction(azimuth, elevation, direction)

    local n = s:simulate(geometry, view.position, direction, depth, flux, time)
    return n, depth[0], flux[0], time[0]
end


function M.configure(opts)
    local function args()
        local events = opts.events ~= nil and opts.events or -1
        local precompute = opts.precompute ~= nil and opts.precompute or -1
        return events, precompute
    end
    ffi.C.simulator_configure(s, args())
end

return M
