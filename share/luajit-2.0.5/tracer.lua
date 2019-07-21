local ffi = require "ffi"

local M = {}


-- Load the C library
local include = function(tag)
    local f = io.open("include/" .. tag .. ".h")
    ffi.cdef(f:read("*a"))
    f:close()

    ffi.load("lib/lib" .. tag .. ".so", true)
end
include "tracer"


-- Wrap the ray tracer
local Tracer = ffi.metatype("struct tracer", {
    __new = function(ct)
        local o = ffi.new(ct)
        ffi.C.tracer_initialise(o)
        return o
    end
})
local t = Tracer()


function M.record(cb) t.record = cb end


function M.trace(geometry, view, azimuth, elevation)
    local depth = ffi.new "double[1]"
    local time = ffi.new "double[1]"
    local direction = ffi.new "double[3]"
    view:direction(azimuth, elevation, direction)

    local n = t:trace(geometry, view.position, direction, depth, time)
    return n, depth[0], time[0]
end


return M
