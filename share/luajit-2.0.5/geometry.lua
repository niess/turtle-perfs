-- Geometry package
local ffi = require "ffi"

local M = {}


-- Load the geometry header
local include = function(tag)
    local f = io.open("include/" .. tag .. ".h")
    local header = f:read("*a")
    header = header:gsub("#[%w%p ]*", "")
    ffi.cdef(header)
    f:close()
end
include "geometry"


local libs = {}
local libname = {
    BVH = "bvh",
    embree = "embree",
    mesh = "mesh",
    opti = "opti"
}


-- Convert named arguments (options) to positional ones
local function Args(kwargs)
    return function(opts)
        local args = {}
        for _, index in pairs(kwargs)
        do
            args[index] = -1
        end
        for k, v in pairs(opts)
        do
            local index = kwargs[k]
            if index == nil then
                error("invalid option " .. k, 1)
            end
            args[index] = v
        end
        return args
    end
end

local libopts = {
    bvh = Args{mode = 1},
    embree = Args{mode = 1},
    mesh = Args{},
    opti = Args{resolution = 1, range = 2, slope = 3}
}


-- Wrap the geometry object
M.Geometry = ffi.metatype("struct geometry", {
    __new = function(ct, algorithm, path, ...)
        local o = ffi.new(ct)
        print("# Creating " ..  algorithm .. " geometry from " .. path)
        local name = libname[algorithm]
        if not libs[algorithm] then
            libs[algorithm] = ffi.load("lib/lib" .. name .. ".so", true)
        end
        name = "geometry_initialise_" .. name
        ffi.C[name](o, path, ...)
        return o
    end,
    __gc = function(o)
        if o.clear ~= nil then o.clear() end
    end,
    __index = {
        configure = function(o, opts)
            local name = libname[ffi.string(o.algorithm)]
            local args = libopts[name](opts)
            name = "geometry_configure_" .. name
            ffi.C[name](o, unpack(args))
        end
    }
})


return M
