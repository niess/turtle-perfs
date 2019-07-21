local ffi = require "ffi"

local M = {}


-- Load the C library
local include = function(tag)
    local f = io.open("include/" .. tag .. ".h")
    ffi.cdef(f:read("*a"))
    f:close()

    ffi.load("lib/lib" .. tag .. ".so", true)
end
include "view"


-- Iterator over the lines of sight of a view
local function iterator (view)
    -- Initialise the iterator
    local naz, nel, az0, daz, del
    if ffi.string(view.site) == "CDC" then
        naz, nel = 601, 301
        az0 = -6.
        daz, del = 0.1, 0.1
    else
        naz, nel = 1801, 241
        az0 = -6.
        daz, del = 0.2, 0.05
    end

    return coroutine.wrap(function ()
        -- Run the loop
        for i = 0,naz-1,1
        do
            for j = 0,nel-1,1
            do
                -- Yield the current loop data
                coroutine.yield {
                    i = i,
                    j = j,
                    azimuth = az0 + daz * i,
                    elevation = del * j
                }
            end
        end
    end)
end


-- Wrap the view object
M.View = ffi.metatype("struct view", {
    __new = function(ct, site)
        local o = ffi.new(ct)
        ffi.C.view_initialise(o, site)
        return o
    end,
    __index = {
        lines_of_sight = iterator
    }
})


return M
