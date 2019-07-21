cmake_minimum_required (VERSION 3.0.0)
project (TURTLE_TEST
    LANGUAGES C CXX
)

include (ExternalProject)


if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING
                "Choose the type of build, options are: Debug Release
                 RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set (CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH
            "Install path prefix." FORCE)
endif ()

set (CMAKE_MACOSX_RPATH 1)
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

link_directories (${CMAKE_INSTALL_PREFIX}/lib)
include_directories (${CMAKE_INSTALL_PREFIX}/include)


set (LUAJIT_VERSION "2.0.5" CACHE STRING "Version of LuaJIT")


# The generic geometry interface
install (FILES src/geometry.h DESTINATION include)


# The BVH geometry
add_library (bvh SHARED src/bvh.cpp)
target_link_libraries (bvh CGAL downsampler tictac turtle)
add_dependencies (bvh CGAL GEOMETRY TURTLE)
target_compile_definitions (bvh PRIVATE
    "-DINSTALL_PREFIX=\"${CMAKE_INSTALL_PREFIX}\"")
install (TARGETS bvh DESTINATION lib)


# The Embree geometry
add_library (embree SHARED src/embree.c)
target_link_libraries (embree embree3 downsampler tbb tbbmalloc tictac turtle)
add_dependencies (embree EMBREE GEOMETRY TURTLE)
target_compile_definitions (embree PRIVATE
    "-DINSTALL_PREFIX=\"${CMAKE_INSTALL_PREFIX}\"")
install (TARGETS embree DESTINATION lib)


# The mesh geometry
add_library (mesh SHARED src/mesh.cpp)
target_link_libraries (mesh CGAL downsampler turtle)
add_dependencies (mesh CGAL GEOMETRY TURTLE)
target_compile_definitions (mesh PRIVATE
    "-DINSTALL_PREFIX=\"${CMAKE_INSTALL_PREFIX}\"")
install (TARGETS mesh DESTINATION lib)


# The optimistic geometry
add_library (opti SHARED src/opti.c)
target_link_libraries (opti downsampler turtle)
add_dependencies (opti GEOMETRY TURTLE)
target_compile_definitions (opti PRIVATE
    "-DINSTALL_PREFIX=\"${CMAKE_INSTALL_PREFIX}\"")
install (TARGETS opti DESTINATION lib)


# The generic geometry
add_custom_target (GEOMETRY SOURCES src/geometry.h)
install (FILES include/geometry.h DESTINATION include)
install (FILES share/luajit-${LUAJIT_VERSION}/geometry.lua DESTINATION share/luajit-${LUAJIT_VERSION})


# The view helper
add_library (view SHARED src/view.c src/view.h)
target_link_libraries (view turtle)
add_dependencies (view TURTLE)
install (TARGETS view DESTINATION lib)
install (FILES src/view.h DESTINATION include)
install (FILES share/luajit-${LUAJIT_VERSION}/view.lua DESTINATION share/luajit-${LUAJIT_VERSION})


# The generic ray tracer algorithm
add_library (tracer SHARED src/tracer.c src/tracer.h)
target_link_libraries (tracer tictac)
install (TARGETS tracer DESTINATION lib)
install (FILES src/tracer.h DESTINATION include)
install (FILES share/luajit-${LUAJIT_VERSION}/tracer.lua DESTINATION share/luajit-${LUAJIT_VERSION})


# The generic Monte Carlo simulator
add_library (simulator SHARED src/simulator.c src/simulator.h)
target_link_libraries (simulator pumas tictac turtle)
add_dependencies (simulator PUMAS TURTLE)
target_compile_definitions (simulator PRIVATE
    "-DINSTALL_PREFIX=\"${CMAKE_INSTALL_PREFIX}\"")
install (TARGETS simulator DESTINATION lib)
install (FILES src/simulator.h DESTINATION include)
install (FILES share/luajit-${LUAJIT_VERSION}/simulator.lua DESTINATION share/luajit-${LUAJIT_VERSION})


# The clock utility
add_library (tictac SHARED src/tictac.c src/tictac.h)
install (TARGETS tictac DESTINATION lib)


# The downsampler utility
add_library (downsampler SHARED src/downsampler.c src/downsampler.h)
target_link_libraries (downsampler turtle)
add_dependencies (downsampler TURTLE)
install (TARGETS downsampler DESTINATION lib)


# Install the TURTLE library
ExternalProject_Add(TURTLE
    URL https://github.com/niess/turtle/archive/v0.7.tar.gz
    CONFIGURE_COMMAND cmake ../TURTLE -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DTURTLE_USE_GEOTIFF=OFF
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)


# Install the CGAL library
ExternalProject_Add(GMP
    URL https://gmplib.org/download/gmp/gmp-6.1.2.tar.xz
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ./configure --prefix=${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)

ExternalProject_Add(MPFR
    URL https://www.mpfr.org/mpfr-current/mpfr-4.0.2.tar.xz
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ./configure --prefix=${CMAKE_INSTALL_PREFIX} --with-gmp=${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)
add_dependencies (MPFR GMP)

set (__env GMP_DIR=${CMAKE_INSTALL_PREFIX} MPFR_DIR=${CMAKE_INSTALL_PREFIX})
ExternalProject_Add(CGAL
    URL https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.13.1/CGAL-4.13.1.tar.xz
    CONFIGURE_COMMAND ${__env} cmake ../CGAL -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DCMAKE_INSTALL_RPATH=${CMAKE_INSTALL_RPATH}
    COMMAND cd ${CMAKE_INSTALL_PREFIX} && ln -fs lib lib64
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    COMMAND rm -f ${CMAKE_INSTALL_PREFIX}/bin/cgal_create_CMakeLists
    COMMAND rm -f ${CMAKE_INSTALL_PREFIX}/bin/cgal_create_cmake_script
    COMMAND rm -f ${CMAKE_INSTALL_PREFIX}/bin/cgal_make_macosx_app
)
add_dependencies (CGAL GMP MPFR)


# Fetch the Embree library
ExternalProject_Add(EMBREE
    URL https://github.com/embree/embree/releases/download/v3.5.2/embree-3.5.2.x86_64.linux.tar.gz
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND cp lib/libembree3.so.3 ${CMAKE_INSTALL_PREFIX}/lib
    COMMAND cp lib/libtbbmalloc.so.2 ${CMAKE_INSTALL_PREFIX}/lib
    COMMAND cp lib/libtbb.so.2 ${CMAKE_INSTALL_PREFIX}/lib
    COMMAND cp -r include/embree3 ${CMAKE_INSTALL_PREFIX}/include
    COMMAND cd ${CMAKE_INSTALL_PREFIX}/lib && ln -fs libembree3.so.3 libembree3.so
    COMMAND cd ${CMAKE_INSTALL_PREFIX}/lib && ln -fs libtbb.so.2 libtbb.so
    COMMAND cd ${CMAKE_INSTALL_PREFIX}/lib && ln -fs libtbbmalloc.so.2 libtbbmalloc.so
)


# Install the (patched) PUMAS library
ExternalProject_Add(PUMAS
    URL https://github.com/niess/pumas/archive/v0.14.tar.gz
    CONFIGURE_COMMAND cmake ../PUMAS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    COMMAND cd ../PUMAS/src && patch -i ${CMAKE_SOURCE_DIR}/src/patch/pumas.patch pumas.c
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)


# Build LuaJIT
set(__env PREFIX=${CMAKE_INSTALL_PREFIX} BUILDMODE="static")
ExternalProject_Add(LuaJIT
    GIT_REPOSITORY "https://github.com/LuaJIT/LuaJIT.git"
    GIT_TAG "v${LUAJIT_VERSION}"
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE TRUE
    BUILD_COMMAND make ${__env}
    INSTALL_COMMAND make ${__env} install
)


# Argparse package for Lua
ExternalProject_Add(ARGPARSE
    URL https://github.com/mpeterv/argparse/archive/0.6.0.tar.gz
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND cp src/argparse.lua ${CMAKE_INSTALL_PREFIX}/share/lua/5.1
)
add_dependencies (ARGPARSE LuaJIT)


# Installer for the topography data
if(NOT EXISTS ${CMAKE_INSTALL_PREFIX}/share/topography/tianshan.png)
    ExternalProject_Add(TOPOGRAPHY_DATA
        URL https://github.com/niess/turtle-perfs/releases/download/topography/topography.tgz
        BUILD_IN_SOURCE TRUE
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND mkdir -p ${CMAKE_INSTALL_PREFIX}/share/topography
        COMMAND mv egm96.png pdd.png tianshan.png ${CMAKE_INSTALL_PREFIX}/share/topography
    )
endif ()


# Extra optional rules for Geant4
option (USE_GEANT4 "Build Geant4 tests." OFF)
if (${USE_GEANT4})
        include (cmake/geant4.cmake)
endif ()
