# Include Geant4
find_package (Geant4 REQUIRED)
include (${Geant4_USE_FILE})


# Build the G4Turtle wrapper
ExternalProject_Add(G4TURTLE
    GIT_REPOSITORY "https://github.com/niess/turtle-geant4.git"
    CONFIGURE_COMMAND cmake ../G4TURTLE -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DTURTLE_PREFIX=${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)
add_dependencies (G4TURTLE TURTLE)


# Build the geant4 processor
add_library (geant4 SHARED src/geant4.cc src/geant4.hh)
target_link_libraries (geant4 downsampler G4turtle ${Geant4_LIBRARIES} tictac)
target_include_directories (geant4 PRIVATE ${Geant4_INCLUDE_DIR})
target_compile_definitions (geant4 PRIVATE
    "-DINSTALL_PREFIX=\"${CMAKE_INSTALL_PREFIX}\"")
add_dependencies (geant4 G4TURTLE)

install (TARGETS geant4 DESTINATION lib)
install (FILES src/geant4.h DESTINATION include)
