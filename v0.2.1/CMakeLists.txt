
set(HEADERS
    bah.hpp
    blockvectors.hpp
    bvh.hpp
    cloth.hpp
    collision.hpp
    collisionutil.hpp
    conf.hpp
    constraint.hpp
    dde.hpp
    display.hpp
    displayplay.hpp
    displayphysics.hpp
    displayreplay.hpp
    displaytesting.hpp
    dynamicremesh.hpp
    geometry.hpp
    handle.hpp
    io.hpp
    magic.hpp
    mesh.hpp
    misc.hpp
    morph.hpp
    mot_parser.hpp
    nearobs.hpp
    obstacle.hpp
    opengl.hpp
    optimization.hpp
    physics.hpp
    plasticity.hpp
    popfilter.hpp
    proximity.hpp
    remesh.hpp
    runphysics.hpp
    separate.hpp
    separateobs.hpp
    simulation.hpp
    sparse.hpp
    spline.hpp
    sparse_solver.hpp
    strainlimiting.hpp
    # taucs.hpp
    tensormax.hpp
    timer.hpp
    transformation.hpp
    util.hpp
    vectors.hpp
    winport.hpp
    tensor.h 
    smpl.h
)

set(SOURCES
    auglag.cpp
    bah.cpp
    bvh.cpp
    cloth.cpp
    collision.cpp
    collisionutil.cpp
    conf.cpp
    constraint.cpp
    dde.cpp
    display.cpp
    displayplay.cpp
    displayphysics.cpp
    displayreplay.cpp
    displaytesting.cpp
    dynamicremesh.cpp
    geometry.cpp
    handle.cpp
    io.cpp
    lbfgs.cpp
    lsnewton.cpp
    magic.cpp
    main.cpp
    mesh.cpp
    misc.cpp
    morph.cpp
    mot_parser.cpp
    nearobs.cpp
    nlcg.cpp
    obstacle.cpp
    physics.cpp
    plasticity.cpp
    popfilter.cpp
    proximity.cpp
    remesh.cpp
    runphysics.cpp
    separate.cpp
    separateobs.cpp
    simulation.cpp
    sparse_solver.cpp
    spline.cpp
    strainlimiting.cpp
    # taucs.cpp
    tensormax.cpp
    timer.cpp
    transformation.cpp
    trustregion.cpp
    util.cpp
    vectors.cpp
    smpl.cpp
    tensor.cpp
)

#HYLC sources
file(GLOB_RECURSE HYLC_SRCS hylc/*.cpp)
file(GLOB_RECURSE HYLC_HDRS hylc/*.hpp hylc/*.h)

# Boost library.
set(Boost_USE_STATIC_LIBS ON)               # Use static libs.
add_definitions(-DBOOST_ALL_NO_LIB)         # Disable auto-linking in MSVC.
find_package(Boost 1.55.0 QUIET REQUIRED COMPONENTS system filesystem)

if (NOT WIN32)
    # png library.
    find_package(PNG QUIET REQUIRED)
    # zlib library.
    find_package(ZLIB QUIET REQUIRED)
    # GLUT library.
    find_package(GLUT QUIET REQUIRED)
endif()

# OpenGL library.
find_package(OpenGL QUIET REQUIRED)

# OpenCV library
find_package(OpenCV REQUIRED)

# Set include directories.
if (WIN32)
include_directories( SYSTEM # system to avoid warnings in libraries
    ../dependencies/alglib/cpp/src
    ../dependencies/jsoncpp/include
    # ../dependencies/fitpackpp/fitpackpp
    ../dependencies/eigen
    ../dependencies/freeglut/include
    ../dependencies/zlib
    ${Boost_INCLUDE_DIRS}
    ${OPENGL_INCLUDE_DIR}
    ${OpenCV_INCLUDE_DIRS}
)
else(WIN32)
include_directories( SYSTEM # system to avoid warnings in libraries
    ../dependencies/alglib/cpp/src
    ../dependencies/jsoncpp/include
    # ../dependencies/fitpackpp/fitpackpp
    ../dependencies/eigen
    ${Boost_INCLUDE_DIRS}
    ${OPENGL_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}
    ${ZLIB_INCLUDE_DIRS}
    ${PNG_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
)
endif(WIN32)

# Disable stupid MSVC compiler warnings.
if (WIN32)
    add_definitions(-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS)
endif()

# Set freeglut to be linked statically.
if (WIN32)
    add_definitions(-DFREEGLUT_STATIC)
    set(CMAKE_DEBUG_POSTFIX "d")
endif(WIN32)

add_executable(arcsim_0.2.1
    ${HEADERS}
    ${HYLC_HDRS}
    ${SOURCES}
    ${HYLC_SRCS}
)

if (WIN32)
# arcsim has the following dependencies.
add_dependencies(arcsim_0.2.1
    alglib
    json
    # fitpackpp
    freeglut_static
    zlibstatic
)
# Set link libraries.
target_link_libraries(arcsim_0.2.1
   alglib
   json
#    fitpackpp
   freeglut_static
   zlibstatic
   ${Boost_LIBRARIES}
   ${OPENGL_LIBRARIES}
   ${OpenCV_LIBRARIES}
)
else(WIN32)
# arcsim has the following dependencies.
add_dependencies(arcsim_0.2.1
    alglib
    json
    # fitpackpp
)
# Set link libraries.
target_link_libraries(arcsim_0.2.1
   alglib
   json
#    fitpackpp
   ${Boost_LIBRARIES}
   ${OPENGL_LIBRARIES}
   ${GLUT_LIBRARIES}
   ${ZLIB_LIBRARIES}
   ${PNG_LIBRARIES}
   ${OpenCV_LIBRARIES}
)
endif(WIN32)
