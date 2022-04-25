
set(HEADERS
    cpp/src/alglibinternal.h
    cpp/src/alglibmisc.h
    cpp/src/ap.h
    cpp/src/dataanalysis.h
    cpp/src/diffequations.h
    cpp/src/fasttransforms.h
    cpp/src/integration.h
    cpp/src/interpolation.h
    cpp/src/linalg.h
    cpp/src/optimization.h
    cpp/src/solvers.h
    cpp/src/specialfunctions.h
    cpp/src/statistics.h
    cpp/src/stdafx.h
)

set(SOURCES
    cpp/src/alglibinternal.cpp
    cpp/src/alglibmisc.cpp
    cpp/src/ap.cpp
    cpp/src/dataanalysis.cpp
    cpp/src/diffequations.cpp
    cpp/src/fasttransforms.cpp
    cpp/src/integration.cpp
    cpp/src/interpolation.cpp
    cpp/src/linalg.cpp
    cpp/src/optimization.cpp
    cpp/src/solvers.cpp
    cpp/src/specialfunctions.cpp
    cpp/src/statistics.cpp
)

# Set include directories.
# include_directories()

# Disable stupid MSVC compiler warnings.
# add_definitions(-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS)

set(CMAKE_DEBUG_POSTFIX "d")

add_library(alglib STATIC
    ${HEADERS}
    ${SOURCES}
)

# propnet has the following dependencies.
#add_dependencies(propnet
#)

# Set link libraries.
#target_link_libraries(propnet
#    ${Boost_LIBRARIES}
#)

# install(TARGETS propnet
#         RUNTIME DESTINATION bin
#         LIBRARY DESTINATION lib
#         ARCHIVE DESTINATION lib)
#
# #install(DIRECTORY include
#         DESTINATION .
#         PATTERN "*.hpp")
