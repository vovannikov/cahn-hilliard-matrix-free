CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

set(CMAKE_CXX_STANDARD 17)

FIND_PACKAGE(deal.II 9.0.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

IF(NOT DEAL_II_WITH_P4EST) # keep in one line
  MESSAGE(FATAL_ERROR "
Error! The deal.II library found at ${DEAL_II_PATH} was not configured with
    DEAL_II_WITH_P4EST = ON
It was set to OFF in your installation but is required for this program."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(cahn-hilliard-matrix-free)

# Define custom targets to easily switch the build type:
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug\nPlease call 'make' to build the debug project"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMENT "Switch CMAKE_BUILD_TYPE to Release\nPlease call 'make' to build the release project"
  )

SET(TARGET_SRC cahn_hilliard_impl.cc)
SET(TARGET "cahn_hilliard_impl")

ADD_EXECUTABLE(${TARGET} ${TARGET_SRC})

DEAL_II_SETUP_TARGET(${TARGET})

TARGET_INCLUDE_DIRECTORIES(${TARGET} PRIVATE ${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR})