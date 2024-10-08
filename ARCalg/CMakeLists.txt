set(PackageName ARCalg)

project(${PackageName})

# Find GSL
find_package(GSL REQUIRED)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Core Graf Graf3d Hist)

# Include ROOT's CMake configuration
include(${ROOT_USE_FILE})

file(GLOB sources
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/src/Reconstruction.cpp
)

file(GLOB headers
    ${PROJECT_SOURCE_DIR}/include/*.h
    ${PROJECT_SOURCE_DIR}/src/Reconstruction.h
)

gaudi_add_module(${PackageName}
  SOURCES ${sources}
  LINK
  k4FWCore::k4FWCore
  k4FWCore::k4Interface
  Gaudi::GaudiKernel
  EDM4HEP::edm4hep
  DD4hep::DDCore
  DD4hep::DDRec
  GSL::gsl          # Link GSL
  GSL::gslcblas     # Link GSL CBLAS
  ${ROOT_LIBRARIES} # Link ROOT libraries
)

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(${PackageName} PROPERTIES PUBLIC_HEADER "${headers}")

file(GLOB scripts
  ${PROJECT_SOURCE_DIR}/test/*.py
)

file(COPY ${scripts} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test)

install(TARGETS ${PackageName}
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/@{CMAKE_PROJECT_NAME}" COMPONENT dev
)

install(FILES ${scripts} DESTINATION test)

# SET(test_name "test_runARCdigitizer")
# ADD_TEST(NAME ${test_name} COMMAND k4run test/runARCdigitizer.py)
# set_test_env(${test_name})
