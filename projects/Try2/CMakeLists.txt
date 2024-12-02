cmake_minimum_required(VERSION 3.10.0)

set(SRCS 
  src/main.cpp
  ../../utils/src/colormap.cpp
)

set(TEST_SRCS
  src/test-curv.cpp
)

add_subdirectory(../../deps/geometry-central deps/geometry-central)
add_subdirectory(../../deps/polyscope deps/polyscope)
add_subdirectory(../../deps/googletest deps/googletest)

# Change the executable name by changing "main" in the lines below.
add_executable(main "${SRCS}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main)
target_include_directories(main PUBLIC "C:/Users/dgrossma/Documents/GitHub/ElasticSim/deps/**")
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# add the args.hxx project which we use for command line args
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_link_libraries(main geometry-central polyscope)

#add_executable(test-curv "${TEST_SRCS}")
#target_include_directories(test-curv PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
#target_link_libraries(test-curv geometry-central polyscope gtest)
