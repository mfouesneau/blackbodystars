cmake_minimum_required(VERSION 3.0.0)
project(blackbodystars VERSION 0.1.0)
set (CMAKE_CXX_STANDARD 17)

# Set the output folder where your program will be created
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY build)

# generate the executables in the root directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")

# side load cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# with xtensor we need some math libraries
find_package(BLAS REQUIRED)
message( STATUS "BLAS found: ${blas_libraries}" )
find_package(LAPACK REQUIRED)
message( STATUS "LAPACK found: ${lapack_libraries}" )

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)

include(conanbuildinfo.cmake)  # Includes the contents of the conanbuildinfo.cmake file.
conan_basic_setup()  # Prepares the CMakeList.txt for Conan.

# add local includes
include_directories(
        "${PROJECT_SOURCE_DIR}/include"
        # "${PROJECT_SOURCE_DIR}/src"
        # "${PROJECT_SOURCE_DIR}/tests"
        )

# Define executables
# ------------------
add_executable(blackbodystars
        ${PROJECT_SOURCE_DIR}/src/main.cpp)

add_executable(test_cpr
        ${PROJECT_SOURCE_DIR}/src/cpr_test.cpp)


# Link the executables to external libraries
# -------------------------------------------
target_link_libraries(blackbodystars
        ${blas_libraries} ${lapack_libraries}
        ${CONAN_LIBS})

target_link_libraries(test_cpr
        ${blas_libraries} ${lapack_libraries}
        ${CONAN_LIBS})

# Where to install the targets --
install(TARGETS blackbodystars
        CONFIGURATIONS runtime
        RUNTIME DESTINATION bin
        )


# Tests
# --------------
# add testing libraries
# include(CTest)
enable_testing()
add_executable(test_main
        ${PROJECT_SOURCE_DIR}/tests/test_main.cpp)

target_link_libraries(test_main
        ${blas_libraries} ${lapack_libraries}
        ${CONAN_LIBS})

add_test(NAME example_tests
         COMMAND test_main)

# Documentation
# first we indicate the documentation build as an option
option(BUILD_DOC "Build documentation" ON)

# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
  if (CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
    # set input and output files
    set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doc_doxygen/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc_doxygen
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
  endif()
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)