# add_unit_test(<target> <sources>...)
#
# Adds test executable, <target>, built from <sources> and
# adds the test so that CTest will run it.
# Both the executable and the test will be named <target>.
function(add_unit_test target)
    # Add executable and libraries
    ADD_EXECUTABLE(${target} ${ARGN})
    # INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/unitTests)
    TARGET_LINK_LIBRARIES(${target} libflexishared)
    # Set properties of target
    SET_TARGET_PROPERTIES(${target} PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                               COMPILE_OPTIONS "-fPIC"
                                               CMAKE_POSITION_INDEPENDENT_CODE ON)
    SET_TARGET_PROPERTIES(${target} PROPERTIES LINKER_LANGUAGE Fortran)
    SET_TARGET_PROPERTIES(${target} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unitTests/lib"
                                               LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unitTests/lib"
                                               RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unitTests/bin")

    # Add this test to the global list of tests
    SET(UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} ${target} PARENT_SCOPE)
    # Add the test
    ADD_TEST(NAME ${target}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/unitTests
    COMMAND  ${CMAKE_BINARY_DIR}/unitTests/bin/${target})
endfunction()

# add_device_unit_test(<target> <sources>...)
#
# Adds test executable, <target>, built from <cpu_source> and
# <device_source> and adds the test so that CTest will run it.
# Both the executable and the test will be named <target>.
function(add_device_unit_test target cpu_source device_source)
    # Add executable and libraries
    ADD_EXECUTABLE(${target} ${cpu_source} ${device_source})
    # INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/unitTests)
    TARGET_LINK_LIBRARIES(${target} libflexishared)
    # Set properties of exe
    IF (LIBS_USE_ACCEL STREQUAL "AMD")
	    SET_SOURCE_FILES_PROPERTIES(${device_source} PROPERTIES LANGUAGE HIP)
    ENDIF()
    SET_TARGET_PROPERTIES(${target} PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                               CMAKE_POSITION_INDEPENDENT_CODE ON)
    TARGET_COMPILE_OPTIONS(${target} PRIVATE ${DEVICE_COMPILE_OPTIONS})
    SET_TARGET_PROPERTIES(${target} PROPERTIES LINKER_LANGUAGE Fortran)
    SET_TARGET_PROPERTIES(${target} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unitTests/lib"
                                               LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unitTests/lib"
                                               RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unitTests/bin")
    SET_TARGET_PROPERTIES(${target} PROPERTIES CUDA_ARCHITECTURES OFF)
    TARGET_INCLUDE_DIRECTORIES(${target} PRIVATE ${CMAKE_SOURCE_DIR}/src/device)

    # Add this test to the global list of tests
    SET(UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} ${target} PARENT_SCOPE)
    # Add the test
    ADD_TEST(NAME ${target}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/unitTests
    COMMAND  ${CMAKE_BINARY_DIR}/unitTests/bin/${target})
endfunction()

# CMake can always request position independent code
# SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Enable testing in general
ENABLE_TESTING()
# ... and make it parallel
IF (LIBS_USE_ACCEL STREQUAL "OFF")
    INCLUDE(ProcessorCount)
    PROCESSORCOUNT(N)
    SET(ENV{CTEST_PARALLEL_LEVEL} "${N}")
    MESSAGE(STATUS "Executing unit tests with $ENV{CTEST_PARALLEL_LEVEL} processes")
ELSE()
    # Actually don't, run in serial to prevent GPU tests clashing over resources
    SET(ENV{CTEST_PARALLEL_LEVEL} "1")
ENDIF()
    

# Add the unit tests
add_unit_test(ReadInTools ./unitTests/ReadInTools.f90)
add_unit_test(NodesAndWeights ./unitTests/NodesAndWeights.f90)
add_unit_test(Vandermonde ./unitTests/Vandermonde.f90)
add_unit_test(DerivativeMatrix ./unitTests/DerivativeMatrix.f90)
add_unit_test(ChangeBasis ./unitTests/ChangeBasis.f90)
add_unit_test(SurfInt ./unitTests/SurfInt.f90)
add_unit_test(ProlongToFace ./unitTests/ProlongToFace.f90)
add_unit_test(MatrixInverse    ./unitTests/MatrixInverse.f90)
add_unit_test(EOS ./unitTests/EOS.f90)
add_unit_test(ApplyDMatrix ./unitTests/ApplyDMatrix.f90)
add_unit_test(ApplyJacobian ./unitTests/ApplyJacobian.f90)
add_unit_test(CalcTimeStep ./unitTests/CalcTimeStep.f90)
add_unit_test(VAX ./unitTests/VAX.f90)

IF (LIBS_USE_ACCEL STREQUAL "OFF")
    # Add unit CPU versions of unit tests that also have .cu file based device tests
    # that can't be built without device support
    IF (FLEXI_SPLIT_DG)
        add_unit_test(SplitFlux ./unitTests/SplitFlux.f90)
    ENDIF()
    add_unit_test(EvalFlux ./unitTests/EvalFlux.f90)
    add_unit_test(Riemann ./unitTests/Riemann.f90)

ELSE()

    IF (LIBS_USE_ACCEL STREQUAL "NVIDIA")
        SET(ACCEL_LIB_TYPE 1)
        SET(DEVICE_LIBRARY "CUDA::cudart")
    ELSEIF (LIBS_USE_ACCEL STREQUAL "AMD")
        SET(ACCEL_LIB_TYPE 2)
        SET(DEVICE_LIBRARY "hip::host")
	    SET_SOURCE_FILES_PROPERTIES(./unitTests/*.cu PROPERTIES LANGUAGE CXX)
    ENDIF()

    IF (USE_APU_ARCH)
       ADD_DEFINITIONS(-DUSE_HYBRID)
    ENDIF()

    # Check for the presence of a supported accelerator using the test program in check_for_device.cpp.
    TRY_RUN (
        DEVICE_CHECK_RESULT
        DEVICE_COMPILE_RESULT
        ${CMAKE_BINARY_DIR}/unitTests/bin
        ${CMAKE_SOURCE_DIR}/unitTests/check_for_device.cpp
        COMPILE_DEFINITIONS -DUSE_ACCEL=${ACCEL_LIB_TYPE}
        LINK_LIBRARIES ${DEVICE_LIBRARY}
        COMPILE_OUTPUT_VARIABLE DEVICE_CHECK_COMPILE_OUT
        RUN_OUTPUT_VARIABLE DEVICE_CHECK_RUN_OUT
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/unitTests/bin
        CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${CMAKE_SOURCE_DIR}/src"
    )
    IF (NOT DEVICE_COMPILE_RESULT)
        MESSAGE (FATAL_ERROR "Compilation of device_check failed with output: ${DEVICE_CHECK_COMPILE_OUT}")
    ENDIF()

    MESSAGE ( STATUS "Device Check Result: ${DEVICE_CHECK_RUN_OUT}")

    SET(HAVE_DEVICE FALSE BOOL "Flag for whether a supported accelerator is present")
    IF (DEVICE_CHECK_RESULT)
        MESSAGE ( STATUS "Found supported devices, using device unit tests.")
        SET(HAVE_DEVICE TRUE)
    ENDIF()

    # If supported devices were found, we can run the tests
    IF (HAVE_DEVICE)
        add_device_unit_test(DeviceMemTools ./unitTests/DeviceMemTools.f90 ./unitTests/DeviceMemTools.cu)
        IF (FLEXI_SPLIT_DG)
            add_device_unit_test(SplitFlux ./unitTests/SplitFlux.f90 ./unitTests/SplitFlux.cu)
        ENDIF()
        add_device_unit_test(EvalFlux ./unitTests/EvalFlux.f90 ./unitTests/EvalFlux.cu)
        add_device_unit_test(Riemann ./unitTests/Riemann.f90 ./unitTests/Riemann.cu)
    ENDIF()
ENDIF()

# Perform unit tests when the last test has been compiled
ADD_CUSTOM_TARGET( all_tests ALL
                   DEPENDS ${UNIT_TEST_TARGETS})
ADD_CUSTOM_COMMAND(TARGET all_tests
                   COMMENT "Run unit tests"
                   POST_BUILD COMMAND ctest --output-on-failure -j${N}
                   WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
