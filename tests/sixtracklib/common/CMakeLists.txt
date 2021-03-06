# sixtracklib/tests/sixtracklib/common/CMakeLists.txt

add_subdirectory( buffer )
add_subdirectory( control )
add_subdirectory( beam_elements )
add_subdirectory( track )

# -----------------------------------------------------------------------------
# ---- Generate testdata:

if( GTEST_FOUND )
    set( SIXTRACKL_COMMON_TESTS_DEPEND_ON_TRACKING_TESTDATA )

    set( UNIT_TEST_TARGETS )
    set( SIXTRACKL_TEST_LIBRARIES ${SIXTRACKL_TEST_LIBRARIES}
         ${SIXTRACKL_GTEST_LIBRARIES} )

    # --------------------------------------------------------------------------
    # test_tools_c99:

    add_executable( test_tools_common_c99 test_tools_c99.cpp )
    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} test_tools_common_c99 )
    add_test( C99_CommonToolTests test_tools_common_c99 )

    # ==========================================================================
    # test_gpu_kernel_tools_common_c99:

    add_executable( test_gpu_kernel_tools_common_c99
                    test_gpu_kernel_tools_c99.cpp )

    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS}
         test_gpu_kernel_tools_common_c99 )

    add_test( C99_CommonGpuKernelTools
              test_gpu_kernel_tools_common_c99 )

    # ==========================================================================
    # test_particles_common_c99:

    add_executable( test_particles_common_c99     test_particles_c99.cpp )
    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS}   test_particles_common_c99 )
    add_test( C99_CommonParticlesTests            test_particles_common_c99 )

    # --------------------------------------------------------------------------
    # test_particles_common_cxx:

    add_executable( test_particles_common_cxx     test_particles_cxx.cpp )
    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS}   test_particles_common_cxx )
    add_test( CXX_CommonParticlesTests            test_particles_common_cxx )

    # ==========================================================================
    # test_particles_aperture_check_enabled_common_c99:

    add_executable( test_particles_aperture_check_enabled_common_c99
                    test_particles_aperture_c99.cpp )

    target_compile_options( test_particles_aperture_check_enabled_common_c99
        BEFORE PRIVATE -DSIXTRL_APERTURE_CHECK_AT_DRIFT=2
                       -DSIXTRL_APERTURE_X_LIMIT=1.0
                       -DSIXTRL_APERTURE_Y_LIMIT=1.0 )

    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS}
         test_particles_aperture_check_enabled_common_c99 )

    add_test( C99_CommonParticlesApertureEnabledTests
              test_particles_aperture_check_enabled_common_c99 )

    # ------------------------------------------------------------------------
    # test_particles_aperture_check_disabled_common_c99:

    add_executable( test_particles_aperture_check_disabled_common_c99
                    test_particles_aperture_c99.cpp )

    target_compile_options( test_particles_aperture_check_disabled_common_c99
        BEFORE PRIVATE -DSIXTRL_APERTURE_CHECK_AT_DRIFT=0 )

    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS}
         test_particles_aperture_check_disabled_common_c99 )

    add_test( C99_CommonParticlesApertureDisabledTests
              test_particles_aperture_check_disabled_common_c99 )

    # --------------------------------------------------------------------------
    # test_elem_by_elem_config_c99:

    add_executable( test_elem_by_elem_config_c99
                    test_elem_by_elem_config_c99.cpp )

    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} test_elem_by_elem_config_c99 )
    add_test( C99_CommonElemByElemConfigTests   test_elem_by_elem_config_c99 )

    # --------------------------------------------------------------------------
    # test_faddeeva_errf:

#     add_executable( test_faddeeva_errf test_faddeeva_errf.cpp )
#     set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} test_faddeeva_errf )
#     add_test( CommonFaddeevaErrfnTests test_faddeeva_errf )

    # --------------------------------------------------------------------------
    # test_beam_beam_element_6d:

#     add_executable( test_beam_beam_element_6d test_beam_beam_element_6d.cpp )
#     set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} test_beam_beam_element_6d )
#     add_test( CommonBeamBeamElement6dTests, test_beam_beam_element_6d )

    # --------------------------------------------------------------------------
    # test_beam_elements:

#     add_executable( test_beam_elements_common test_beam_elements.cpp )
#     set(   UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} test_beam_elements_common )
#     add_test( CommonBeamElementsTests test_beam_elements_common )

    # ==========================================================================
    # test_output_buffer_common_c99:

    add_executable( test_output_buffer_common_c99 test_output_buffer_c99.cpp )
    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} test_output_buffer_common_c99 )
    add_test( C99_CommonOutputBufferTests test_output_buffer_common_c99 )

    # ==========================================================================
    # test_track_common_c99:

    add_executable( test_track_common_c99 test_track_c99.cpp )
    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} test_track_common_c99 )
    add_test( C99_CommonTrackTests test_track_common_c99 )

    # ==========================================================================
    # test_track_line_common_c99:

    add_executable( test_track_line_common_c99 test_track_line_c99.cpp )
    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} test_track_line_common_c99 )
    add_test( C99_CommonTrackLineTests test_track_line_common_c99 )

    # ==========================================================================
    # test_track_job_common_c99:

#     add_executable( test_track_job_common_c99 test_track_job_c99.cpp )
#     set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} test_track_job_common_c99 )
#     add_test( C99_CommonTrackJobTests, test_track_job_common_c99 )

    # --------------------------------------------------------------------------
    # test_track_job_common_cxx:

    add_executable( test_track_job_common_cxx test_track_job_cxx.cpp )
    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} test_track_job_common_cxx )
    add_test( CXX_CommonTrackJobTests, test_track_job_common_cxx )

    # ==========================================================================
    # test_track_job_cpu_common_c99:

    add_executable( test_track_job_cpu_common_c99 test_track_job_cpu_c99.cpp )

    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} test_track_job_cpu_common_c99 )
    add_test( C99_CommonTrackJobCpuTests test_track_job_cpu_common_c99 )

    # --------------------------------------------------------------------------
    # test_track_job_cpu_cxx:

    add_executable( test_track_job_cpu_common_cxx test_track_job_cpu_cxx.cpp )

    set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} test_track_job_cpu_common_cxx )
    add_test( CXX_CommonTrackJobCpuTests test_track_job_cpu_common_cxx )

    # --------------------------------------------------------------------------

    if( SIXTRACKL_ENABLE_AUTOVECTORIZATION )
        add_executable( test_track_common_autovec_c99 test_track_c99.cpp )
        set( UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS}
             test_track_common_autovec_c99 )

        target_compile_options( test_track_common_autovec_c99 BEFORE PRIVATE
            -ftree-vectorize -ftree-vectorizer-verbose=6 -fopt-info-loop
            -fno-fast-math --param vect-max-version-for-alias-checks=150 )

        add_test( C99_CommonTrackTestsAutoVec test_track_common_autovec_c99 )

    endif()

    # *************************************************************************
    # Set all properties:

    set_property(
        TARGET ${UNIT_TEST_TARGETS}
        APPEND PROPERTY INCLUDE_DIRECTORIES
        $<BUILD_INTERFACE:${SIXTRACKL_GTEST_INCLUDE_DIRS}>
    )

    set_property( TARGET ${UNIT_TEST_TARGETS}
                  PROPERTY LINK_LIBRARIES ${SIXTRACKL_TEST_LIBRARIES} )

    set_property( TARGET ${UNIT_TEST_TARGETS} PROPERTY CXX_STANDARD 11 )
    set_property( TARGET ${UNIT_TEST_TARGETS} PROPERTY CXX_STANDARD_REQUIRED ON )

    set_property( TARGET ${UNIT_TEST_TARGETS}
                  APPEND PROPERTY COMPILE_OPTIONS
                  ${SIXTRACKLIB_CXX_FLAGS} ${SIXTRACKLIB_CPU_FLAGS} )

endif()

# end: tests/sixtracklib/common/CMakeLists.txt
