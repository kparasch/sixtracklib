#define _USE_MATH_DEFINES 

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include <gtest/gtest.h>

#if defined( __NAMESPACE )
    #define __SAVED_NAMESPACE __NAMESPACE
    #undef  __NAMESPACE     
#endif /* !defiend( __NAMESPACE ) */

#if !defined( __NAMESPACE )
    #define __NAMESPACE st_    
#endif /* !defiend( __NAMESPACE ) */

#include "sixtracklib/_impl/definitions.h"
#include "sixtracklib/_impl/path.h"
#include "sixtracklib/_impl/testdata_files.h"

#include "sixtracklib/common/blocks.h"
#include "sixtracklib/common/beam_elements.h"
#include "sixtracklib/common/particles.h"
#include "sixtracklib/common/track.h"
#include "sixtracklib/common/details/random.h"

#include "sixtracklib/common/tests/test_particles_tools.h"
#include "sixtracklib/common/tests/test_track_tools.h"

#if defined( __SAVED_NAMESPACE )
    #undef __NAMESPACE
    #define __NAMESPACE __SAVED_NAMESPACE
#endif /* defined( __SAVED_NAMESPACE ) */

/* ************************************************************************* */
/* *****                TESTS FOR BEAM_ELEMENT DRIFT                   ***** */
/* ************************************************************************* */

/* ========================================================================= */
/* ====   Test using the drift-only testdata file and call              ==== */
/* ====   the minimal st_Track_drift_particle() function from track.h   ==== */
/* ====   manually -> i.e. do everything by hand, incl. elem_by_elem io ==== */
/* ====                                                                 ==== */
/* ====   Since this has to be repeated for each and every BeamElement, ==== */
/* ====   and to avoid code duplication even in the unit-tests, cf.     ==== */
/* ====   the sixtracklib/common/tests/test_track_tools.h file for the  ==== */
/* ====   declaration and definition of the helper function which       ==== */
/* ====   contains the whole unit-test!                                 ==== */
/* ========================================================================= */

TEST( CommonTrackTests, TrackDriftParticle )
{
    bool const success = 
        ::st_TestData_test_tracking_single_particle_over_specific_be_type< 
            st_Drift >( st_PATH_TO_TEST_TRACKING_DRIFT_DATA, 
                        st_Track_drift_particle );
        
    ASSERT_TRUE( success );
}

/* ========================================================================= */
/* ====   Test using the drift-only testdata file and call              ==== */
/* ====   st_Track_drift() function from track_api.h. Expected are same ==== */
/* ====   results as from st_Track_drift()                              ==== */
/* ====                                                                 ==== */
/* ====   Since this has to be repeated for each and every BeamElement, ==== */
/* ====   and to avoid code duplication even in the unit-tests, cf.     ==== */
/* ====   the sixtracklib/common/tests/test_track_tools.h file for the  ==== */
/* ====   declaration and definition of the helper function which       ==== */
/* ====   contains the whole unit-test!                                 ==== */
/* ========================================================================= */

TEST( CommonTrackTests, TrackDrift )
{
    bool const success =
        ::st_TestData_test_tracking_particles_over_specific_be_type< 
            st_Drift >( st_PATH_TO_TEST_TRACKING_DRIFT_DATA, 
                        st_Track_drift );
        
    ASSERT_TRUE( success );
}

/* ========================================================================= */
/* ====   Test using the drift-only testdata file and call              ==== */
/* ====   st_Track_beam_elements_particle() function from track_api.h.  ==== */
/* ====   Expected are same   results as from st_Track_drift_particle() ==== */
/* ====                                                                 ==== */
/* ====   Since this has to be repeated for each and every BeamElement, ==== */
/* ====   and to avoid code duplication even in the unit-tests, cf.     ==== */
/* ====   the sixtracklib/common/tests/test_track_tools.h and           ==== */
/* ====   sixtracklib/common/tests/test_track_tools.c files for the     ==== */
/* ====   declaration and definition of the helper functionm which      ==== */
/* ====   contains the whole unit-test, respectively!                   ==== */
/* ========================================================================= */

TEST( CommonTrackTests, TrackBeamElementParticleForDrift )
{
    bool const success = ::st_TestData_test_tracking_single_particle(
        st_PATH_TO_TEST_TRACKING_DRIFT_DATA );
    
    ASSERT_TRUE( success );
}

/* ========================================================================= */
/* ====   Test using the drift-only testdata file and call              ==== */
/* ====   st_Track_beam_elements_particle() function from track_api.h.  ==== */
/* ====   Expected are same   results as from st_Track_drift_particle() ==== */
/* ====                                                                 ==== */
/* ====   Since this has to be repeated for each and every BeamElement, ==== */
/* ====   and to avoid code duplication even in the unit-tests, cf.     ==== */
/* ====   the sixtracklib/common/tests/test_track_tools.h and           ==== */
/* ====   sixtracklib/common/tests/test_track_tools.c files for the     ==== */
/* ====   declaration and definition of the helper functionm which      ==== */
/* ====   contains the whole unit-test, respectively!                   ==== */
/* ========================================================================= */

TEST( CommonTrackTests, TrackBeamElementForDrift )
{
    bool const success = ::st_TestData_test_tracking_particles(
        st_PATH_TO_TEST_TRACKING_DRIFT_DATA );
    
    ASSERT_TRUE( success );
}

/* ************************************************************************* */
/* *****             TESTS FOR BEAM_ELEMENT DRIFT_EXACT                ***** */
/* ************************************************************************* */

/* ========================================================================= */
/* ====   Test using the drift_exact-only testdata file and call        ==== */
/* ====   the minimal st_Track_drift_particle() function from track.h   ==== */
/* ====   manually -> i.e. do everything by hand, incl. elem_by_elem io ==== */
/* ====                                                                 ==== */
/* ====   Since this has to be repeated for each and every BeamElement, ==== */
/* ====   and to avoid code duplication even in the unit-tests, cf.     ==== */
/* ====   the sixtracklib/common/tests/test_track_tools.h file for the  ==== */
/* ====   declaration and definition of the helper function which       ==== */
/* ====   contains the whole unit-test!                                 ==== */
/* ========================================================================= */

TEST( CommonTrackTests, TrackDriftExactParticle )
{
    bool const success = 
        ::st_TestData_test_tracking_single_particle_over_specific_be_type< 
            st_DriftExact >( st_PATH_TO_TEST_TRACKING_DRIFT_EXACT_DATA, 
                        st_Track_drift_exact_particle );
        
    ASSERT_TRUE( success );
}

/* ========================================================================= */
/* ====   Test using the drift_exact-only testdata file and call        ==== */
/* ====   st_Track_drift() function from track_api.h. Expected are same ==== */
/* ====   results as from st_Track_drift()                              ==== */
/* ====                                                                 ==== */
/* ====   Since this has to be repeated for each and every BeamElement, ==== */
/* ====   and to avoid code duplication even in the unit-tests, cf.     ==== */
/* ====   the sixtracklib/common/tests/test_track_tools.h file for the  ==== */
/* ====   declaration and definition of the helper function which       ==== */
/* ====   contains the whole unit-test!                                 ==== */
/* ========================================================================= */

TEST( CommonTrackTests, TrackDriftExact )
{
    bool const success =
        ::st_TestData_test_tracking_particles_over_specific_be_type< 
            st_DriftExact >( st_PATH_TO_TEST_TRACKING_DRIFT_EXACT_DATA, 
                        st_Track_drift_exact );
        
    ASSERT_TRUE( success );
}

/* ========================================================================= */
/* ====   Test using the drift_exact-only testdata file and call        ==== */
/* ====   st_Track_beam_elements_particle() function from track_api.h.  ==== */
/* ====   Expected are same   results as from st_Track_drift_particle() ==== */
/* ====                                                                 ==== */
/* ====   Since this has to be repeated for each and every BeamElement, ==== */
/* ====   and to avoid code duplication even in the unit-tests, cf.     ==== */
/* ====   the sixtracklib/common/tests/test_track_tools.h and           ==== */
/* ====   sixtracklib/common/tests/test_track_tools.c files for the     ==== */
/* ====   declaration and definition of the helper functionm which      ==== */
/* ====   contains the whole unit-test, respectively!                   ==== */
/* ========================================================================= */

TEST( CommonTrackTests, TrackBeamElementParticleForDriftExact )
{
    bool const success = ::st_TestData_test_tracking_single_particle(
        st_PATH_TO_TEST_TRACKING_DRIFT_EXACT_DATA );
    
    ASSERT_TRUE( success );
}

/* ========================================================================= */
/* ====   Test using the drift_exact-only testdata file and call        ==== */
/* ====   st_Track_beam_elements_particle() function from track_api.h.  ==== */
/* ====   Expected are same   results as from st_Track_drift_particle() ==== */
/* ====                                                                 ==== */
/* ====   Since this has to be repeated for each and every BeamElement, ==== */
/* ====   and to avoid code duplication even in the unit-tests, cf.     ==== */
/* ====   the sixtracklib/common/tests/test_track_tools.h and           ==== */
/* ====   sixtracklib/common/tests/test_track_tools.c files for the     ==== */
/* ====   declaration and definition of the helper functionm which      ==== */
/* ====   contains the whole unit-test, respectively!                   ==== */
/* ========================================================================= */

TEST( CommonTrackTests, TrackBeamElementForDriftExact )
{
    bool const success = ::st_TestData_test_tracking_particles(
        st_PATH_TO_TEST_TRACKING_DRIFT_EXACT_DATA );
    
    ASSERT_TRUE( success );
}

/* ************************************************************************* */
/* *****             TESTS FOR BEAM_ELEMENT MULTIPOLE                  ***** */
/* ************************************************************************* */

/* ========================================================================= */
/* ====   Test using the multipole-only testdata file and call          ==== */
/* ====   the minimal st_Track_drift_particle() function from track.h   ==== */
/* ====   manually -> i.e. do everything by hand, incl. elem_by_elem io ==== */
/* ====                                                                 ==== */
/* ====   Since this has to be repeated for each and every BeamElement, ==== */
/* ====   and to avoid code duplication even in the unit-tests, cf.     ==== */
/* ====   the sixtracklib/common/tests/test_track_tools.h file for the  ==== */
/* ====   declaration and definition of the helper function which       ==== */
/* ====   contains the whole unit-test!                                 ==== */
/* ========================================================================= */

TEST( CommonTrackTests, TrackMultiPoleParticle )
{
    bool const success = 
        ::st_TestData_test_tracking_single_particle_over_specific_be_type< 
            st_MultiPole >( st_PATH_TO_TEST_TRACKING_MULTIPOLE_DATA, 
                        st_Track_multipole_particle );
        
    ASSERT_TRUE( success );
}

/* ========================================================================= */
/* ====   Test using the multipole-only testdata file and call        ==== */
/* ====   st_Track_drift() function from track_api.h. Expected are same ==== */
/* ====   results as from st_Track_drift()                              ==== */
/* ====                                                                 ==== */
/* ====   Since this has to be repeated for each and every BeamElement, ==== */
/* ====   and to avoid code duplication even in the unit-tests, cf.     ==== */
/* ====   the sixtracklib/common/tests/test_track_tools.h file for the  ==== */
/* ====   declaration and definition of the helper function which       ==== */
/* ====   contains the whole unit-test!                                 ==== */
/* ========================================================================= */

TEST( CommonTrackTests, TrackMultiPole )
{
    bool const success =
        ::st_TestData_test_tracking_particles_over_specific_be_type< 
            st_MultiPole >( st_PATH_TO_TEST_TRACKING_MULTIPOLE_DATA, 
                        st_Track_multipole );
        
    ASSERT_TRUE( success );
}

/* ========================================================================= */
/* ====   Test using the multipole-only testdata file and call        ==== */
/* ====   st_Track_beam_elements_particle() function from track_api.h.  ==== */
/* ====   Expected are same   results as from st_Track_drift_particle() ==== */
/* ====                                                                 ==== */
/* ====   Since this has to be repeated for each and every BeamElement, ==== */
/* ====   and to avoid code duplication even in the unit-tests, cf.     ==== */
/* ====   the sixtracklib/common/tests/test_track_tools.h and           ==== */
/* ====   sixtracklib/common/tests/test_track_tools.c files for the     ==== */
/* ====   declaration and definition of the helper functionm which      ==== */
/* ====   contains the whole unit-test, respectively!                   ==== */
/* ========================================================================= */

TEST( CommonTrackTests, TrackBeamElementParticleForMultiPole )
{
    bool const success = ::st_TestData_test_tracking_single_particle(
        st_PATH_TO_TEST_TRACKING_MULTIPOLE_DATA );
    
    ASSERT_TRUE( success );
}

/* ========================================================================= */
/* ====   Test using the multipole-only testdata file and call        ==== */
/* ====   st_Track_beam_elements_particle() function from track_api.h.  ==== */
/* ====   Expected are same   results as from st_Track_drift_particle() ==== */
/* ====                                                                 ==== */
/* ====   Since this has to be repeated for each and every BeamElement, ==== */
/* ====   and to avoid code duplication even in the unit-tests, cf.     ==== */
/* ====   the sixtracklib/common/tests/test_track_tools.h and           ==== */
/* ====   sixtracklib/common/tests/test_track_tools.c files for the     ==== */
/* ====   declaration and definition of the helper functionm which      ==== */
/* ====   contains the whole unit-test, respectively!                   ==== */
/* ========================================================================= */

TEST( CommonTrackTests, TrackBeamElementForMultiPole )
{
    bool const success = ::st_TestData_test_tracking_particles(
        st_PATH_TO_TEST_TRACKING_MULTIPOLE_DATA );
    
    ASSERT_TRUE( success );
}

/* end: sixtracklib/common/tests/test_particles.cpp */
