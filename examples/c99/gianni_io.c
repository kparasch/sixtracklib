#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "sixtracklib/testlib.h"
#include "sixtracklib/sixtracklib.h"

int main( int argc, char* argv[] )
{

  printf("Hello world!\n");

  // Load particles
  st_Buffer* input_pb = st_Buffer_new_from_file( "particles.buffer" );
  SIXTRL_ASSERT( input_pb != SIXTRL_NULLPTR );

  st_Particles* particles =
      st_Particles_buffer_get_particles( input_pb, 0u );

  st_Particles_print_out( particles );





  st_Buffer* outp_buffer = st_Buffer_new( 0u );
  st_Particles_add_copy( outp_buffer, particles );

  st_Buffer_write_to_file( outp_buffer, "outp.buffer" );




}

/* end: examples/c99/track_io.c */
