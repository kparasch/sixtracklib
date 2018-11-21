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

  st_Buffer* input_pb = st_Buffer_new_from_file( "particles.buffer" );
  SIXTRL_ASSERT( input_pb != SIXTRL_NULLPTR );

  st_Particles const* in_particles =
      st_Particles_buffer_get_const_particles( input_pb, 0u );

  st_Particles_print_out( in_particles );




}

/* end: examples/c99/track_io.c */
