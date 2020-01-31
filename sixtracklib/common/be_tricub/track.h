#ifndef SIXTRACKLIB_COMMON_BE_TRICUB_TRACK_C99_H__
#define SIXTRACKLIB_COMMON_BE_TRICUB_TRACK_C99_H__

#if !defined( SIXTRL_NO_INCLUDES )
    #include "sixtracklib/common/definitions.h"
    #include "sixtracklib/common/track/definitions.h"
    #include "sixtracklib/common/internal/beam_elements_defines.h"
    #include "sixtracklib/common/internal/objects_type_id.h"
    #include "sixtracklib/common/particles.h"
#endif /* !defined( SIXTRL_NO_INCLUDES ) */

#if !defined( _GPUCODE ) && defined( __cplusplus )
extern "C" {
#endif /* !defined(  _GPUCODE ) && defined( __cplusplus ) */

struct NS(TriCub);

SIXTRL_STATIC SIXTRL_FN NS(track_status_t) NS(Track_particle_tricub)(
    SIXTRL_PARTICLE_ARGPTR_DEC NS(Particles)* SIXTRL_RESTRICT particles,
    NS(particle_num_elements_t) const particle_index,
    SIXTRL_BE_ARGPTR_DEC const struct NS(TriCub) *const SIXTRL_RESTRICT tricub );

SIXTRL_STATIC SIXTRL_FN void NS(tricub_construct_b1_vector)(
    SIXTRL_BUFFER_DATAPTR_DEC const NS(TriCubData) *const SIXTRL_RESTRICT tricub_data,
    NS(be_tricub_int_t) const ix, NS(be_tricub_int_t) const iy, NS(be_tricub_int_t) const iz, 
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* SIXTRL_RESTRICT b_vector);

SIXTRL_STATIC SIXTRL_FN void NS(tricub_construct_bn_vector)(
    SIXTRL_BUFFER_DATAPTR_DEC const NS(TriCubData) *const SIXTRL_RESTRICT tricub_data,
    NS(be_tricub_int_t) const ix, NS(be_tricub_int_t) const iy, NS(be_tricub_int_t) const iz, 
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* SIXTRL_RESTRICT b_vector,
    NS(be_tricub_int_t) const l);//ow, NS(be_tricub_int_t) const high);

SIXTRL_STATIC SIXTRL_FN void NS(tricub_construct_b_vector)(
    SIXTRL_BUFFER_DATAPTR_DEC const NS(TriCubData) *const SIXTRL_RESTRICT tricub_data,
    NS(be_tricub_int_t) const ix, NS(be_tricub_int_t) const iy, NS(be_tricub_int_t) const iz, 
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* SIXTRL_RESTRICT b_vector);

#if !defined( _GPUCODE ) && defined( __cplusplus )
}
#endif /* !defined(  _GPUCODE ) && defined( __cplusplus ) */

/* ========================================================================= */
/* =====        Implementation of Inline functions and methods         ===== */
/* ========================================================================= */

#if !defined( SIXTRL_NO_SYSTEM_INCLUDES )
    #include <stddef.h>
    #include <stdint.h>
    #include <stdlib.h>
    #include <math.h>
#endif /* #if !defined( SIXTRL_NO_SYSTEM_INCLUDES ) */

#if !defined( SIXTRL_NO_INCLUDES )
    #include "sixtracklib/common/constants.h"
    #include "sixtracklib/common/particles.h"
    #include "sixtracklib/common/be_tricub/be_tricub.h"
    #include "sixtracklib/common/be_tricub/coefficients.h"
#endif /* !defined( SIXTRL_NO_INCLUDES ) */

#if !defined( _GPUCODE ) && defined( __cplusplus )
extern "C" {
#endif /* !defined(  _GPUCODE ) && defined( __cplusplus ) */

SIXTRL_INLINE NS(track_status_t) NS(Track_particle_tricub)(
    SIXTRL_PARTICLE_ARGPTR_DEC NS(Particles)* SIXTRL_RESTRICT particles,
    NS(particle_num_elements_t) const ii,
    SIXTRL_BE_ARGPTR_DEC const struct NS(TriCub) *const SIXTRL_RESTRICT tricub )
{
    typedef NS(be_tricub_real_t) real_t;
    typedef NS(be_tricub_int_t)  int_t;

    SIXTRL_BUFFER_DATAPTR_DEC NS(TriCubData) const* tricub_data =
        NS(TriCub_const_data)( tricub );

    real_t const length = NS(TriCub_length)( tricub );
    
    real_t const x_shift = NS(TriCub_x_shift)( tricub );
    real_t const y_shift = NS(TriCub_y_shift)( tricub );
    real_t const z_shift = NS(TriCub_tau_shift)( tricub );

    // method = 1 -> Finite Differences for derivatives (Do not use)
    // method = 2 -> Exact derivatives
    // method = 3 -> Exact derivatives and mirrored in X, Y
                 
    // real_t const inv_dx = 1./( NS(TriCubData_dx)( tricub_data ) );
    // real_t const inv_dy = 1./( NS(TriCubData_dy)( tricub_data ) );
    // real_t const inv_dz = 1./( NS(TriCubData_dz)( tricub_data ) );

    real_t const x0 = NS(TriCubData_x0)( tricub_data );
    real_t const y0 = NS(TriCubData_y0)( tricub_data );
    real_t const z0 = NS(TriCubData_z0)( tricub_data );

    real_t const zeta  = NS(Particles_get_zeta_value)( particles, ii );
    real_t const rvv   = NS(Particles_get_rvv_value)( particles, ii );
    real_t const beta0 = NS(Particles_get_beta0_value)( particles, ii );

    real_t const x = NS(Particles_get_x_value)( particles, ii );
    real_t const y = NS(Particles_get_y_value)( particles, ii );
    real_t const z = zeta / ( beta0 * rvv);

    real_t const fx = ( (x - x_shift ) - x0 )/( NS(TriCubData_dx)( tricub_data ) );
    real_t const fy = ( (y - y_shift ) - y0 )/( NS(TriCubData_dy)( tricub_data ) );
    real_t const fz = ( (z - z_shift ) - z0 )/( NS(TriCubData_dz)( tricub_data ) );

    real_t const sign_x = ( NS(TriCubData_mirror_x)( tricub_data ) == 1 && fx < 0.0 ) ? -1. : 1.;
    real_t const sign_y = ( NS(TriCubData_mirror_y)( tricub_data ) == 1 && fy < 0.0 ) ? -1. : 1.;
    real_t const sign_z = ( NS(TriCubData_mirror_z)( tricub_data ) == 1 && fz < 0.0 ) ? -1. : 1.;

    real_t const sfx = sign_x * fx;
    real_t const sfy = sign_y * fy;
    real_t const sfz = sign_z * fz;
                 
    real_t const ixf = floor(sfx);
    real_t const iyf = floor(sfy);
    real_t const izf = floor(sfz);

    int_t const ix = (int_t)ixf;
    int_t const iy = (int_t)iyf;
    int_t const iz = (int_t)izf;

    real_t const xn = sfx - ixf;
    real_t const yn = sfy - iyf;
    real_t const zn = sfz - izf;

    // const int_t inside_box = 
    //   ( ( ( ix < 0 || ix > NS(TriCubData_nx)( tricub_data ) - 2 )   ||
    //       ( iy < 0 || iy > NS(TriCubData_ny)( tricub_data ) - 2 ) ) ||
    //       ( iz < 0 || iz > NS(TriCubData_nz)( tricub_data ) - 2 ) ) ? 0 : 1;

    // SIXTRL_ASSERT( inside_box == 1 );
    // Check that coordinates are inside bounding box
    SIXTRL_ASSERT( ix >= 0 && ix <= NS(TriCubData_nx)( tricub_data ) - 2 );
    SIXTRL_ASSERT( iy >= 0 && iy <= NS(TriCubData_ny)( tricub_data ) - 2 ); 
    SIXTRL_ASSERT( iz >= 0 && iz <= NS(TriCubData_nz)( tricub_data ) - 2 );
    // =========================================================================
    
    real_t kicks[3] = {0.};

    // kicks[0] = sign_x*(NS(TriCubData_dx)( tricub_data ) / length )*NS(TriCub_dipolar_kick_px)( tricub );
    // kicks[1] = sign_y*(NS(TriCubData_dy)( tricub_data ) / length )*NS(TriCub_dipolar_kick_py)( tricub );
    // kicks[2] = sign_z*(NS(TriCubData_dz)( tricub_data ) / length )*NS(TriCub_dipolar_kick_ptau)( tricub );

    //real_t b_vector[64];
    //NS(tricub_construct_bn_vector)(tricub_data, ix, iy, iz, b_vector, 0, 8);
    //NS(tricub_unscaled_kicks_from_b)(b_vector, xn, yn, zn, kicks);
    //NS(tricub_construct_bn_vector)(tricub_data, ix, iy, iz, b_vector, 2, 4);
    //NS(tricub_unscaled_kicks_from_b2)(b_vector, xn, yn, zn, kicks);
    //NS(tricub_construct_bn_vector)(tricub_data, ix, iy, iz, b_vector, 4, 6);
    //NS(tricub_unscaled_kicks_from_b3)(b_vector, xn, yn, zn, kicks);
    //NS(tricub_construct_bn_vector)(tricub_data, ix, iy, iz, b_vector, 6, 8);
    //NS(tricub_unscaled_kicks_from_b4)(b_vector, xn, yn, zn, kicks);
    
    real_t b_vector[8];
    NS(tricub_construct_bn_vector)(tricub_data, ix, iy, iz, b_vector, 0);
    NS(tricub_unscaled_kicks_from_b1)(b_vector, xn, yn, zn, kicks);
    NS(tricub_construct_bn_vector)(tricub_data, ix, iy, iz, b_vector, 1);
    NS(tricub_unscaled_kicks_from_b2)(b_vector, xn, yn, zn, kicks);
    NS(tricub_construct_bn_vector)(tricub_data, ix, iy, iz, b_vector, 2);
    NS(tricub_unscaled_kicks_from_b3)(b_vector, xn, yn, zn, kicks);
    NS(tricub_construct_bn_vector)(tricub_data, ix, iy, iz, b_vector, 3);
    NS(tricub_unscaled_kicks_from_b4)(b_vector, xn, yn, zn, kicks);
    NS(tricub_construct_bn_vector)(tricub_data, ix, iy, iz, b_vector, 4);
    NS(tricub_unscaled_kicks_from_b5)(b_vector, xn, yn, zn, kicks);
    NS(tricub_construct_bn_vector)(tricub_data, ix, iy, iz, b_vector, 5);
    NS(tricub_unscaled_kicks_from_b6)(b_vector, xn, yn, zn, kicks);
    NS(tricub_construct_bn_vector)(tricub_data, ix, iy, iz, b_vector, 6);
    NS(tricub_unscaled_kicks_from_b7)(b_vector, xn, yn, zn, kicks);
    NS(tricub_construct_bn_vector)(tricub_data, ix, iy, iz, b_vector, 7);
    NS(tricub_unscaled_kicks_from_b8)(b_vector, xn, yn, zn, kicks);

   // real_t b_vector[64];
   // NS(tricub_construct_b_vector)(tricub_data, ix, iy, iz, b_vector);

   // NS(tricub_unscaled_kicks_from_b)(b_vector, xn, yn, zn, kicks);


    // real_t const sign2_x = ( NS(TriCubData_mirror_x)( tricub_data ) == 1 && fx < 0.0 ) ? -1. : 1.;
    // real_t const sign2_y = ( NS(TriCubData_mirror_y)( tricub_data ) == 1 && fy < 0.0 ) ? -1. : 1.;
    // real_t const sign2_z = ( NS(TriCubData_mirror_z)( tricub_data ) == 1 && fz < 0.0 ) ? -1. : 1.;


    kicks[0] *= ( length /( NS(TriCubData_dx)( tricub_data ) ) );
    kicks[0] *= -sign_x;
    kicks[0] -= NS(TriCub_dipolar_kick_px)( tricub );

    kicks[1] *= ( length /( NS(TriCubData_dy)( tricub_data ) ) );
    kicks[1] *= -sign_y;
    kicks[1] -= NS(TriCub_dipolar_kick_py)( tricub );

    kicks[2] *= ( length /( NS(TriCubData_dz)( tricub_data ) ) );
    kicks[2] *= -sign_z;
    kicks[2] -= NS(TriCub_dipolar_kick_ptau)( tricub );

    real_t const q = NS(Particles_get_q0_value)( particles, ii ) * 
                     NS(Particles_get_charge_ratio_value)( particles, ii );
    real_t const p0c = NS(Particles_get_p0c_value)( particles, ii );
    real_t const energy_kick = q * ( p0c * kicks[2] );

    NS(Particles_add_to_px_value)( particles, ii, kicks[0] );
    NS(Particles_add_to_py_value)( particles, ii, kicks[1] );
    NS(Particles_add_to_energy_value)( particles, ii, energy_kick );


    /* how to mark a particle in a particle set as lost */
    /*
    NS(Particles_mark_as_lost_value)( particles, ii );
    */

    ( void )particles;
    ( void )ii;
    //( void )lookup_table_begin;
    //( void )lookup_table_size;

    return SIXTRL_TRACK_SUCCESS;
}

SIXTRL_INLINE void NS(tricub_construct_b1_vector)(
    SIXTRL_BUFFER_DATAPTR_DEC const NS(TriCubData) *const SIXTRL_RESTRICT tricub_data,
    NS(be_tricub_int_t) const ix, NS(be_tricub_int_t) const iy, NS(be_tricub_int_t) const iz, 
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* SIXTRL_RESTRICT b_vector)
{

    SIXTRL_BUFFER_DATAPTR_DEC NS(be_tricub_real_t) const* lookup_table_begin =
        NS(TriCubData_const_table_begin)( tricub_data );

    //NS(be_tricub_int_t) const lookup_table_size = NS(TriCubData_table_size)( tricub_data );
    NS(be_tricub_int_t) const nx = NS(TriCubData_nx)( tricub_data );
    NS(be_tricub_int_t) const ny = NS(TriCubData_ny)( tricub_data );
    //NS(be_tricub_int_t) const nz = NS(TriCubData_nz)( tricub_data );

    // NS(be_tricub_real_t) const dx = NS(TriCubData_dx)( tricub_data );
    // NS(be_tricub_real_t) const dy = NS(TriCubData_dy)( tricub_data );
    // NS(be_tricub_real_t) const dz = NS(TriCubData_dz)( tricub_data );

    // NS(be_tricub_real_t) const scale[8] = { 1., dx, dy, dz, 
    //                                       dx * dy, dx * dz, dy * dz, 
    //                                       (dx * dy) * dz };
    for(int l = 0; l < 4; l++)
    {
        b_vector[8 * l    ] = lookup_table_begin[ l + 8 * ( (ix  ) + nx * ( (iy  ) + ny * ( (iz  ) ) ) ) ] ;//* scale[l];
        b_vector[8 * l + 1] = lookup_table_begin[ l + 8 * ( (ix+1) + nx * ( (iy  ) + ny * ( (iz  ) ) ) ) ] ;//* scale[l];
        b_vector[8 * l + 2] = lookup_table_begin[ l + 8 * ( (ix  ) + nx * ( (iy+1) + ny * ( (iz  ) ) ) ) ] ;//* scale[l];
        b_vector[8 * l + 3] = lookup_table_begin[ l + 8 * ( (ix+1) + nx * ( (iy+1) + ny * ( (iz  ) ) ) ) ] ;//* scale[l];
        b_vector[8 * l + 4] = lookup_table_begin[ l + 8 * ( (ix  ) + nx * ( (iy  ) + ny * ( (iz+1) ) ) ) ] ;//* scale[l];
        b_vector[8 * l + 5] = lookup_table_begin[ l + 8 * ( (ix+1) + nx * ( (iy  ) + ny * ( (iz+1) ) ) ) ] ;//* scale[l];
        b_vector[8 * l + 6] = lookup_table_begin[ l + 8 * ( (ix  ) + nx * ( (iy+1) + ny * ( (iz+1) ) ) ) ] ;//* scale[l];
        b_vector[8 * l + 7] = lookup_table_begin[ l + 8 * ( (ix+1) + nx * ( (iy+1) + ny * ( (iz+1) ) ) ) ] ;//* scale[l];
    }

    return;
}

SIXTRL_INLINE void NS(tricub_construct_bn_vector)(
    SIXTRL_BUFFER_DATAPTR_DEC const NS(TriCubData) *const SIXTRL_RESTRICT tricub_data,
    NS(be_tricub_int_t) const ix, NS(be_tricub_int_t) const iy, NS(be_tricub_int_t) const iz, 
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* SIXTRL_RESTRICT b_vector, 
    NS(be_tricub_int_t) const l)//ow, NS(be_tricub_int_t) const high)
{

    SIXTRL_BUFFER_DATAPTR_DEC NS(be_tricub_real_t) const* lookup_table_begin =
        NS(TriCubData_const_table_begin)( tricub_data );

    //NS(be_tricub_int_t) const lookup_table_size = NS(TriCubData_table_size)( tricub_data );
    NS(be_tricub_int_t) const nx = NS(TriCubData_nx)( tricub_data );
    NS(be_tricub_int_t) const ny = NS(TriCubData_ny)( tricub_data );
//    NS(be_tricub_int_t) const nz = NS(TriCubData_nz)( tricub_data );

    // NS(be_tricub_real_t) const dx = NS(TriCubData_dx)( tricub_data );
    // NS(be_tricub_real_t) const dy = NS(TriCubData_dy)( tricub_data );
    // NS(be_tricub_real_t) const dz = NS(TriCubData_dz)( tricub_data );

    //NS(be_tricub_real_t) const scale[8] = { 1., dx, dy, dz, 
    //                                      dx * dy, dx * dz, dy * dz, 
    //                                      (dx * dy) * dz };

    b_vector[ 0 ] = lookup_table_begin[ l + 8 * ( (ix  ) + nx * ( (iy  ) + ny * ( (iz  ) ) ) ) ] ;//* scale[l];
    b_vector[ 1 ] = lookup_table_begin[ l + 8 * ( (ix+1) + nx * ( (iy  ) + ny * ( (iz  ) ) ) ) ] ;//* scale[l];
    b_vector[ 2 ] = lookup_table_begin[ l + 8 * ( (ix  ) + nx * ( (iy+1) + ny * ( (iz  ) ) ) ) ] ;//* scale[l];
    b_vector[ 3 ] = lookup_table_begin[ l + 8 * ( (ix+1) + nx * ( (iy+1) + ny * ( (iz  ) ) ) ) ] ;//* scale[l];
    b_vector[ 4 ] = lookup_table_begin[ l + 8 * ( (ix  ) + nx * ( (iy  ) + ny * ( (iz+1) ) ) ) ] ;//* scale[l];
    b_vector[ 5 ] = lookup_table_begin[ l + 8 * ( (ix+1) + nx * ( (iy  ) + ny * ( (iz+1) ) ) ) ] ;//* scale[l];
    b_vector[ 6 ] = lookup_table_begin[ l + 8 * ( (ix  ) + nx * ( (iy+1) + ny * ( (iz+1) ) ) ) ] ;//* scale[l];
    b_vector[ 7 ] = lookup_table_begin[ l + 8 * ( (ix+1) + nx * ( (iy+1) + ny * ( (iz+1) ) ) ) ] ;//* scale[l];

    // NS(be_tricub_int_t) const off = 8*low;
    // for(int l = low; l < high; l++)
    // {
    //     b_vector[ 0 + 8*l - off] = lookup_table_begin[ l + 8 * ( (ix  ) + nx * ( (iy  ) + ny * ( (iz  ) ) ) ) ] ;//* scale[l];
    //     b_vector[ 1 + 8*l - off] = lookup_table_begin[ l + 8 * ( (ix+1) + nx * ( (iy  ) + ny * ( (iz  ) ) ) ) ] ;//* scale[l];
    //     b_vector[ 2 + 8*l - off] = lookup_table_begin[ l + 8 * ( (ix  ) + nx * ( (iy+1) + ny * ( (iz  ) ) ) ) ] ;//* scale[l];
    //     b_vector[ 3 + 8*l - off] = lookup_table_begin[ l + 8 * ( (ix+1) + nx * ( (iy+1) + ny * ( (iz  ) ) ) ) ] ;//* scale[l];
    //     b_vector[ 4 + 8*l - off] = lookup_table_begin[ l + 8 * ( (ix  ) + nx * ( (iy  ) + ny * ( (iz+1) ) ) ) ] ;//* scale[l];
    //     b_vector[ 5 + 8*l - off] = lookup_table_begin[ l + 8 * ( (ix+1) + nx * ( (iy  ) + ny * ( (iz+1) ) ) ) ] ;//* scale[l];
    //     b_vector[ 6 + 8*l - off] = lookup_table_begin[ l + 8 * ( (ix  ) + nx * ( (iy+1) + ny * ( (iz+1) ) ) ) ] ;//* scale[l];
    //     b_vector[ 7 + 8*l - off] = lookup_table_begin[ l + 8 * ( (ix+1) + nx * ( (iy+1) + ny * ( (iz+1) ) ) ) ] ;//* scale[l];
    // }

    return;
}

SIXTRL_INLINE void NS(tricub_construct_b_vector)(
    SIXTRL_BUFFER_DATAPTR_DEC const NS(TriCubData) *const SIXTRL_RESTRICT tricub_data,
    NS(be_tricub_int_t) const ix, NS(be_tricub_int_t) const iy, NS(be_tricub_int_t) const iz, 
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* SIXTRL_RESTRICT b_vector)
{

    SIXTRL_BUFFER_DATAPTR_DEC NS(be_tricub_real_t) const* lookup_table_begin =
        NS(TriCubData_const_table_begin)( tricub_data );

    //NS(be_tricub_int_t) const lookup_table_size = NS(TriCubData_table_size)( tricub_data );
    NS(be_tricub_int_t) const nx = NS(TriCubData_nx)( tricub_data );
    NS(be_tricub_int_t) const ny = NS(TriCubData_ny)( tricub_data );
    NS(be_tricub_int_t) const nz = NS(TriCubData_nz)( tricub_data );

    NS(be_tricub_real_t) const dx = NS(TriCubData_dx)( tricub_data );
    NS(be_tricub_real_t) const dy = NS(TriCubData_dy)( tricub_data );
    NS(be_tricub_real_t) const dz = NS(TriCubData_dz)( tricub_data );

    NS(be_tricub_real_t) const scale[8] = { 1., dx, dy, dz, 
                                          dx * dy, dx * dz, dy * dz, 
                                          (dx * dy) * dz };
    for(int l = 0; l < 8; l++)
    {
        b_vector[8 * l    ] = lookup_table_begin[ (ix  ) + nx * ( (iy  ) + ny * ( (iz  ) + nz * l) ) ] * scale[l];
        b_vector[8 * l + 1] = lookup_table_begin[ (ix+1) + nx * ( (iy  ) + ny * ( (iz  ) + nz * l) ) ] * scale[l];
        b_vector[8 * l + 2] = lookup_table_begin[ (ix  ) + nx * ( (iy+1) + ny * ( (iz  ) + nz * l) ) ] * scale[l];
        b_vector[8 * l + 3] = lookup_table_begin[ (ix+1) + nx * ( (iy+1) + ny * ( (iz  ) + nz * l) ) ] * scale[l];
        b_vector[8 * l + 4] = lookup_table_begin[ (ix  ) + nx * ( (iy  ) + ny * ( (iz+1) + nz * l) ) ] * scale[l];
        b_vector[8 * l + 5] = lookup_table_begin[ (ix+1) + nx * ( (iy  ) + ny * ( (iz+1) + nz * l) ) ] * scale[l];
        b_vector[8 * l + 6] = lookup_table_begin[ (ix  ) + nx * ( (iy+1) + ny * ( (iz+1) + nz * l) ) ] * scale[l];
        b_vector[8 * l + 7] = lookup_table_begin[ (ix+1) + nx * ( (iy+1) + ny * ( (iz+1) + nz * l) ) ] * scale[l];
    }

    return;
}

#if !defined( _GPUCODE ) && defined( __cplusplus )
}
#endif /* !defined(  _GPUCODE ) && defined( __cplusplus ) */

#endif /* SIXTRACKLIB_COMMON_BE_TRICUB_TRACK_C99_H__ */

/*end: sixtracklib/common/be_tricub/track.h */
