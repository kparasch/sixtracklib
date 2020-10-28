#ifndef SIXTRACKLIB_COMMON_BE_TRICUB_COEFFICIENTS_C99_H__
#define SIXTRACKLIB_COMMON_BE_TRICUB_COEFFICIENTS_C99_H__

/*
#if !defined(SIXTRL_NO_SYSTEM_INCLUDES)
    #include <stddef.h>
    #include <math.h>
#endif */
/* !defined(SIXTRL_NO_SYSTEM_INCLUDES) */

#if !defined(SIXTRL_NO_INCLUDES)
    #include "sixtracklib/common/definitions.h"
    #include "sixtracklib/common/be_tricub/be_tricub.h"
#endif /* !defined(SIXTRL_NO_INCLUDES) */

#if !defined( _GPUCODE ) && defined( __cplusplus )
extern "C" {
#endif /* !defined( _GPUCODE ) && defined( __cplusplus ) */

SIXTRL_STATIC SIXTRL_FN void NS(tricub_construct_coefs)(
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t) const* SIXTRL_RESTRICT b,
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* coefs);

#if !defined( _GPUCODE ) && defined( __cplusplus )
}
#endif /* !defined( _GPUCODE ) && defined( __cplusplus ) */

/* ************************************************************************* */
/* ************************************************************************* */

#if !defined( _GPUCODE ) && defined( __cplusplus )
extern "C" {
#endif /* !defined( _GPUCODE ) && defined( __cplusplus ) */

SIXTRL_INLINE void NS(tricub_construct_coefs)(
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t) const* SIXTRL_RESTRICT b,
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* coefs)
{

    coefs[0] = b[0];
    coefs[1] = b[1] - coefs[0];
    coefs[2] = b[2] - coefs[0];
    coefs[3] = b[4] - coefs[0];
    coefs[4] = b[3] - ( coefs[0] - ( coefs[1] - coefs[2] ) );
    coefs[5] = b[5] - ( coefs[0] - ( coefs[1] - coefs[3] ) );
    coefs[6] = b[6] - ( coefs[0] - ( coefs[2] - coefs[2] ) );
    coefs[7] = b[7] - ( coefs[0] - ( coefs[1] - ( coefs[2] - ( coefs[3] - ( coefs[4] - ( coefs[5] - ( coefs[6] ) ) ) ) ) ) );

    return;
}


#if !defined( _GPUCODE ) && defined( __cplusplus )
}
#endif /* !defined( _GPUCODE ) && defined( __cplusplus ) */

#endif /* SIXTRACKLIB_COMMON_BE_TRICUB_COEFFICIENTS_C99_H__ */


