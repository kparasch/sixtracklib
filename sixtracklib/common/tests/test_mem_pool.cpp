#if !defined( __NAMESPACE )
    #define __NAMESPACE st_
    #define __UNDEF_NAMESPACE_AT_END 1
#endif /* !defiend( __NAMESPACE ) */

#include "sixtracklib/_impl/namespace_begin.h"
#include "sixtracklib/common/mem_pool.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

#include <gtest/gtest.h>

#if defined( __NAMESPACE ) && defined( __UNDEF_NAMESPACE_AT_END )
    #undef __NAMESPACE
    #undef __UNDEF_NAMESPACE_AT_END
#endif /* !defined( __NAMESPACE ) && defined( __UNDEF_NAMESPACE_AT_END ) */

/* ========================================================================== */
/* ====  Test basic usage of MemPool: init and free operations */

TEST( MemPoolTests, InitFreeBasic )
{
    NS( MemPool ) mem_pool;

    /* --------------------------------------------------------------------- */
    /* Test the "happy path": capacity and chunk_size are compatible */

    std::size_t const capacity = std::size_t{64};
    std::size_t const chunk_size = std::size_t{8};

    NS( MemPool_init )( &mem_pool, capacity, chunk_size );

    ASSERT_TRUE( NS( MemPool_get_buffer )( &mem_pool ) != nullptr );
    ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) == capacity );
    ASSERT_TRUE( NS( MemPool_get_chunk_size )( &mem_pool ) == chunk_size );
    ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == std::size_t{0} );

    NS( MemPool_free )( &mem_pool );

    ASSERT_TRUE( NS( MemPool_get_buffer )( &mem_pool ) == nullptr );
    ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) == std::size_t{0} );
    ASSERT_TRUE( NS( MemPool_get_chunk_size )( &mem_pool ) == std::size_t{0} );
    ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == std::size_t{0} );
}

/* ========================================================================== */
/* ====  Test handling of odd-sized but admissible capacities */

TEST( MemPoolTests, InitFreeNonIntegerNumChunks )
{
    NS( MemPool ) mem_pool;

    /* --------------------------------------------------------------------- */
    /* Test how the MemoryPool operates if the capacity is not an integer
     * multiple of the chunk size: */

    std::size_t const capacity = std::size_t{62};
    std::size_t const chunk_size = std::size_t{8};

    static std::size_t const ZERO_SIZE = std::size_t{0};

    NS( MemPool_init )( &mem_pool, capacity, chunk_size );

    ASSERT_TRUE( NS( MemPool_get_buffer )( &mem_pool ) != nullptr );
    ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) >= capacity );

    ASSERT_TRUE( ZERO_SIZE ==
                 ( NS( MemPool_get_capacity )( &mem_pool ) % chunk_size ) );

    ASSERT_TRUE( NS( MemPool_get_chunk_size )( &mem_pool ) == chunk_size );
    ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == std::size_t{0} );

    NS( MemPool_free )( &mem_pool );
}

/* ========================================================================== */
/* ====  Test handling of pathological zero-sized capacities*/

TEST( MemPoolTests, InitFreeZeroCapacityNonZeroChunk )
{
    NS( MemPool ) mem_pool;

    /* --------------------------------------------------------------------- */
    /* Test how the MemoryPool operates if the capacity is not an integer
     * multiple of the chunk size: */

    std::size_t const capacity = std::size_t{0};
    std::size_t const chunk_size = std::size_t{8};

    NS( MemPool_init )( &mem_pool, capacity, chunk_size );

    ASSERT_TRUE( NS( MemPool_get_buffer )( &mem_pool ) != nullptr );
    ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) >= capacity );
    ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) == chunk_size );
    ASSERT_TRUE( NS( MemPool_get_chunk_size )( &mem_pool ) == chunk_size );
    ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == std::size_t{0} );

    NS( MemPool_free )( &mem_pool );
}

/* ========================================================================== */
/* ====  Happy-Path-Testing of adding blocks aligned and non-aligned */

TEST( MemPoolTests, AppendSuccess )
{
    NS( MemPool ) mem_pool;

    std::size_t const chunk_size = std::size_t{8u};
    std::size_t const capacity = 12 * chunk_size;

    NS( MemPool_init )( &mem_pool, capacity, chunk_size );

    /* --------------------------------------------------------------------- */
    std::size_t num_bytes_to_add = std::size_t{2} * chunk_size;
    std::size_t expected_length = num_bytes_to_add;
    uint64_t expected_offset = uint64_t{0};

    NS( AllocResult )
    result = NS( MemPool_append )( &mem_pool, num_bytes_to_add );

    ASSERT_TRUE( NS( AllocResult_valid )( &result ) );
    ASSERT_TRUE( NS( AllocResult_get_pointer )( &result ) ==
                 NS( MemPool_get_buffer )( &mem_pool ) );

    ASSERT_TRUE( NS( AllocResult_get_offset )( &result ) == expected_offset );
    ASSERT_TRUE( NS( AllocResult_get_length )( &result ) == expected_length );

    /* --------------------------------------------------------------------- */

    num_bytes_to_add = ( chunk_size >> 1 );
    expected_length = chunk_size;
    expected_offset += NS( AllocResult_get_length )( &result );

    result = NS( MemPool_append )( &mem_pool, num_bytes_to_add );

    ASSERT_TRUE( NS( AllocResult_valid )( &result ) );
    ASSERT_TRUE( NS( AllocResult_get_offset( &result ) ) == expected_offset );
    ASSERT_TRUE( NS( AllocResult_get_length( &result ) ) == expected_length );
    ASSERT_TRUE( NS( AllocResult_get_pointer( &result ) ) ==
                 ( NS( MemPool_get_buffer )( &mem_pool ) + expected_offset ) );

    /* --------------------------------------------------------------------- */

    std::size_t alignment = chunk_size << 1;
    num_bytes_to_add = chunk_size << 2;
    expected_length = num_bytes_to_add;
    expected_offset += NS( AllocResult_get_length )( &result ) + chunk_size;

    result =
        NS( MemPool_append_aligned )( &mem_pool, num_bytes_to_add, alignment );

    ASSERT_TRUE( NS( AllocResult_valid )( &result ) );
    ASSERT_TRUE( NS( AllocResult_get_offset( &result ) ) == expected_offset );
    ASSERT_TRUE( NS( AllocResult_get_length( &result ) ) == expected_length );
    ASSERT_TRUE( NS( AllocResult_get_pointer( &result ) ) ==
                 ( NS( MemPool_get_buffer )( &mem_pool ) + expected_offset ) );

    /* --------------------------------------------------------------------- */

    alignment = chunk_size;
    num_bytes_to_add = chunk_size << 1;
    expected_length = num_bytes_to_add;
    expected_offset += NS( AllocResult_get_length )( &result );

    result =
        NS( MemPool_append_aligned )( &mem_pool, num_bytes_to_add, alignment );

    ASSERT_TRUE( NS( AllocResult_valid )( &result ) );
    ASSERT_TRUE( NS( AllocResult_get_offset( &result ) ) == expected_offset );
    ASSERT_TRUE( NS( AllocResult_get_length( &result ) ) == expected_length );
    ASSERT_TRUE( NS( AllocResult_get_pointer( &result ) ) ==
                 ( NS( MemPool_get_buffer )( &mem_pool ) + expected_offset ) );

    /* --------------------------------------------------------------------- */

    expected_offset += NS( AllocResult_get_length )( &result );

    ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == expected_offset );
    ASSERT_TRUE( NS( MemPool_get_remaining_bytes )( &mem_pool ) ==
                 capacity - expected_offset );

    /* --------------------------------------------------------------------- */

    NS( MemPool_clear )( &mem_pool );
    ASSERT_TRUE( NS( MemPool_get_buffer )( &mem_pool ) != nullptr );
    ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) == capacity );
    ASSERT_TRUE( NS( MemPool_get_chunk_size )( &mem_pool ) == chunk_size );
    ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == std::size_t{0} );

    /* --------------------------------------------------------------------- */

    NS( MemPool_free )( &mem_pool );
}

/* ========================================================================== */
/* ====  Test the failing of adding blocks with problematic properties        */

TEST( MemPoolTests, AppendFailures )
{
    NS( MemPool ) mem_pool;

    std::size_t const chunk_size = std::size_t{8u};
    std::size_t const capacity = 8 * chunk_size;

    static std::size_t const ZERO_SIZE = std::size_t{0};

    NS( MemPool_init )( &mem_pool, capacity, chunk_size );

    /* --------------------------------------------------------------------- */
    /* Asked to add a block into an empty MemPool that would exceed capacity */

    std::size_t num_bytes_to_add = std::size_t{10} * chunk_size;

    NS( AllocResult )
    result = NS( MemPool_append )( &mem_pool, num_bytes_to_add );

    ASSERT_TRUE( !NS( AllocResult_valid )( &result ) );
    ASSERT_TRUE( NS( AllocResult_get_pointer )( &result ) == nullptr );
    ASSERT_TRUE( NS( AllocResult_get_offset )( &result ) == uint64_t{0} );
    ASSERT_TRUE( NS( AllocResult_get_length )( &result ) == uint64_t{0} );

    ASSERT_TRUE( NS( MemPool_get_buffer )( &mem_pool ) != nullptr );
    ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) == capacity );
    ASSERT_TRUE( NS( MemPool_get_chunk_size )( &mem_pool ) == chunk_size );
    ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == ZERO_SIZE );

    /* --------------------------------------------------------------------- */
    /* Add a block successfully - so we can check whether the MemPool keeps
     * its properties if we insert a block - nonsuccessfully */

    result = NS( MemPool_append )( &mem_pool, chunk_size );
    ASSERT_TRUE( NS( AllocResult_valid )( &result ) );

    std::size_t const current_size = NS( MemPool_get_size )( &mem_pool );

    /* --------------------------------------------------------------------- */
    /* Try to add a block with zero bytes length */

    result = NS( MemPool_append )( &mem_pool, ZERO_SIZE );

    ASSERT_TRUE( !NS( AllocResult_valid )( &result ) );
    ASSERT_TRUE( NS( AllocResult_get_pointer )( &result ) == nullptr );
    ASSERT_TRUE( NS( AllocResult_get_offset )( &result ) == uint64_t{0} );
    ASSERT_TRUE( NS( AllocResult_get_length )( &result ) == uint64_t{0} );

    ASSERT_TRUE( NS( MemPool_get_buffer )( &mem_pool ) != nullptr );
    ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) == capacity );
    ASSERT_TRUE( NS( MemPool_get_chunk_size )( &mem_pool ) == chunk_size );
    ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == current_size );

    /* --------------------------------------------------------------------- */
    /* Try to add a non-zero-length block with a too short alignment -
     * alignments have to be multiples of the MemPool's chunk_size */

    std::size_t const half_alignment = chunk_size >> 1;
    result =
        NS( MemPool_append_aligned )( &mem_pool, chunk_size, half_alignment );

    ASSERT_TRUE( !NS( AllocResult_valid )( &result ) );
    ASSERT_TRUE( NS( AllocResult_get_pointer )( &result ) == nullptr );
    ASSERT_TRUE( NS( AllocResult_get_offset )( &result ) == uint64_t{0} );
    ASSERT_TRUE( NS( AllocResult_get_length )( &result ) == uint64_t{0} );

    ASSERT_TRUE( NS( MemPool_get_buffer )( &mem_pool ) != nullptr );
    ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) == capacity );
    ASSERT_TRUE( NS( MemPool_get_chunk_size )( &mem_pool ) == chunk_size );
    ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == current_size );

    /* --------------------------------------------------------------------- */
    /* Try to add a non-zero-length block with a non-integer-multiple
     * relationship with the MemPool's chunk_size: */

    std::size_t const one_and_a_half_alignment =
        chunk_size + ( chunk_size >> 1 );

    result = NS( MemPool_append_aligned )(
        &mem_pool, chunk_size, one_and_a_half_alignment );

    ASSERT_TRUE( !NS( AllocResult_valid )( &result ) );
    ASSERT_TRUE( NS( AllocResult_get_pointer )( &result ) == nullptr );
    ASSERT_TRUE( NS( AllocResult_get_offset )( &result ) == uint64_t{0} );
    ASSERT_TRUE( NS( AllocResult_get_length )( &result ) == uint64_t{0} );

    ASSERT_TRUE( NS( MemPool_get_buffer )( &mem_pool ) != nullptr );
    ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) == capacity );
    ASSERT_TRUE( NS( MemPool_get_chunk_size )( &mem_pool ) == chunk_size );
    ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == current_size );

    /* --------------------------------------------------------------------- */
    /* Try to add a valid block exceeding the number of remaining bytes      */

    std::size_t remaining_bytes =
        NS( MemPool_get_remaining_bytes )( &mem_pool );

    std::size_t const too_large = remaining_bytes + std::size_t{1};

    result = NS( MemPool_append )( &mem_pool, too_large );

    ASSERT_TRUE( !NS( AllocResult_valid )( &result ) );
    ASSERT_TRUE( NS( AllocResult_get_pointer )( &result ) == nullptr );
    ASSERT_TRUE( NS( AllocResult_get_offset )( &result ) == uint64_t{0} );
    ASSERT_TRUE( NS( AllocResult_get_length )( &result ) == uint64_t{0} );

    ASSERT_TRUE( NS( MemPool_get_buffer )( &mem_pool ) != nullptr );
    ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) == capacity );
    ASSERT_TRUE( NS( MemPool_get_chunk_size )( &mem_pool ) == chunk_size );
    ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == current_size );

    /* --------------------------------------------------------------------- */
    /* Try to add a valid block fitting in memory but exceeding the capacity
     * due to too ambitious alignment requirements: */

    /* Step one: find the "too ambitious alignment" values */

    uint64_t const current_offset =
        NS( MemPool_get_next_begin_offset )( &mem_pool, chunk_size );

    ASSERT_TRUE( current_offset != UINT64_MAX );

    std::size_t alignment = ( chunk_size << 1 );
    std::size_t end_alignment = capacity;

    bool found_alignment_for_test = false;

    for( ; alignment < end_alignment; alignment += chunk_size )
    {
        if( NS( MemPool_get_next_begin_offset )( &mem_pool, alignment ) >
            current_offset )
        {
            result = NS( MemPool_append_aligned(
                &mem_pool, remaining_bytes, alignment ) );

            found_alignment_for_test = true;

            break;
        }
    }

    if( found_alignment_for_test )
    {
        ASSERT_TRUE( !NS( AllocResult_valid )( &result ) );
        ASSERT_TRUE( NS( AllocResult_get_pointer )( &result ) == nullptr );
        ASSERT_TRUE( NS( AllocResult_get_offset )( &result ) == uint64_t{0} );
        ASSERT_TRUE( NS( AllocResult_get_length )( &result ) == uint64_t{0} );

        ASSERT_TRUE( NS( MemPool_get_buffer )( &mem_pool ) != nullptr );
        ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) == capacity );
        ASSERT_TRUE( NS( MemPool_get_chunk_size )( &mem_pool ) == chunk_size );
        ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == current_size );
    }

    /* Verify that non-aligned insert would work, however: */

    result = NS( MemPool_append )( &mem_pool, remaining_bytes );

    ASSERT_TRUE( NS( AllocResult_valid )( &result ) );
    ASSERT_TRUE( NS( AllocResult_get_pointer )( &result ) != nullptr );
    ASSERT_TRUE( NS( AllocResult_get_offset )( &result ) == current_size );
    ASSERT_TRUE( NS( AllocResult_get_length )( &result ) == remaining_bytes );

    ASSERT_TRUE( NS( MemPool_get_buffer )( &mem_pool ) != nullptr );
    ASSERT_TRUE( NS( MemPool_get_capacity )( &mem_pool ) == capacity );
    ASSERT_TRUE( NS( MemPool_get_chunk_size )( &mem_pool ) == chunk_size );
    ASSERT_TRUE( NS( MemPool_get_size )( &mem_pool ) == capacity );
    ASSERT_TRUE( NS( MemPool_get_remaining_bytes )( &mem_pool ) == ZERO_SIZE );

    /* --------------------------------------------------------------------- */

    NS( MemPool_free )( &mem_pool );
}

/* end: sixtracklib/common/tests/test_mem_pool.cpp */
