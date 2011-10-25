#include <gtest/gtest.h>

#include <gsl/gsl_sf_gamma.h>

#include <ObjectSet.h>
#include <BasicTypedefs.h>

#include <SemindexedMap.h>

/*************************************************/
TEST( SemindexedMap, generic )
{
    object_set_t    objs56;
    objs56.insert( 5 );
    objs56.insert( 6 );
    
    typedef EntityIndexing<object_set_t> obj_indexing_type;
    typedef ::SemindexedMap<probe_index_t, object_set_t, double, ::ObjectSetDistance> semimap_type;
    
    obj_indexing_type   objIndexing;
    semimap_type  smap( objIndexing );

    // missing value
    EXPECT_TRUE( smap.find( 1, objs56 ) == smap.end() );
    // inserting
    EXPECT_TRUE( smap.put( 1, objs56, 5 ) );
    EXPECT_EQ( smap.key2index( objs56 )->serial(), 0 );
    // fetching
    EXPECT_EQ( smap.find( 1, objs56 )->value(), 5 );
    EXPECT_EQ( smap.find( 1, 0 )->value(), 5 );
    EXPECT_TRUE( smap.find( 2, objs56 ) == smap.end() );
    EXPECT_TRUE( smap.put( 1, objs56, 6 ) );
    // replacing
    EXPECT_EQ( smap.find( 1, objs56 )->value(), 6 );

    // finding non-exact
    object_set_t    objs5;
    objs5.insert( 5 );
    EXPECT_TRUE( smap.findClosest( 1, objs56, 1 ) == smap.find( 1, objs56 ) );
    EXPECT_TRUE( smap.findClosest( 1, objs5, 1 ) == smap.find( 1, objs56 ) );
    EXPECT_TRUE( smap.findClosest( 1, objs5, 2 ) == smap.find( 1, objs56 ) );
    EXPECT_TRUE( smap.findClosest( 2, objs5, 2 ) == smap.end() );
    EXPECT_TRUE( smap.put( 1, objs5, 7 ) );
    EXPECT_TRUE( smap.findClosest( 1, objs5, 1 ) != smap.find( 1, objs56 ) );

    object_set_t    objs567;
    objs567.insert( 5 );
    objs567.insert( 6 );
    objs567.insert( 7 );
    EXPECT_TRUE( smap.findClosest( 1, objs567, 1 ) == smap.find( 1, objs56 ) );
    EXPECT_TRUE( smap.put( 1, objs567, 9 ) );

    object_set_t    objs67;
    objs67.insert( 6 );
    objs67.insert( 7 );
    EXPECT_TRUE( smap.findClosest( 1, objs67, 1 ) == smap.find( 1, objs567 ) );

    object_set_t    objs4;
    objs4.insert( 4 );
    EXPECT_TRUE( smap.findClosest( 1, objs4, 1 ) == smap.end() );
    EXPECT_TRUE( smap.findClosest( 1, objs4, 2 ) == smap.find( 1, objs5 ) );
}
