#include <gtest/gtest.h>

#include <OrderedDistanceMatrix.h>

#include "TestCommon.h"

class OrderedDistanceMatrixTest: public ::testing::Test {
protected:
    typedef size_t dist_type;
    typedef size_t ix_type;
    typedef std::set<ix_type> ix_set_type;
    typedef OrderedDistanceMatrix<size_t, dist_type> odm_type;
    typedef odm_type::dist_to_elm_type dte_type;
    symmetric_array2d<dist_type> dist;
    odm_type arr;

    struct ElmContext
    {
        const ix_set_type& elms;

        ElmContext( const ix_set_type& elms )
        : elms( elms )
        {}

        bool operator()( ix_type elm ) const {
            return ( elms.count( elm ) > 0 );
        }
    };

    OrderedDistanceMatrixTest()
    : dist( 5 )
    {
        dist( 0, 0 ) = 0;
        dist( 0, 1 ) = 1;
        dist( 0, 2 ) = 2;
        dist( 0, 3 ) = 3;
        dist( 0, 4 ) = 4;
        dist( 1, 1 ) = 0;
        dist( 1, 2 ) = 3;
        dist( 1, 3 ) = 2;
        dist( 1, 4 ) = 4;
        dist( 2, 2 ) = 0;
        dist( 2, 3 ) = 3;
        dist( 2, 4 ) = 1;
        dist( 3, 3 ) = 0;
        dist( 3, 4 ) = 1;
        dist( 4, 4 ) = 0;
        arr = odm_type( dist );
    }

    virtual ~OrderedDistanceMatrixTest()
    {
    }
};

TEST_F( OrderedDistanceMatrixTest, calculation )
{
    ix_set_type e12;
    e12.insert( 1 );
    e12.insert( 2 );

    EXPECT_EQ( dte_type( 3, 3 ), arr.rankedDistance( 0, ElmContext( e12 ), false, 0, true ) );
    EXPECT_EQ( dte_type( 4, 4 ), arr.rankedDistance( 0, ElmContext( e12 ), false, 1, true ) );
    EXPECT_EQ( dte_type( 4, 4 ), arr.rankedDistance( 0, ElmContext( e12 ), false, 0, false ) );
    EXPECT_EQ( dte_type( 3, 3 ), arr.rankedDistance( 0, ElmContext( e12 ), false, 1, false ) );

    EXPECT_EQ( dte_type( 1, 1 ), arr.rankedDistance( 0, ElmContext( e12 ), true, 0, true ) );
    EXPECT_EQ( dte_type( 2, 2 ), arr.rankedDistance( 0, ElmContext( e12 ), true, 1, true ) );
    EXPECT_EQ( dte_type( 2, 2 ), arr.rankedDistance( 0, ElmContext( e12 ), true, 0, false ) );
    EXPECT_EQ( dte_type( 1, 1 ), arr.rankedDistance( 0, ElmContext( e12 ), true, 1, false ) );

    ix_set_type e0;
    e0.insert( 0 );

    EXPECT_EQ( dte_type( 0, 1 ), arr.rankedDistance( 1, ElmContext( e12 ), false, 0, true ) );
    EXPECT_EQ( dte_type( 3, 2 ), arr.rankedDistance( 1, ElmContext( e12 ), false, 1, true ) );
    EXPECT_EQ( dte_type( 2, 3 ), arr.rankedDistance( 1, ElmContext( e0 ), false, 1, true ) );
    EXPECT_EQ( dte_type( 4, 4 ), arr.rankedDistance( 1, ElmContext( e0 ), false, 0, false ) );
    EXPECT_EQ( dte_type( 2, 3 ), arr.rankedDistance( 1, ElmContext( e0 ), false, 1, false ) );
}
