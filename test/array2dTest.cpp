#include <gtest/gtest.h>

#include <array2d.h>

#include "TestCommon.h"

class array2dTest: public ::testing::Test {
protected:
    array2d<size_t> arr;

    array2dTest()
    : arr( 3, 5 )
    {
    }

    virtual ~array2dTest()
    {
    }
};

TEST_F( array2dTest, insert_row )
{
    for ( size_t i = 0; i < arr.size1(); i++ ) {
        for ( size_t j = 0; j < arr.size2(); j++ ) {
            arr( i, j ) = i;
        }
    }
    array2d<size_t> arr2 = arr;
    arr2.insert1( 1, 98 );
    EXPECT_EQ( 4, arr2.size1() );
    for ( size_t i = 0; i < arr2.size1(); i++ ) {
        for ( size_t j = 0; j < arr2.size2(); j++ ) {
            EXPECT_EQ( i < 1 ? i : ( i > 1 ? i-1 : 98 ), arr2( i, j ) ); 
        }
    }
    arr2 = arr;
    arr2.insert1( 3, 98 );
    EXPECT_EQ( arr2.size1(), 4 );
    for ( size_t i = 0; i < arr2.size1(); i++ ) {
        for ( size_t j = 0; j < arr2.size2(); j++ ) {
            EXPECT_EQ( i < 3 ? i : 98, arr2( i, j ) ); 
        }
    }
}

TEST_F( array2dTest, insert_col )
{
    for ( size_t i = 0; i < arr.size1(); i++ ) {
        for ( size_t j = 0; j < arr.size2(); j++ ) {
            arr( i, j ) = j;
        }
    }
    array2d<size_t> arr2 = arr;
    arr2.insert2( 1, 98 );
    EXPECT_EQ( 6, arr2.size2() );
    for ( size_t i = 0; i < arr2.size1(); i++ ) {
        for ( size_t j = 0; j < arr2.size2(); j++ ) {
            EXPECT_EQ( j < 1 ? j : ( j > 1 ? j-1 : 98 ), arr2( i, j ) ); 
        }
    }
    arr2 = arr;
    arr2.insert2( 5, 98 );
    EXPECT_EQ( 6, arr2.size2() );
    for ( size_t i = 0; i < arr2.size1(); i++ ) {
        for ( size_t j = 0; j < arr2.size2(); j++ ) {
            EXPECT_EQ( j < 5 ? j : 98, arr2( i, j ) ); 
        }
    }
}

TEST_F( array2dTest, remove_row )
{
    for ( size_t i = 0; i < arr.size1(); i++ ) {
        for ( size_t j = 0; j < arr.size2(); j++ ) {
            arr( i, j ) = i;
        }
    }
    array2d<size_t> arr2 = arr;
    arr2.remove1( 0, 98 );
    EXPECT_EQ( 2, arr2.size1() );
    for ( size_t i = 0; i < arr2.size1(); i++ ) {
        for ( size_t j = 0; j < arr2.size2(); j++ ) {
            EXPECT_EQ( i+1, arr2( i, j ) ); 
        }
    }
    arr2 = arr;
    arr2.remove1( 1, 98 );
    EXPECT_EQ( 2, arr2.size1() );
    for ( size_t i = 0; i < arr2.size1(); i++ ) {
        for ( size_t j = 0; j < arr2.size2(); j++ ) {
            EXPECT_EQ( i < 1 ? i : i+1, arr2( i, j ) ); 
        }
    }
    arr2 = arr;
    arr2.remove1( 2, 98 );
    EXPECT_EQ( 2, arr2.size1() );
    for ( size_t i = 0; i < arr2.size1(); i++ ) {
        for ( size_t j = 0; j < arr2.size2(); j++ ) {
            EXPECT_EQ( i, arr2( i, j ) ); 
        }
    }
}

TEST_F( array2dTest, remove_col )
{
    for ( size_t i = 0; i < arr.size1(); i++ ) {
        for ( size_t j = 0; j < arr.size2(); j++ ) {
            arr( i, j ) = j;
        }
    }
    array2d<size_t> arr2 = arr;
    arr2.remove2( 0, 98 );
    EXPECT_EQ( 4, arr2.size2() );
    for ( size_t i = 0; i < arr2.size1(); i++ ) {
        for ( size_t j = 0; j < arr2.size2(); j++ ) {
            EXPECT_EQ( j+1, arr2( i, j ) ); 
        }
    }
    arr2 = arr;
    arr2.remove2( 1, 98 );
    EXPECT_EQ( 4, arr2.size2() );
    for ( size_t i = 0; i < arr2.size1(); i++ ) {
        for ( size_t j = 0; j < arr2.size2(); j++ ) {
            EXPECT_EQ( j < 1 ? j : j+1, arr2( i, j ) ); 
        }
    }
    arr2 = arr;
    arr2.remove2( 4, 98 );
    EXPECT_EQ( 4, arr2.size2() );
    for ( size_t i = 0; i < arr2.size1(); i++ ) {
        for ( size_t j = 0; j < arr2.size2(); j++ ) {
            EXPECT_EQ( j, arr2( i, j ) ); 
        }
    }
}
