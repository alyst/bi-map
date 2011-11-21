#include <gtest/gtest.h>

#include <symmetric_array2d.h>

#include "TestCommon.h"

class symmetric_array2dTest: public ::testing::Test {
protected:
    symmetric_array2d<size_t> arr;

    symmetric_array2dTest()
    : arr( 5 )
    {
    }

    virtual ~symmetric_array2dTest()
    {
    }
};

TEST_F( symmetric_array2dTest, assignment )
{
    for ( size_t i = 0; i < arr.size(); i++ ) {
        for ( size_t j = 0; j <= i; j++ ) {
            arr( i, j ) = i + j;
        }
    }
    for ( size_t i = 0; i < arr.size(); i++ ) {
        for ( size_t j = 0; j <= i; j++ ) {
            EXPECT_EQ( i+j, arr( i, j ) );
        }
    }
    for ( size_t i = 0; i < arr.size(); i++ ) {
        for ( size_t j = i; j < arr.size(); j++ ) {
            EXPECT_EQ( i+j, arr( i, j ) );
        }
    }
}
