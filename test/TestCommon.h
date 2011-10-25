#pragma once

#include <gtest/gtest.h>

#include <gsl/gsl_rng.h>

/**
    GTest fixture to initialize gsl_rng.
 */
class GslRngTestF: public ::testing::Test {
protected:
    gsl_rng* rndGen;

    GslRngTestF()
    : rndGen( NULL )
    {
        rndGen = gsl_rng_alloc( gsl_rng_default );
    }

    virtual ~GslRngTestF()
    {
        if ( rndGen ) {
            gsl_rng_free( rndGen );
        }
    }
};

#define TEST_F2(test_fixture, test_case_name, test_name)\
  GTEST_TEST_(test_case_name, test_name, test_fixture, \
              ::testing::internal::GetTypeId<test_fixture>())

extern std::string BIMAPtest_data_path;
