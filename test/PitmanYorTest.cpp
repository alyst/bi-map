#include <gsl/gsl_sf_zeta.h>

#include <math/PitmanYorProcess.h>

#include "TestCommon.h"

/*************************************************/

#define PY_EPS ( 20 * std::numeric_limits<log_prob_t>::epsilon() )

TEST_F2( GslRngTestF, PitmanYor, random_sample )
{
    PitmanYorProcess    py( 0.1, 0 );

    PitmanYorSample sample = py.random( rndGen, 100 );
    EXPECT_EQ( sample.samplesCount(), 100 );
    EXPECT_GT( sample.clustersCount(), 0 );
}

struct PitmanYorRandomTestParam {
    PitmanYorProcess    py;
    const size_t        samplesCount;
    const size_t        iterationsCount;

    PitmanYorRandomTestParam( const PitmanYorProcess& py, const size_t samplesCount, const size_t iterationsCount )
    : py( py )
    , samplesCount( samplesCount )
    , iterationsCount( iterationsCount )
    {
    }
};

struct PitmanYorRandomTest: public ::testing::TestWithParam<PitmanYorRandomTestParam> {
    gsl_rng*            rndGen;

    PitmanYorRandomTest()
    : rndGen( gsl_rng_alloc( gsl_rng_default ) )
    {
    }

    virtual ~PitmanYorRandomTest()
    {
        gsl_rng_free( rndGen );
    }
};

TEST_P( PitmanYorRandomTest, mean_expected_clusters )
{
    const PitmanYorRandomTestParam& param = GetParam();

    size_t  totalClusters = 0;
    for ( size_t i = 0; i < param.iterationsCount; i++ ) {
        PitmanYorSample samplei = param.py.random( rndGen, param.samplesCount );
        totalClusters += samplei.clustersCount();
    }
    EXPECT_NEAR( (double)totalClusters / param.iterationsCount, param.py.expectedClustersCount( param.samplesCount ), 
                 0.0125 * param.py.expectedClustersCount( param.samplesCount ) );
}

INSTANTIATE_TEST_CASE_P( mean_DP_PYP_test,
                         PitmanYorRandomTest,
                         ::testing::Values( PitmanYorRandomTestParam( PitmanYorProcess( 0.1, 0 ), 100, 5000 ),
                                            PitmanYorRandomTestParam( PitmanYorProcess( 0.1, 0 ), 1000, 5000 ) ,
                                            PitmanYorRandomTestParam( PitmanYorProcess( 0.1, 0.2 ), 100, 5000 ), 
                                            PitmanYorRandomTestParam( PitmanYorProcess( 0.1, 0.2 ), 100, 8000 ),
                                            PitmanYorRandomTestParam( PitmanYorProcess( 0.1, 0.5 ), 1000, 5000 ) ,
                                            PitmanYorRandomTestParam( PitmanYorProcess( 0.2, 0.2 ), 100, 5000 ) 
                         ) 
);

struct PitmanYorLnPTestParam {
    PitmanYorProcess    py;

    PitmanYorLnPTestParam( const PitmanYorProcess& py )
    : py( py )
    {
    }
};

struct PitmanYorLnPTest: public ::testing::TestWithParam<PitmanYorLnPTestParam> {
    gsl_rng*            rndGen;

    PitmanYorLnPTest()
    : rndGen( gsl_rng_alloc( gsl_rng_default ) )
    {
    }

    virtual ~PitmanYorLnPTest()
    {
        gsl_rng_free( rndGen );
    }
};

TEST_P( PitmanYorLnPTest, iterative_priors )
{
    const PitmanYorLnPTestParam& param = GetParam();

    double lnPIterative1 = log( param.py.clusterAssignmentPrior( 0, 0, 0 )
                         * param.py.clusterAssignmentPrior( 0, 1, 1 )
                         * param.py.clusterAssignmentPrior( 1, 2, 2 )
                         * param.py.clusterAssignmentPrior( 1, 2, 3 ) );
    double lnPIterative2 = log( param.py.clusterAssignmentPrior( 0, 0, 0 )
                         * param.py.clusterAssignmentPrior( 1, 1, 1 )
                         * param.py.clusterAssignmentPrior( 0, 1, 2 )
                         * param.py.clusterAssignmentPrior( 1, 2, 3 ) );
    EXPECT_NEAR( lnPIterative1, lnPIterative2, PY_EPS );

    std::vector<size_t> clusterSizes;
    clusterSizes.push_back( 2 );
    clusterSizes.push_back( 2 );
    double lnPInstant = param.py.lnP( clusterSizes );
    EXPECT_NEAR( lnPIterative1, lnPInstant, PY_EPS );
    EXPECT_NEAR( lnPIterative1, lnPInstant, PY_EPS );
}

TEST_P( PitmanYorLnPTest, element_reassignment )
{
    const PitmanYorLnPTestParam& param = GetParam();

    std::vector<size_t> clusterSizesBefore;
    clusterSizesBefore.push_back( 2 );
    clusterSizesBefore.push_back( 2 );

    std::vector<size_t> clusterSizesAfter;
    clusterSizesAfter.push_back( 3 );
    clusterSizesAfter.push_back( 1 );

    double lnPBefore = param.py.lnP( clusterSizesBefore );
    double lnPAfter = param.py.lnP( clusterSizesAfter );

    double lnAssBefore = log( param.py.clusterAssignmentPrior( 1, 2, 3 ) );
    double lnAssAfter = log( param.py.clusterAssignmentPrior( 2, 2, 3 ) );

    EXPECT_NEAR( lnPAfter - lnPBefore, lnAssAfter - lnAssBefore, PY_EPS );
}

INSTANTIATE_TEST_CASE_P( iterative_priors_test,
                         PitmanYorLnPTest,
                         ::testing::Values( PitmanYorLnPTestParam( PitmanYorProcess( 0.1, 0 ) ),
                                            PitmanYorLnPTestParam( PitmanYorProcess( 0.1, 0.2 ) ), 
                                            PitmanYorLnPTestParam( PitmanYorProcess( 0.1, 0.5 ) ) ,
                                            PitmanYorLnPTestParam( PitmanYorProcess( 0.2, 0.2 ) ) 
                         ) 
);
