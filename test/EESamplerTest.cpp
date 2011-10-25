#include <BasicTypedefs.h>

#include "EESampler2D.h"

#include "ConsolePTCExecutionMonitor.h"

#include "TestCommon.h"

/*************************************************/
#if 0
TEST( EquiEnergySampler, EnergyLadderStepRatio )
{
    EXPECT_NEAR( TurbineCascadeParams::EnergyLadderStepRatio( 5, 30, 0 ), 2, 0.00001 );
}
#endif


class Particle2dEquiEnergySampler: public ::testing::TestWithParam<Particle2dEquiEnergySamplerTestParam> {
protected:
    gsl_rng* rndGen;

    Particle2dEquiEnergySampler()
    : rndGen( gsl_rng_alloc( gsl_rng_default ) )
    {
    }

    virtual ~Particle2dEquiEnergySampler()
    {
        gsl_rng_free( rndGen );
    }
};

TEST_F2( GslRngTestF, Particle2dEquiEnergySampler, empty_cascade_throw )
{
    typedef TurbineCascadeUnit<Particle2dCollector, DynamicParticle2dFactory, Particle2dInterpolationGenerator> test_turbine_cascade;

    Particle2dEquiEnergySamplerTestParam params( Particle2d( 3, 2 ), 
                                                 Particle2d( 0.1, 0.1 ), 0.5, 
                                                 10000,
                                                 1, 1000, 0.1 );
    DynamicParticle2dFactory    dynFry( rndGen, params.distrParams );
    Particle2dInterpolationGenerator pgen;
    Particle2dCollector pcoll;

    boost::scoped_ptr<test_turbine_cascade> ees( CreateSingleUnitCascade( params.eesParams, rndGen,
                                                        dynFry, &pgen, pcoll ) );
    EXPECT_THROW( ees->iterate(), std::runtime_error );
}

TEST_P( Particle2dEquiEnergySampler, Particle2d_sampling )
{
    typedef TurbineCascadeUnit<Particle2dCollector, DynamicParticle2dFactory, Particle2dInterpolationGenerator> test_turbine_cascade;

    LOG_DEBUG1( "Samples count " << GetParam().samplesCount );
    DynamicParticle2dFactory    dynFry( rndGen, GetParam().distrParams );
    Particle2dInterpolationGenerator pgen;
    Particle2dCollector pcoll( GetParam().samplesCount );

    // sampling
    StdOutPTCExecutionMonitor   mon( 1 );
    boost::scoped_ptr<test_turbine_cascade> ees( CreateSingleUnitCascade( GetParam().eesParams, rndGen,
                                                        dynFry, &pgen, pcoll ) );
    Particle2d  sum;
    Particle2d  sqrsum;
    double      mult = 0;
    ees->init( StaticParticle2d() );
    LOG_DEBUG1( "Generate rate " << GetParam().eesParams.turbineParams.generateRate );
    ees->run();

    for ( size_t ix = 0; ix < pcoll.samples.size(); ix++ ) {
        const Particle2d& cur = pcoll.samples[ ix ];
        sum.first += cur.first;
        sum.second += cur.second;
        sqrsum.first += gsl_pow_2( cur.first );
        sqrsum.second += gsl_pow_2( cur.second );
        mult += cur.first * cur.second;
    }

    // checking
    Particle2d mean;
    mean.first = sum.first / pcoll.samples.size();
    mean.second = sum.second / pcoll.samples.size();
    EXPECT_NEAR( mean.first, GetParam().distrParams.center.first, 0.01 );
    EXPECT_NEAR( mean.second, GetParam().distrParams.center.second, 0.01 );

    Particle2d sd;
    sd.first = sqrt( sqrsum.first / (pcoll.samples.size() - 1) - gsl_pow_2( mean.first ) );
    sd.second = sqrt( sqrsum.second / (pcoll.samples.size() - 1) - gsl_pow_2( mean.second ) );
    EXPECT_NEAR( sd.first, GetParam().distrParams.sigma.first, 0.01 );
    EXPECT_NEAR( sd.second, GetParam().distrParams.sigma.second, 0.01 );

    double ratio = mult / ( pcoll.samples.size() - 1 ) - mean.first * mean.second;
    ratio /= sd.first;
    ratio /= sd.second;
    EXPECT_NEAR( ratio, GetParam().distrParams.ratio, 0.05 );
}

INSTANTIATE_TEST_CASE_P( singleTurbine,
                         Particle2dEquiEnergySampler,
                         ::testing::Values( Particle2dEquiEnergySamplerTestParam( Particle2d( 3, 2 ), 
                                                                                  Particle2d( 0.1, 0.1 ), 0.5, 
                                                                                  10000,
                                                                                  1, 1 ) ) );

INSTANTIATE_TEST_CASE_P( threeTurbines,
    Particle2dEquiEnergySampler,
    ::testing::Values( Particle2dEquiEnergySamplerTestParam( Particle2d( 3, 2 ), 
                                                             Particle2d( 0.1, 0.1 ), 0.5, 
                                                             6000,
                                                             3, 3, 0.05, 0.0 ) ) );

