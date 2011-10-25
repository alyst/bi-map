#include <BasicTypedefs.h>

#include "../EESampler2D.h"

#include <mpisampler/MPIUnitCommunicator.h>

#include <ConsolePTCExecutionMonitor.h>

#include "../TestCommon.h"

BOOST_IS_MPI_DATATYPE( Particle2d )
BOOST_MPI_SERIALIZATION_TRAITS( Particle2d )

class MPI_Particle2dEquiEnergySampler: public ::testing::TestWithParam<Particle2dEquiEnergySamplerTestParam> {
protected:
    gsl_rng* rndGen;

    MPI_Particle2dEquiEnergySampler()
    : rndGen( gsl_rng_alloc( gsl_rng_default ) )
    {
    }

    virtual ~MPI_Particle2dEquiEnergySampler()
    {
        gsl_rng_free( rndGen );
    }
};

int main(int argc, char* argv[])
{
    Particle2dEquiEnergySamplerTestParam testParams( Particle2d( 3, 2 ), 
                                                     Particle2d( 0.1, 0.1 ), 0.5, 
                                                     100,
                                                     2, 2 );
    typedef MPIUnitCommunicator<StaticParticle2d> particle2d_unit_communicator;
    typedef TurbineCascadeUnit<Particle2dCollector, DynamicParticle2dFactory, Particle2dInterpolationGenerator,
                               particle2d_unit_communicator> particle2d_cascade_unit;
    typedef particle2d_cascade_unit test_turbine_cascade;

    mpi::environment env( argc, argv, true );

    gsl_rng* rndGen = gsl_rng_alloc( gsl_rng_default );

    particle2d_unit_communicator::mpi_communicator_type world;

    LOG_DEBUG1( "#" << world.rank() << ": samples count " << testParams.samplesCount );
    LOG_DEBUG1( "#" << world.rank() << ": generate rate " << testParams.eesParams.turbineParams.generateRate );

    Particle2dEval              eval( testParams.distrParams );
    DynamicParticle2dFactory    dynFry( rndGen, testParams.distrParams );
    Particle2dInterpolationGenerator pgen;
    Particle2dCollector pcoll( testParams.samplesCount );

    // sampling
    StdOutPTCExecutionMonitor   mon( 1 );
    tcascade_structure_type cascadeStructure = testParams.eesParams.turbinesPerLevel.size() > 0
        ? createTurbineCascadeStructure( testParams.eesParams.turbinesPerLevel, world.size() )
        : createTurbineCascadeStructure( testParams.eesParams.levelsCount, testParams.eesParams.turbinesCount,
                                         world.size() );

    int collectorRank = MPIGetCollectorRank( world, testParams.eesParams,
                                             world.rank() == 0 ? &pcoll : NULL );
    if ( collectorRank == - 1 ) THROW_RUNTIME_ERROR( "No collecting units found" );

    particle2d_unit_communicator unitComm( world, collectorRank );

    test_turbine_cascade ees( world.rank(), testParams.eesParams, cascadeStructure,
                rndGen, dynFry,
                &pgen, world.rank() == collectorRank ? &pcoll : NULL, &unitComm );
    ees.init( StaticParticle2d() );
    ees.run();

    LOG_INFO( "Process #" << world.rank() << " finished with " << pcoll.samples.size() << " samples collected" );

    if ( pcoll.samples.size() > 0 ) {
        Particle2d  sum;
        Particle2d  sqrsum;
        double      mult = 0;

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
        //EXPECT_NEAR( mean.first, GetParam().distrParams.center.first, 0.01 );
        //EXPECT_NEAR( mean.second, GetParam().distrParams.center.second, 0.01 );

        Particle2d sd;
        sd.first = sqrt( sqrsum.first / (pcoll.samples.size() - 1) - gsl_pow_2( mean.first ) );
        sd.second = sqrt( sqrsum.second / (pcoll.samples.size() - 1) - gsl_pow_2( mean.second ) );
        //EXPECT_NEAR( sd.first, GetParam().distrParams.sigma.first, 0.01 );
        //EXPECT_NEAR( sd.second, GetParam().distrParams.sigma.second, 0.01 );

        double ratio = mult / ( pcoll.samples.size() - 1 ) - mean.first * mean.second;
        ratio /= sd.first;
        ratio /= sd.second;
        //EXPECT_NEAR( ratio, GetParam().distrParams.ratio, 0.05 );
    }

    gsl_rng_free( rndGen );
}
