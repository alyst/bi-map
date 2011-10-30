#include "../BasicTypedefs.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "../ConsolePTCExecutionMonitor.h"

#include "../BIMAPResultsSerialize.h"
#include "../ParametersReader.h"
#include "../OPADataImportCSV.h"
#include "BIMAP-sampler-mpi.h"

boost::optional<BIMAPWalk> MPI_BIMAPSampler_run(
    const BIMAPSamplerHelper&      helper,
    mpi::communicator&              comm,
    const TurbineCascadeParams&     eeCascadeParams,
    const ChessboardBiclustering&          iniClus,            /** initial clustering */
    ChessboardBiclusteringsIndexing*       ccIndexing,
    const BIMAPSampleCollectorParams*      collectorParams,
    const TurbineCascadeExecutionMonitor*   pMon
){
    typedef TurbineCascadeUnit<BIMAPSampleCollector,
                               DynamicChessboardBiclusteringFactory, ChessboardBiclusteringCrossoverGenerator,
                               MPIUnitCommunicator<StaticChessboardBiclustering> > equi_energy_sampler_type;

    boost::scoped_ptr<BIMAPSampleCollector> collector;
    if ( ccIndexing != NULL && collectorParams != NULL ) {
        collector.reset( new BIMAPSampleCollector( *ccIndexing, *collectorParams ) );
    }

    LOG_INFO( comm.rank() << ": initializing " 
                          << eeCascadeParams.turbinesCount << " node(s), "
                          << eeCascadeParams.levelsCount << " level(s) cascade..." );

    LOG_DEBUG1( "MPI Turbine Cascade Creation for communicator of size " << comm.size() );

    int collectorRank = MPIGetCollectorRank( comm, eeCascadeParams, collector.get() );
    if ( collectorRank == - 1 ) THROW_RUNTIME_ERROR( "No collecting units found" );

    tcascade_structure_type cascadeStructure;
    if ( collectorRank == comm.rank() ) {
        cascadeStructure = eeCascadeParams.turbinesPerLevel.size() > 0
            ? createTurbineCascadeStructure( eeCascadeParams.turbinesPerLevel, comm.size() )
            : createTurbineCascadeStructure( eeCascadeParams.levelsCount, eeCascadeParams.turbinesCount,
                                             comm.size() );
    }
    mpi::broadcast( comm, cascadeStructure, collectorRank );


    typedef MPIUnitCommunicator<StaticChessboardBiclustering> unit_communicator_type;
    unit_communicator_type unitComm( comm, collectorRank );

    equi_energy_sampler_type eeSampler(
                comm.rank(), eeCascadeParams, cascadeStructure,
                helper.rndNumGen, helper.dynamicCrossCluFactory,
                &helper.crossoverGenerator, collector.get(), &unitComm );
    {
        ChessboardBiclustering iniClu( iniClus );
        iniClu.check();
        eeSampler.init( StaticChessboardBiclustering( iniClu ) );
    }
    LOG_INFO( comm.rank() << ": running sampling..." );
    eeSampler.run();
    LOG_INFO( comm.rank() << ": finished sampling in "
              << boost::format( "%0.f" ) % eeSampler.elapsed() << " second(s)" );
    return ( collector ? collector->walk() : boost::optional<BIMAPWalk>() );
}

template<class Object>
void broadcast_statically_tracked( boost::mpi::communicator& comm, const char* name, boost::scoped_ptr<Object>& object, int root )
{
    LOG_INFO( "#" << comm.rank() << ": " 
                << ( comm.rank() == root ? "broadcasting " : "receiving " ) << name << "..." );
    const Object* pObj = object.get();
    if ( ( pObj == NULL ) && ( comm.rank() == root ) ) {
        THROW_RUNTIME_ERROR( "Root does not have data for broadcasting" );
    }
    else if ( ( pObj != NULL ) && ( comm.rank() != root ) ) {
        THROW_RUNTIME_ERROR( "Non-root has non-empty data before broadcasting" );
    }
    statically_tracked<Object> pstObj( name, pObj );
    boost::mpi::broadcast( comm, pstObj, root );
    LOG_INFO( "#" << comm.rank() << ": " 
                << ( comm.rank() == root ? "broadcast " : "receive " ) << name << " complete" );
    if ( comm.rank() != root ) {
        object.reset( const_cast<Object*>( pObj ) );
    }
}

int main( int argc, char* argv[] )
{
    boost::mpi::environment env( argc, argv, true );
    boost::mpi::communicator world;

    LOG_INFO( "#" << world.rank() << ": MPI environment initialized, "
              << world.size() << " process(es) in total..." );
    ChessboardBiclusteringHyperPriors  hyperpriors;
    hyperpriors.signalHyperprior.meanVarScale = 2;
    GibbsSamplerParams          gibbsParams;
    TurbineCascadeParams        cascadeParams;
    cascadeParams.turbinesCount = world.size();
    //cascadeParams.maxTemperature = 5;

    boost::scoped_ptr<CellSignalParams>         signalParams;
    boost::scoped_ptr<PrecomputedDataParams>    precomputedDataParams;
    boost::scoped_ptr<OPAData> data;
    boost::scoped_ptr<ChessboardBiclusteringPriors> priors;
    boost::scoped_ptr<PrecomputedData> precomputed;
    boost::scoped_ptr<BIMAPSampleCollectorParams>    collectorParams;
    boost::scoped_ptr<BIMAPIOParams>  ioParams;
    bool params_res = false;
    if ( world.rank() == 0 ) {
        priors.reset( new ChessboardBiclusteringPriors() );
        priors->probeClustering.concentration = 0.01;
        priors->cellEnablementProb = 0.5;

        signalParams.reset( new CellSignalParams() );
        signalParams->scShape = 0.3;
        signalParams->sequenceLengthFactor = 1.0;

        precomputedDataParams.reset( new PrecomputedDataParams() );

        ioParams.reset( new BIMAPIOParams() );
        collectorParams.reset( new BIMAPSampleCollectorParams() );

        params_res = BIMAPParamsRead( argc, argv,
                          hyperpriors, gibbsParams, cascadeParams,
                          *signalParams, *precomputedDataParams, *priors,
                          *collectorParams, *ioParams );

        if ( params_res ) {
        if ( !ioParams->dataFilename.empty() ) {
            boost::filesystem::path data_file_path( ioParams->dataFilename );
            LOG_DEBUG1( "#" << world.rank() << ": loading data from " << data_file_path << "..." );
            if ( !boost::filesystem::exists( data_file_path ) ) {
                THROW_RUNTIME_ERROR( "OPAData file not found: " << data_file_path );
            }
            LOG_INFO( "Loading data..." );
            data.reset( new OPAData( OPAData::load( data_file_path.string().c_str() ) ) );
        }
        else if ( !ioParams->proteinsFilename.empty() ) {
            boost::filesystem::path proteins_file_path( ioParams->proteinsFilename );
            if ( !boost::filesystem::exists( proteins_file_path ) ) {
                THROW_RUNTIME_ERROR( "Proteins file not found: " << proteins_file_path );
            }
            boost::filesystem::path exp_design_file_path( ioParams->expDesignFilename );
            if ( !boost::filesystem::exists( exp_design_file_path ) ) {
                THROW_RUNTIME_ERROR( "Experimental design file not found: " << exp_design_file_path );
            }
            boost::filesystem::path measurements_file_path( ioParams->measurementsFilename );
            if ( !boost::filesystem::exists( measurements_file_path ) ) {
                THROW_RUNTIME_ERROR( "Measurements file not found: " << measurements_file_path );
            }
            data.reset( new OPAData( OPADataImportCSV( proteins_file_path.string().c_str(),
                                    exp_design_file_path.string().c_str(),
                                    measurements_file_path.string().c_str(),
                                    ioParams->mapBaitsToObjects,
                                    ioParams->csvColumnSeparator ) ) );
        } else {
            THROW_RUNTIME_ERROR( "No input data" );
        }

        LOG_INFO( "Precomputing..." );
        precomputed.reset( new PrecomputedData( *data, *precomputedDataParams, *signalParams ) );
        }
    }
    LOG_INFO( "#" << world.rank() << ": " 
              << ( world.rank() == 0 ? "broadcasting" : "receiving" )
              << " input data" );
    mpi::broadcast( world, params_res, 0 );
    if ( !params_res ) return ( 0 ); // --help option, no computation

    mpi::broadcast( world, hyperpriors, 0 );
    mpi::broadcast( world, gibbsParams, 0 );
    mpi::broadcast( world, cascadeParams, 0 );
    broadcast_statically_tracked( world, "data", data, 0 );
    broadcast_statically_tracked( world, "signalParams", signalParams, 0 );
    broadcast_statically_tracked( world, "priors", priors, 0 );
    broadcast_statically_tracked( world, "precomputedDataParams", precomputedDataParams, 0 );
    broadcast_statically_tracked( world, "precomputed", precomputed, 0 );
    LOG_DEBUG1_IF( world.rank() != 0, "#" << world.rank() << ": received OPAData (" 
                    << data->objectsCount() << " objects, " << data->probesCount() << " probes)" );
    LOG_INFO( "#" << world.rank() << ": initializing sampler..." );
    BIMAPSamplerHelper         helper( *precomputed,
                                        hyperpriors, *priors, gibbsParams,
                                        world.rank() * 1981 );

    LOG_DEBUG1( "#" << world.rank() << ": setting initial clustering..." );
    ChessboardBiclustering iniClus;
    if ( world.rank() == 0 ) iniClus = helper.trivialClustering();
    mpi::broadcast( world, iniClus, 0 );

    StdOutPTCExecutionMonitor   mon( 1 );
    boost::scoped_ptr<ChessboardBiclusteringsIndexing> ccIndexing;
    if ( world.rank() == 0 ) ccIndexing.reset( new ChessboardBiclusteringsIndexing() );

    boost::optional<BIMAPWalk> res = MPI_BIMAPSampler_run( helper, world, cascadeParams, iniClus, 
                                                             ccIndexing.get(), collectorParams.get(), &mon );

    if ( res ) {
        LOG_INFO( "Process #" << world.rank()
                  << " finished with " << res->stepsCount() << " samples collected" );
        BOOST_ASSERT( res->check() );
        if ( !ioParams->outputFilename.empty() ) {
            LOG_INFO( "Saving walk to file " << ioParams->outputFilename << "..." );
            size_t removedSteps = res->filterSteps( ioParams->minCrossClusRefCount,
                                                    ioParams->minObjectsPtnRefCount,
                                                    ioParams->minProbesPtnRefCount );
            LOG_INFO_IF( removedSteps > 0,
                         "Removed " << removedSteps << " samples"
                         << " of rarely encountered chessboard biclusterings" );
            ccIndexing->remove_unreferenced();
            BIMAPResultsSave( ioParams->outputFilename.c_str(),
                               data->objectsLabelMap(), data->probesLabelMap(),
                               *res, NULL );
            LOG_INFO( "Saving done\n" );
        } else {
            LOG_INFO( "Walk was not saved: no walk filename specified\n" );
        }
    }
    else {
        LOG_INFO( "Process #" << world.rank() << " finished (not collector)" );
    }
}
