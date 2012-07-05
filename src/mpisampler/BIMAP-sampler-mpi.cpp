#include "cemm/bimap/BasicTypedefs.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#if defined(USE_BOOST_LOG)
#include <boost/asio/ip/host_name.hpp>
#include <boost/log/utility/init/common_attributes.hpp>
#include <boost/log/utility/init/to_console.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
#include <boost/log/filters.hpp>
#include <boost/log/formatters.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/severity_channel_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/empty_deleter.hpp>
#include <boost/log/attributes/constant.hpp>
#include <boost/log/attributes/timer.hpp>
#endif

#include <cemm/eesampler/ConsolePTCExecutionMonitor.h>
#include <cemm/eesampler/MPITrace.h>

#include "cemm/bimap/BIMAPResultsSerialize.h"
#include "cemm/bimap/OPADataImportCSV.h"

#include "../ParametersReader.h"
#include "BIMAP-sampler-mpi.h"

namespace log_src = boost::log::sources;
namespace log_fmt = boost::log::formatters;
namespace log_flt = boost::log::filters;
namespace log_sinks = boost::log::sinks;
namespace log_attrs = boost::log::attributes;
namespace log_keywords = boost::log::keywords;

namespace cemm { namespace bimap {

BIMAPWalk* MPI_BIMAPSampler_run(
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
                               MPIUnitCommunicator<ChessboardBiclusteringEnergyEval> > equi_energy_sampler_type;

    boost::scoped_ptr<BIMAPSampleCollector> collector;
    if ( ccIndexing != NULL && collectorParams != NULL ) {
        collector.reset( new BIMAPSampleCollector( *ccIndexing, *collectorParams ) );
    }

    LOG_INFO( "Initializing " 
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
    {
        VT_TRACER( "Broadcast cascade structure" );
        mpi::broadcast( comm, cascadeStructure, collectorRank );
    }


    typedef MPIUnitCommunicator<ChessboardBiclusteringEnergyEval> unit_communicator_type;
    VT_USER_START( "MPIUnitCommunicator construct" );
   	unit_communicator_type unitComm( comm, collectorRank );
    VT_USER_END( "MPIUnitCommunicator construct" );

    VT_USER_START( "EE sampler construct" );
    equi_energy_sampler_type eeSampler(
                comm.rank(), eeCascadeParams, cascadeStructure,
                helper.rndNumGen, helper.dynamicCrossCluFactory,
                &helper.crossoverGenerator, collector.get(), &unitComm );
    VT_USER_END( "EE sampler construct" );
    {
        ChessboardBiclustering iniClu( iniClus );
        iniClu.check();
        eeSampler.init( StaticChessboardBiclustering( iniClu ) );
    }
    LOG_INFO( "Running sampling..." );
    VT_USER_START( "EE sampler run" );
    eeSampler.run();
    VT_USER_END( "EE sampler run" );
    LOG_INFO( "Finished sampling in "
              << boost::format( "%0.f" ) % eeSampler.elapsed() << " second(s)" );
    return ( collector ? new BIMAPWalk( collector->walk() )
                       : nullptr );
}

template<class Object>
void broadcast_statically_tracked( boost::mpi::communicator& comm, const char* name, boost::scoped_ptr<Object>& object, int root )
{
#if defined(VTRACE)
	std::ostringstream msg;
	msg << "broadcast_statically_tracked('" << name << "')";
    VT_TRACER( msg.str().c_str() );
#endif
    LOG_INFO( ( comm.rank() == root ? "broadcasting " : "receiving " ) << name << "..." );
    const Object* pObj = object.get();
    if ( ( pObj == NULL ) && ( comm.rank() == root ) ) {
        THROW_RUNTIME_ERROR( "Root does not have data for broadcasting" );
    }
    else if ( ( pObj != NULL ) && ( comm.rank() != root ) ) {
        THROW_RUNTIME_ERROR( "Non-root has non-empty data before broadcasting" );
    }
    statically_tracked<Object> pstObj( name, pObj );
    boost::mpi::broadcast( comm, pstObj, root );
    LOG_INFO( ( comm.rank() == root ? "broadcast " : "receive " ) << name << " complete" );
    if ( comm.rank() != root ) {
        object.reset( const_cast<Object*>( pObj ) );
    }
}

} }

int main( int argc, char* argv[] )
{
    using namespace cemm::bimap;
    boost::scoped_ptr<ChessboardBiclusteringsIndexing> ccIndexing;
    boost::scoped_ptr<BIMAPIOParams>  ioParams;
    boost::scoped_ptr<OPAData> data;
    boost::scoped_ptr<BIMAPWalk> res;
{
    std::cout << "Starting program" << std::endl;
    boost::mpi::environment env( argc, argv, true );
    boost::mpi::communicator world;

#if defined(USE_BOOST_LOG)
    typedef log_attrs::constant<std::string> hostname_attr_t;
    typedef log_attrs::constant<int>  mpi_rank_attr_t;
    typedef log_attrs::timer         timer_attr_t;

    boost::shared_ptr< boost::log::core > log_core = boost::log::core::get();

    log_core->add_global_attribute( "hostname",
                                    hostname_attr_t( boost::asio::ip::host_name() ) );
    log_core->add_global_attribute( "mpi_rank",
                                    mpi_rank_attr_t( world.rank() ) );
    log_core->add_global_attribute( "timer",
                                    timer_attr_t() );

    typedef boost::mpl::vector<std::time_t, std::tm, boost::posix_time::ptime> date_time_types;

    boost::log::add_common_attributes();
    boost::log::init_log_to_console
    (
        std::cout,
        log_keywords::format =
        (
            log_fmt::stream
                << "<" << log_fmt::attr< boost::log::trivial::severity_level >("Severity") << "> "
                << "(" << log_fmt::attr<std::string>("hostname") << "-"
                    << log_fmt::attr<int>("mpi_rank") << ") "
                << log_fmt::date_time<date_time_types>("TimeStamp", log_keywords::format = "%Y.%m.%d %H:%M:%S" ) << " "
                // << log_fmt::attr< unsigned int >("LineID", log_keywords::format = "%08x")
                << log_fmt::if_(log_flt::has_attr("timer"))
                   [  log_fmt::stream << "[" << log_fmt::time_duration("timer", log_keywords::format = "%h:%M:%s" ) << "] " ]
                << log_fmt::message()
        )
    );
#endif

    LOG_INFO( "MPI environment initialized, "
              << world.size() << " process(es) in total..." );
    ChessboardBiclusteringHyperPriors  hyperpriors;
    hyperpriors.signalHyperprior.meanVarScale = 2;
    GibbsSamplerParams          gibbsParams;
    TurbineCascadeParams        cascadeParams;
    cascadeParams.turbinesCount = world.size();
    //cascadeParams.maxTemperature = 5;

    boost::scoped_ptr<CellSignalParams>         signalParams;
    boost::scoped_ptr<PrecomputedDataParams>    precomputedDataParams;
    boost::scoped_ptr<ChessboardBiclusteringPriors> priors;
    boost::scoped_ptr<PrecomputedData> precomputed;
    boost::scoped_ptr<BIMAPSampleCollectorParams>    collectorParams;
    bool params_res = false;
    bool is_collector = world.rank() == 0;
    if ( is_collector ) {
        VT_TRACER( "Read parameters" );
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
        VT_TRACER( "Read OPA Data" );
        if ( !ioParams->dataFilename.empty() ) {
            boost::filesystem::path data_file_path( ioParams->dataFilename );
            LOG_INFO( "Loading data from " << data_file_path << "..." );
            if ( !boost::filesystem::exists( data_file_path ) ) {
                THROW_RUNTIME_ERROR( "OPAData file not found: " << data_file_path );
            }
            data.reset( new OPAData( OPAData::load( data_file_path.string().c_str() ) ) );
        }
        else if ( !ioParams->proteinsFilename.empty() ) {
            data.reset( new OPAData( OPADataImportCSV( *ioParams ) ) );
        } else {
            THROW_RUNTIME_ERROR( "No input data" );
        }

        {
        LOG_INFO( "Precomputing..." );
        VT_TRACER( "Precomputing" );
        precomputed.reset( new PrecomputedData( *data, *precomputedDataParams, *signalParams ) );
        }
        }
    }
    LOG_INFO( ( is_collector ? "Broadcasting" : "Receiving" ) << " input data" );
    {
    VT_TRACER( "broadcast parameter parsing results" );
    mpi::broadcast( world, params_res, 0 );
    }
    if ( !params_res ) return ( 0 ); // --help option, no computation

    {
    VT_TRACER( "broadcast params" );
    mpi::broadcast( world, hyperpriors, 0 );
    mpi::broadcast( world, gibbsParams, 0 );
    mpi::broadcast( world, cascadeParams, 0 );
    }
    broadcast_statically_tracked( world, "data", data, 0 );
    broadcast_statically_tracked( world, "signalParams", signalParams, 0 );
    broadcast_statically_tracked( world, "priors", priors, 0 );
    broadcast_statically_tracked( world, "precomputedDataParams", precomputedDataParams, 0 );
    broadcast_statically_tracked( world, "precomputed", precomputed, 0 );
    LOG_DEBUG1_IF( !is_collector, "Received OPAData (" 
                    << data->objectsCount() << " objects, " << data->probesCount() << " probes)" );
    LOG_INFO( "Initializing sampler..." );
    BIMAPSamplerHelper         helper( *precomputed,
                                        hyperpriors, *priors, gibbsParams,
                                        world.rank() * 1981 );

    LOG_DEBUG1( "Setting initial clustering..." );
    ChessboardBiclustering iniClus;
    if ( is_collector ) iniClus = helper.trivialClustering();
    {
    VT_TRACER( "broadcast initial clustering" );
    mpi::broadcast( world, iniClus, 0 );
    }
    StdOutPTCExecutionMonitor   mon( 1 );
    if ( is_collector ) ccIndexing.reset( new ChessboardBiclusteringsIndexing() );

    res.reset( MPI_BIMAPSampler_run( helper, world, cascadeParams, iniClus, 
                                ccIndexing.get(), collectorParams.get(), &mon ) );

    LOG_INFO( "Waiting for the other process(es) to finish..." );
    world.barrier();
    LOG_INFO( "All processes finished" );
	VT_BUFFER_FLUSH();
}
    if ( res ) {
        LOG_INFO( "Process finished with " << res->stepsCount() << " samples collected" );
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
        LOG_INFO( "Process finished (not collector)" );
    }
}
