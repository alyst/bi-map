#include "cemm/bimap/BasicTypedefs.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <cemm/eesampler/ConsolePTCExecutionMonitor.h>

#include "cemm/bimap/BIMAPResultsSerialize.h"
#include "cemm/bimap/OPADataImportCSV.h"

#include "ParametersReader.h"

using namespace cemm::bimap;

int main( int argc, char* argv[] )
{
    ChessboardBiclusteringHyperPriors      hyperpriors;
    hyperpriors.signalHyperprior.meanVarScale = 2;
    GibbsSamplerParams              gibbsParams;
    TurbineCascadeParams            cascadeParams;
    CellSignalParams                signalParams;
    PrecomputedDataParams           precomputedDataParams;
    ChessboardBiclusteringPriors           priors;
    BIMAPSampleCollectorParams     collectorParams;
    BIMAPIOParams                  ioParams;

    if ( !BIMAPParamsRead( argc, argv,
                      hyperpriors, gibbsParams, cascadeParams,
                      signalParams, precomputedDataParams, priors,
                      collectorParams, ioParams ) )
    {
        // --help option was given, no computation
        return ( 0 );
    }

    OPAData data;
    if ( !ioParams.dataFilename.empty() ) {
        boost::filesystem::path data_file_path( ioParams.dataFilename );
        if ( !boost::filesystem::exists( data_file_path ) ) {
            THROW_RUNTIME_ERROR( "OPAData file not found: " << data_file_path );
        }
        LOG_INFO( "Loading data from file " << data_file_path << "..." );
        data = OPAData::load( data_file_path.string().c_str() );
    }
    else if ( !ioParams.proteinsFilename.empty() ) {
        data = OPADataImportCSV( ioParams );
    } else {
        THROW_RUNTIME_ERROR( "No input data" );
    }

    LOG_INFO( "Calculating preliminary signals" );
    PrecomputedData precomputed( data, precomputedDataParams, signalParams );

    LOG_DEBUG1( "Initializing sampler..." );
    BIMAPSamplerHelper         helper( precomputed,
                                       hyperpriors, priors, gibbsParams, 1981 );

    LOG_DEBUG1( "Setting initial clustering..." );
    ChessboardBiclustering iniClus = helper.trivialClustering();

    StdOutPTCExecutionMonitor   mon( 1 );
    ChessboardBiclusteringsIndexing    ccIndexing;

    LOG_INFO( "Running sampling..." );
    BIMAPWalk res = BIMAPSampler_run( helper, ccIndexing, gibbsParams,
                                        cascadeParams, collectorParams,
                                        iniClus, &mon );

    LOG_INFO( "Sampler finished with " << res.stepsCount() << " samples collected" );
    BOOST_ASSERT( res.check() );
    if ( !ioParams.outputFilename.empty() ) {
        LOG_INFO( "Saving walk to file " << ioParams.outputFilename << "..." );
        size_t removedSteps = res.filterSteps( ioParams.minCrossClusRefCount,
                                               ioParams.minObjectsPtnRefCount,
                                               ioParams.minProbesPtnRefCount );
        LOG_INFO_IF( removedSteps > 0, "Removed " << removedSteps << " samples of rarely encountered chessboard biclusterings" );
        ccIndexing.remove_unreferenced();
        BIMAPResultsSave( ioParams.outputFilename.c_str(),
                            data.objectsLabelMap(), data.probesLabelMap(),
                            res, NULL );
        LOG_INFO( "Saving done\n" );
    } else {
        LOG_INFO( "Walk was not saved: no walk filename specified\n" );
    }
}
