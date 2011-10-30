#include "BasicTypedefs.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "ConsolePTCExecutionMonitor.h"

#include "BIMAPResultsSerialize.h"
#include "ParametersReader.h"
#include "OPADataImportCSV.h"

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
    bool                           mapBaitsToObjects = true;

    BIMAPParamsRead( argc, argv,
                      hyperpriors, gibbsParams, cascadeParams,
                      signalParams, precomputedDataParams, priors,
                      collectorParams, ioParams, mapBaitsToObjects );

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
        boost::filesystem::path proteins_file_path( ioParams.proteinsFilename );
        if ( !boost::filesystem::exists( proteins_file_path ) ) {
            THROW_RUNTIME_ERROR( "Proteins file not found: " << proteins_file_path );
        }
        boost::filesystem::path exp_design_file_path( ioParams.expDesignFilename );
        if ( !boost::filesystem::exists( exp_design_file_path ) ) {
            THROW_RUNTIME_ERROR( "Experimental design file not found: " << exp_design_file_path );
        }
        boost::filesystem::path measurements_file_path( ioParams.measurementsFilename );
        if ( !boost::filesystem::exists( measurements_file_path ) ) {
            THROW_RUNTIME_ERROR( "Measurements file not found: " << measurements_file_path );
        }
        data = OPADataImportCSV( proteins_file_path.string().c_str(),
                                 exp_design_file_path.string().c_str(),
                                 measurements_file_path.string().c_str(),
                                 ioParams.csvColumnSeparator );
    } else {
        THROW_RUNTIME_ERROR( "No input data" );
    }

    LOG_INFO( "Calculating preliminary signals" );
    PrecomputedData precomputed( data, precomputedDataParams, signalParams );

    LOG_DEBUG1( "Initializing sampler..." );
    BIMAPSamplerHelper         helper( precomputed,
                                       hyperpriors, priors, gibbsParams, 1981 );

    LOG_DEBUG1( "Setting initial clustering..." );
    ChessboardBiclustering iniClus = helper.trivialClustering( mapBaitsToObjects );

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
