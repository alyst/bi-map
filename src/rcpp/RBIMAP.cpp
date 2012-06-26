#include <cemm/bimap/BasicTypedefs.h>

#include <cemm/debug.h>

#include <stdarg.h>
#include <cmath>

#include <Rcpp/iostream/Rostream.h>
#include <cemm/R_log.h>

#include <R_ext/Print.h>
#include <R_ext/Rdynload.h>
#include <Rcpp/S4.h>
#include <Rcpp/Vector.h>

#include <cemm/RUtils.h>
#include <cemm/containers/dynamic_bitset_foreach.h>
#include <cemm/eesampler/ConsolePTCExecutionMonitor.h>
#include <cemm/math/RngPtr.h>

#ifndef NDEBUG
//#define DYNLOAD_DEBUG
#endif

#include <cemm/bimap/OPAData.h>
#include <cemm/bimap/BIMAPSampler.h>
#include <cemm/bimap/ChessboardBiclusteringsPDFEval.h>
#include <cemm/bimap/BlocksScoring.h>
#include <cemm/bimap/BIMAPResultsSerialize.h>

#define BIMAP_PARAM_WALK_SAMPLES                    "walk.samples"
#define BIMAP_PARAM_WALK_CREATE_ROBJECT             "walk.create.RObject"
#define BIMAP_PARAM_WALK_FILE                       "walk.file"

#define BIMAP_PARAM_OBJECTS_COMPONENTS_THRESHOLD    "objects.components.threshold"
#define BIMAP_PARAM_PROBES_COMPONENTS_THRESHOLD     "probes.components.threshold"
#define BIMAP_PARAM_OBJECTS_CLOT_THRESHOLD          "objects.clot.threshold"
#define BIMAP_PARAM_STABLE_OBJECTS_CLUSTERS_THRESHOLD   "stable.objects.clusters.threshold"
#define BIMAP_PARAM_STABLE_PROBES_CLUSTERS_THRESHOLD    "stable.probes.clusters.threshold"

#define BIMAP_PARAM_EESAMPLER_LEVELS_COUNT          "eesampler.levelsCount"
#define BIMAP_PARAM_EESAMPLER_TURBINES_COUNT        "eesampler.turbinesCount"
#define BIMAP_PARAM_EESAMPLER_TEMPERATURE_MULT      "eesampler.temperatureMultiplier"
#define BIMAP_PARAM_EESAMPLER_BURNIN_ITERATIONS     "eesampler.burnInIterations"
#define BIMAP_PARAM_EESAMPLER_LADDER_ADJUST_PERIOD  "eesampler.ladderAdjustPeriod"
#define BIMAP_PARAM_EESAMPLER_ENERGY_TOLERANCE      "eesampler.energyTolerance"
#define BIMAP_PARAM_EESAMPLER_JUMP_RATE             "eesampler.jumpRate"
#define BIMAP_PARAM_EESAMPLER_GENERATE_RATE         "eesampler.generateRate"
#define BIMAP_PARAM_EESAMPLER_DETENTION_ITERATIONS  "eesampler.detentionIterations"

#define BIMAP_PARAM_INITIAL_OBJECTS_PARTITION       "ini.objects.partition"
#define BIMAP_PARAM_INITIAL_PROBES_PARTITION        "ini.probes.partition"
#define BIMAP_PARAM_INITIAL_BLOCKS                  "ini.blocks"

#define BIMAP_PARAM_SAMPLE_RATE_OBJECT_MEMBERSHIP      "sampleRate.object.membership"
#define BIMAP_PARAM_SAMPLE_RATE_OBJECTS_SPLIT_MERGE    "sampleRate.objects.splitMerge"
#define BIMAP_PARAM_SAMPLE_RATE_PROBE_MEMBERSHIP       "sampleRate.probe.membership"
#define BIMAP_PARAM_SAMPLE_RATE_PROBES_SPLIT_MERGE     "sampleRate.probes.splitMerge"
#define BIMAP_PARAM_SAMPLE_RATE_BLOCK_FLIP             "sampleRate.block.flip"
#define BIMAP_PARAM_SAMPLE_RATE_OBJECT_MULTIPLE        "sampleRate.object.multiple"
#define BIMAP_PARAM_SAMPLE_RATE_SIGNAL                 "sampleRate.signal"
#define BIMAP_PARAM_SAMPLE_PERIOD_PRIORS               "samplePeriod.priors"
#define BIMAP_PARAM_SAMPLE_PERIOD_BICLUSTERING         "samplePeriod.chessboardBiclustering"
#define BIMAP_PARAM_BLOCK_RESAMPLES                    "block.resamples"

#define BIMAP_PARAM_PRIOR_OBJECTS_CLUSTERING_CONCENTRATION      "prior.objects.clustering.concentration"
#define BIMAP_PARAM_PRIOR_OBJECTS_CLUSTERING_DISCOUNT           "prior.objects.clustering.discount"
#define BIMAP_PARAM_PRIOR_PROBES_CLUSTERING_CONCENTRATION       "prior.probes.clustering.concentration"
#define BIMAP_PARAM_PRIOR_PROBES_CLUSTERING_DISCOUNT            "prior.probes.clustering.discount"

#define BIMAP_PARAM_SIGNAL_SEQUENCE_LENGTH_FACTOR      "signal.sequence.length.factor"

#define BIMAP_PARAM_SIGNAL_SHAPE            "signal.shape"

#define BIMAP_PARAM_PRIOR_FALSE_HITS        "prior.false_hits"
#define BIMAP_PARAM_PRIOR_TRUE_MISSES       "prior.true_misses"

#define BIMAP_PARAM_PRIOR_BLOCK_ENABLED     "prior.block.enabled"

#define BIMAP_PARAM_HPRIOR_BASELINE_SIGNAL       "hyperprior.baseline"
#define BIMAP_PARAM_HPRIOR_BASELINE_SCALE        "hyperprior.baseline.scale"
#define BIMAP_PARAM_HPRIOR_SIGNAL_VAR_SHAPE      "hyperprior.signal.variance.shape"
#define BIMAP_PARAM_HPRIOR_SIGNAL_VAR_SCALE      "hyperprior.signal.variance.scale"

#define BIMAP_PARAM_LOG_PARTICLES_FILE           "log.particles.file"
#define BIMAP_PARAM_LOG_EEJUMPS_FILE             "log.eeJumps.file"

#define PRECOMP_PARAM_OBJECT_FREQ_THRESHOLD "precomp.object.freq.threshold"
#define PRECOMP_PARAM_PROBE_FREQ_THRESHOLD  "precomp.probe.freq.threshold"

#define R_CLASS_BIMAP_WALK                 "BIMAPWalk"
#define R_CLASS_BIMAP_PRIOR_PARAMS         "BIMAPPriorParams"
#define R_CLASS_BIMAP_CLUSTERING           "BIMAPClustering"
#define R_CLASS_BIMAP_CLUSTER              "BIMAPCluster"

#define R_SLOT_PROTEINS                     "proteins"
#define R_SLOT_SAMPLES                      "samples"
#define R_SLOT_CLUSTERINGS                  "clusterings"

#define R_SLOT_OBJECTS_PARTITIONS           "objects.partitions"
#define R_SLOT_OBJECTS_CLUSTERS             "objects.clusters"
#define R_SLOT_OBJECTS_CLUSTERS_INFO        "objects.clusters.info"

#define R_SLOT_PROBES_PARTITIONS            "probes.partitions"
#define R_SLOT_PROBES_CLUSTERS              "probes.clusters"
#define R_SLOT_PROBES_CLUSTERS_INFO         "probes.clusters.info"

#define R_SLOT_BLOCKS                       "blocks"
#define R_SLOT_BLOCKS_FREQUENCY             "blocks.freq"
#define R_SLOT_STABLE_BLOCKS_SCORING        "stable.blocks.scores"
#define R_SLOT_SIGNALS                      "signals"
#define R_SLOT_OBJECTS_DATA                 "objects.data"

#define R_SLOT_CLUSTERINGS_WALK             "clusterings.walk"
#define R_SLOT_PRIORS_WALK                  "priors.walk"

#define R_SLOT_OBJECTS_COMPONENTS           "objects.components"
#define R_SLOT_PROBES_COMPONENTS            "probes.components"

#define R_SLOT_OBJECTS_SUBPARTITIONS        "objects.subpartitions"
#define R_SLOT_PROBES_SUBPARTITIONS         "probes.subpartitions"

#define R_SLOT_PROB_WEIGHTS                 "prob.weights"

#define R_COLUMN_COMPONENT_INDEX            "component.index"
#define R_COLUMN_SUBPARTITION_SERIAL        "subpartition.serial"

#define R_COLUMN_CLUSTERING_SERIAL          "clustering.serial"

#define R_COLUMN_OBJECTS_PARTITION_SERIAL   "objects.partition.serial"
#define R_COLUMN_OBJECTS_CLUSTER_SERIAL     "objects.cluster.serial"
#define R_COLUMN_OBJECT_MULTIPLE            "object.multiple"
#define R_COLUMN_OBJECT                     "object"

#define R_COLUMN_PROBES_PARTITION_SERIAL    "probes.partition.serial"
#define R_COLUMN_PROBES_CLUSTER_SERIAL      "probes.cluster.serial"
#define R_COLUMN_PROBE                      "probe"

#define R_COLUMN_TOTAL_LN_PDF               "total.lnpdf"
#define R_COLUMN_OBJECTS_LN_PDF             "objects.lnpdf"
#define R_COLUMN_PROBES_LN_PDF              "probes.lnpdf"
#define R_COLUMN_CELLS_LN_PDF               "cells.lnpdf"
#define R_COLUMN_TOTAL_LN_PARTS_PDF         "total.parts.lnpdf"
#define R_COLUMN_OBJECTS_LN_PARTS_PDF       "objects.parts.lnpdf"
#define R_COLUMN_PROBES_LN_PARTS_PDF        "probes.parts.lnpdf"

#define R_COLUMN_SIGNAL                     "signal"

#define R_COLUMN_TOTAL                      "total"
#define R_COLUMN_ENABLED                    "enabled"

#define R_COLUMN_SIZE                       "size"
#define R_COLUMN_AVG_PAIRS_COOCCUR          "avg.pairs.cooccur"
#define R_COLUMN_NSTEPS                     "nsteps"
#define R_COLUMN_NSTEPS_INCLUDED            "nsteps.included"

#define R_COLUMN_BASELINE_PEAK              "baseline.peak"
#define R_COLUMN_BASELINE_SHAPE             "baseline.shape"
#define R_COLUMN_TRUE_HIT_RATE              "trueHit.rate"
#define R_COLUMN_TRUE_MISS_RATE             "trueMiss.rate"

#define R_COLUMN_STEP                       "step"
#define R_COLUMN_TIME                       "time"
#define R_COLUMN_TURBINE                    "turbine"
#define R_COLUMN_LPP                        "lpp"
#define R_COLUMN_LPP_BLOCKS                 "lpp.blocks"
#define R_COLUMN_LPP_OBJECTS_PTN            "lpp.objects.ptn"
#define R_COLUMN_LPP_PROBES_PTN             "lpp.probes.ptn"
#define R_COLUMN_LLH                        "llh"
#define R_COLUMN_LLH_QUANT                  "llh.quantitative"
#define R_COLUMN_LLH_OBJECTS_TOPO           "llh.objects.topological"
#define R_COLUMN_LLH_OBJECTS_EXP_DESIGN     "llh.objects.exp_design"
#define R_COLUMN_LLH_PROBES_TOPO            "llh.probes.topological"
#define R_COLUMN_LLH_PROBES_EXP_DESIGN      "llh.probes.exp_design"
#define R_COLUMN_BASELINE_SIGNAL            "baseline.signal"
#define R_COLUMN_BASELINE_SIGNAL_SIGMA      "baseline.signal.sigma"
#define R_SLOT_NOISE_SIGNAL                 "noise.signal"

#define R_STRINGS_AS_FACTORS                "stringsAsFactors"

#if defined(USE_BOOST_LOG)
// Called when the library is loaded and before dlopen() returns
void __attribute__ ((constructor)) RBIMAP_load()
{
    cemm::R_logging_init();

    LOG_INFO( "R BI-MAP library loaded" );
}

// Called when the library is unloaded and before dlclose() returns
void __attribute__ ((destructor)) RBIMAP_unload()
{
    cemm::R_logging_cleanup();
}
#endif

namespace cemm { namespace bimap {

/**
    Reads-in Mass-Spec data for clustering.

    @return OPAData object
 */
OPAData ReadMassSpecData(
    SEXP    proteinsDataExp,        /**< @param[in] proteins data frame: 1 - protein label, 2 - sequence length */
    SEXP    samplesDataExp,         /**< @param[in] samples data frame: 1 - sample label, 2 - bait protein label */
    SEXP    msRunsDataExp,          /**< @param[in] MS-runs data frame: 1 - MS run label, 2 - examined sample label */
    SEXP    measurementsDataExp     /**< @param[in] MS measurements data frame: 1 - MS run label, 2 - protein label, 3 - TSC */
){
    Rcpp::DataFrame proteinsFrame( proteinsDataExp );
    Rcpp::DataFrame samplesFrame( samplesDataExp );
    Rcpp::DataFrame msRunsFrame( msRunsDataExp );
    Rcpp::DataFrame measurementsFrame( measurementsDataExp );

    Rprintf( "Reading Mass-Spec data...\n" );
    // print actual data frames passed
    Rprintf_columns( proteinsFrame, "Input Proteins" );
    Rprintf_columns( samplesFrame, "Input Samples" );
    Rprintf_columns( msRunsFrame, "Input Mass-Spec runs" );
    Rprintf_columns( measurementsFrame, "Input Mass-Spec measurements" );

    OPAData    data;

    // fill-in proteins
    {
        Rprintf( "Reading proteins...\n" );
        Rcpp::StringVector proteinAc = proteinsFrame[ 0 ];
        Rcpp::IntegerVector seqLength = proteinsFrame[ 1 ];
        // protein - length
        for ( int i = 0; i < proteinAc.size(); i++ ) {
            data.addObject( (std::string)( proteinAc[ i ] ), seqLength[ i ] );
        }
    }
    // fill-in samples
    {
        Rprintf( "Reading samples...\n" );
        Rcpp::StringVector sampleCodeVec = samplesFrame[ 0 ];
        Rcpp::StringVector baitAcVec = samplesFrame[ 1 ];
        // sample-id bait-id
        for ( int i = 0; i < sampleCodeVec.size(); i++ ) {
            std::string sampleCode = (std::string)sampleCodeVec[i];
            std::string baitAc = (std::string)baitAcVec[i];
            data.addProbe( sampleCode, baitAc != "nobait" ? baitAc : std::string() );
        }
    }
    // fill-in MS-runs
    // msrun-id sample-id multiplier
    {
        Rprintf( "Reading MS runs...\n" );
        Rcpp::StringVector msrunCodeVec = msRunsFrame[ 0 ];
        Rcpp::StringVector sampleCodeVec = msRunsFrame[ 1 ];
        Rcpp::NumericVector multiplier = msRunsFrame.size() >= 3
                                       ? msRunsFrame[ 2 ] : Rcpp::NumericVector();
        for ( int i = 0; i < msrunCodeVec.size(); i++ ) {
            std::string msrunCode = (std::string)msrunCodeVec[i];
            std::string sampleCode = (std::string)sampleCodeVec[i];
            data.addAssay( msrunCode, sampleCode,
                           i <  multiplier.size() ? multiplier[i] : 1.0 );
        }
    }
    // fill-in measurements
    // protein-id msrun-id sc pc
    {
        Rprintf( "Reading measurements...\n" );
        Rcpp::StringVector msrunCode = measurementsFrame[ 0 ];
        Rcpp::StringVector proteinAc = measurementsFrame[ 1 ];
        Rcpp::IntegerVector scVec = measurementsFrame[ 2 ];
        Rcpp::IntegerVector pcVec = measurementsFrame.size() >= 4
                               ? measurementsFrame[ 3 ] : Rcpp::IntegerVector();
        for ( int i = 0; i < msrunCode.size(); i++ ) {
            int sc_value = scVec[ i ];
            int pc_value = i < pcVec.size() ? pcVec[ i ] : 0;
            data.addMeasurement( (object_label_t)proteinAc[ i ],
                                 (assay_label_t)msrunCode[ i ],
                                 OPAData::MassSpectraData( sc_value, pc_value ) );
        }
    }
    Rprintf( "MS-data read: %d proteins, %d samples, %d MS-runs\n", 
                data.objectsCount(), data.probesCount(), data.assaysCount() );

    return ( data );
}

typedef boost::unordered_map<std::string, object_clundex_t> objects_clu_code_map_t;
typedef boost::unordered_map<std::string, probe_clundex_t> probes_clu_code_map_t;

objects_clu_code_map_t ReadObjectsClusters(
    ChessboardBiclustering&    clu,
    SEXP                objsCluExp,
    const OPAData&      data
){
    Rcpp::DataFrame   cluFrame( objsCluExp );

    Rprintf( "Reading Objects clusters...\n" );
    // print actual data frames passed
    Rprintf_columns( cluFrame, "Objects clusters" );

    objects_clu_code_map_t codeMap;
    Rcpp::StringVector clusterCode = cluFrame[ 0 ];
    Rcpp::StringVector objectCode = cluFrame[ 1 ];
    Rcpp::IntegerVector multipleVec = cluFrame[ 2 ];

    for ( int i = 0; i < clusterCode.size(); i++ ) {
        std::string cluCode = (std::string)clusterCode[ i ];
        std::string objLabel = (std::string)objectCode[ i ];
        int         objMult = multipleVec[ i ];

        OPAData::const_object_ptr_t pObj = data.object( objLabel );
        if ( !pObj ) {
            throw entity_not_found( "Object", objLabel );
        }
        objects_clu_code_map_t::const_iterator  cluIt = codeMap.find( cluCode );
        if ( cluIt == codeMap.end() ) {
            object_clundex_t cluIx = clu.addObjectCluster( pObj->index() );
            codeMap.insert( cluIt, std::make_pair( cluCode, cluIx ) );
        }
        else {
            object_clundex_t cluIx = cluIt->second;
            clu.setObjectCluster( pObj->index(), cluIx );
        }
        clu.setObjectMultiple( pObj->index(), objMult );
    }
    clu.checkObjectsPartition();
    return ( codeMap );
}

probes_clu_code_map_t ReadProbesClusters(
    ChessboardBiclustering&    clu,
    SEXP                probesCluExp,
    const OPAData&      data
){
    Rcpp::DataFrame   cluFrame( probesCluExp );

    Rprintf( "Reading Probes clusters...\n" );
    // print actual data frames passed
    Rprintf_columns( cluFrame, "Probes clusters" );

    probes_clu_code_map_t codeMap;
    Rcpp::StringVector clusterCode = cluFrame[ 0 ];
    Rcpp::StringVector probeCode = cluFrame[ 1 ];

    for ( int i = 0; i < clusterCode.size(); i++ ) {
        std::string cluCode = (std::string)clusterCode[i];
        std::string probeLabel = (std::string)probeCode[i];

        OPAData::const_probe_ptr_t pProbe = data.probe( probeLabel );
        if ( !pProbe ) {
            throw entity_not_found( "Probe", probeLabel );
        }
        probes_clu_code_map_t::const_iterator  cluIt = codeMap.find( cluCode );
        if ( cluIt == codeMap.end() ) {
            probe_clundex_t cluIx = clu.addProbeCluster( pProbe->index() );
            codeMap.insert( cluIt, std::make_pair( cluCode, cluIx ) );
        }
        else {
            probe_clundex_t cluIx = cluIt->second;
            clu.setProbeCluster( pProbe->index(), cluIx );
        }
    }
    clu.checkProbesPartition();
    return ( codeMap );
}

void ReadChessboardBiclusteringBlocks(
    ChessboardBiclustering&         clu,
    const objects_clu_code_map_t&   objsCluCodeMap,
    const objects_clu_code_map_t&   probesCluCodeMap,
    SEXP                            blocksExp
){
    Rcpp::DataFrame   blocksFrame( blocksExp );

    Rprintf( "Reading chessboard biclustering blocks...\n" );
    // print actual data frames passed
    Rprintf_columns( blocksFrame, "Blocks" );
    Rcpp::StringVector objsCluCodeVec = blocksFrame[ 0 ];
    Rcpp::StringVector probesCluCodeVec = blocksFrame[ 1 ];

    for ( int i = 0; i < objsCluCodeVec.size(); i++ ) {
        std::string objsCluCode = (std::string)objsCluCodeVec[i];
        std::string probesCluCode = (std::string)probesCluCodeVec[i];

        objects_clu_code_map_t::const_iterator objsCluIt = objsCluCodeMap.find( objsCluCode );
        if ( objsCluIt == objsCluCodeMap.end() ) {
            throw entity_not_found( "Objects cluster", objsCluCode );
        }
        probes_clu_code_map_t::const_iterator probesCluIt = probesCluCodeMap.find( probesCluCode );
        if ( probesCluIt == probesCluCodeMap.end() ) {
            throw entity_not_found( "Probes cluster", probesCluCode );
        }
        clu.setBlock( objsCluIt->second, probesCluIt->second, true );
    }
}

#if 0
SEXP ConvertClusterToRObject(
    const BlockIndexed&    cluster,
    const OPAData&             data
){
    RObject rClu( R_CLASS_BIMAP_CLUSTER );

    // make proteins vector
    SEXP    proteinLabelsVec = PROTECT( Rf_allocVector( STRSXP, cluster.objectsCount() ) );
    SEXP    proteinMultiplesVec = PROTECT( Rf_allocVector( VECSXP, cluster.objectsCount() ) );
    int     proteinIx = 0;
    for ( object_set_t::const_iterator objIt = cluster.objects().begin(); objIt != cluster.objects().end(); ++objIt ){
        object_index_t  objIx = *objIt;
        SET_STRING_ELT( proteinLabelsVec, proteinIx, Rf_mkChar( data.object( objIx ).label().c_str() ) );
        SET_VECTOR_ELT( proteinMultiplesVec, proteinIx, Rf_ScalarInteger( cluster.objectMultiple( objIx ) ) );
        proteinIx++;
    }
    Rf_setAttrib( proteinMultiplesVec, R_NamesSymbol, proteinLabelsVec );

    // make samples vector
    SEXP    sampleLabelsVec = PROTECT( Rf_allocVector( STRSXP, cluster.probesCount() ) );
    SEXP    sampleSignalsVec = PROTECT( Rf_allocVector( VECSXP, cluster.probesCount() ) );
    int     sampleIx = 0;
    foreach_bit( probe_index_t, probeIx, cluster.probes() ){
        SET_STRING_ELT( sampleLabelsVec, sampleIx, Rf_mkChar( data.probe( probeIx ).label().c_str() ) );
        SET_VECTOR_ELT( sampleSignalsVec, sampleIx, Rf_ScalarReal( cluster.probeSignal( probeIx ) ) );
        sampleIx++;
    }
    Rf_setAttrib( sampleSignalsVec, R_NamesSymbol, sampleLabelsVec );

    rClu.setSlot( R_SLOT_SERIAL, Rf_ScalarInteger( cluster.serial() ) );
    rClu.setSlot( R_SLOT_PROTEINS, proteinMultiplesVec );
    rClu.setSlot( R_SLOT_SAMPLES, sampleSignalsVec );

    UNPROTECT( 4 );

    return ( rClu );
}
#endif

/**
    Converts BIMAP Priors object to R-compatible (S4) object.
 */
Rcpp::DataFrame ConvertPriorParamsWalkToRDataFrame(
    const BIMAPWalk&       walk   /** @param[in] MCMC biclustering walk to convert */
){
    Rcpp::IntegerVector priorStep( walk.priorParamsStepsCount() );
    Rcpp::NumericVector priorTime( walk.priorParamsStepsCount() );
    Rcpp::IntegerVector priorTurbine( walk.priorParamsStepsCount() );
    Rcpp::NumericVector baseline( walk.priorParamsStepsCount() );
    Rcpp::NumericVector sigma( walk.priorParamsStepsCount() );
    //Rcpp::NumericVector noise( walk.priorParamsStepsCount() );

    int elmIx = 0;

    for ( BIMAPWalk::const_priors_step_iterator stepIt = walk.priorParamsStepsBegin(); stepIt != walk.priorParamsStepsEnd(); ++stepIt ){
        const ChessboardBiclusteringDerivedPriors& priors = stepIt->priors;
        priorStep[ elmIx ] = elmIx;
        priorTime[ elmIx ] = stepIt->time;
        priorTurbine[ elmIx ] = stepIt->turbineIx;
        baseline[ elmIx ] = priors.signalPrior.mean;
        sigma[ elmIx ] = priors.signalPrior.sigma;
        //noise[ elmIx ] = priors.signalPrior.sigma;
        //rPriors.slot( R_SLOT_NOISE_SIGNAL ) = priors..posShape() ) );

        //rPriors.setSlot( R_SLOT_NOISE_SIGMA, Rf_ScalarReal( priors.noiseSignalSigma() ) );
        //rPriors.setSlot( R_SLOT_NEG_SHAPE, Rf_ScalarReal( priors.negShape() ) );
        elmIx++;
    }

    return ( Rcpp::DataFrame::create( Rcpp::Named( R_COLUMN_STEP, priorStep ),
                                      Rcpp::Named( R_COLUMN_TIME, priorTime ),
                                      Rcpp::Named( R_COLUMN_TURBINE, priorTurbine ),
                                      Rcpp::Named( R_COLUMN_BASELINE_SIGNAL, baseline ),
                                      Rcpp::Named( R_COLUMN_BASELINE_SIGNAL_SIGMA, sigma ),
                                      Rcpp::Named( R_STRINGS_AS_FACTORS, false )
 ) );
}

/**
    Converts BIMAP walk to R-compatible (S4) object.
 */
SEXP ConvertBIMAPWalkToRObject(
    const objects_label_map_type&           objects,
    const probes_label_map_type&            probes,
    const BIMAPWalk&                        walk,       /** @param[in] MCMC biclustering walk to convert */
    const ChessboardBiclusteringsPDFEval*   pdfAdjust,  /** @param[in] PDF adjustment data */
    const BlocksScoring*                    stableBlocksScoring   /** @param[in] optional scoring for on-blocks for selected objects and probes clusters */
){
    Rprintf( "Exporting BI-MAP walk...\n" );
    Rcpp::S4 rWalk( R_CLASS_BIMAP_WALK );

    {
        Rprintf( "Exporting chessboard biclusterings walk...\n" );
        Rcpp::IntegerVector stepVec( walk.stepsCount() );
        Rcpp::IntegerVector turbineVec( walk.stepsCount() );
        Rcpp::IntegerVector clusSerialVec( walk.stepsCount() );
        Rcpp::NumericVector timeVec( walk.stepsCount() );
        Rcpp::NumericVector llhVec( walk.stepsCount() );
        Rcpp::NumericVector llhQuantVec( walk.stepsCount() );
        Rcpp::NumericVector llhObjsTopoVec( walk.stepsCount() );
        Rcpp::NumericVector llhObjsExpDesignVec( walk.stepsCount() );
        Rcpp::NumericVector llhProbesTopoVec( walk.stepsCount() );
        Rcpp::NumericVector llhProbesExpDesignVec( walk.stepsCount() );
        Rcpp::NumericVector lppVec( walk.stepsCount() );
        Rcpp::NumericVector lppBlocksVec( walk.stepsCount() );
        Rcpp::NumericVector lppObjsPtnVec( walk.stepsCount() );
        Rcpp::NumericVector lppProbesPtnVec( walk.stepsCount() );
        Rcpp::NumericVector baselinePeak( walk.stepsCount() );
        Rcpp::NumericVector baselineShape( walk.stepsCount() );
        Rcpp::NumericVector trueMissRate( walk.stepsCount() );

        int     elmIx = 0;
        for ( BIMAPWalk::const_step_iterator stepIt = walk.stepsBegin(); stepIt != walk.stepsEnd(); ++stepIt ){
            const ChessboardBiclusteringIndexed& clus = stepIt->clustering;

            stepVec[ elmIx ] = elmIx;
            turbineVec[ elmIx ] = stepIt->turbineIx;

            llhVec[ elmIx ] = stepIt->metrics.llh( walk.probWeights() );
            llhQuantVec[ elmIx ] = stepIt->metrics.llhObjs.quant;
            llhObjsTopoVec[ elmIx ] = stepIt->metrics.llhObjs.topo;
            llhObjsExpDesignVec[ elmIx ] = stepIt->metrics.llhObjs.conf;
            llhProbesTopoVec[ elmIx ] = stepIt->metrics.llhProbes.topo;
            llhProbesExpDesignVec[ elmIx ] = stepIt->metrics.llhProbes.conf;

            lppVec[ elmIx ] = stepIt->metrics.lpp();
            lppBlocksVec[ elmIx ] = stepIt->metrics.lppBlocks;
            lppObjsPtnVec[ elmIx ] = stepIt->metrics.lppObjsPtn;
            lppProbesPtnVec[ elmIx ] = stepIt->metrics.lppProbesPtn;
            timeVec[ elmIx ] = stepIt->time;
            clusSerialVec[ elmIx ] = clus.serial();
            baselinePeak[ elmIx ] = clus.baselineSignalParams().lnScRate();
            baselineShape[ elmIx ] = clus.baselineSignalParams().scShape();
            trueMissRate[ elmIx ] = clus.noiseParams().successRate;
            elmIx++;
        }
        rWalk.slot( R_SLOT_CLUSTERINGS_WALK ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_STEP, stepVec ),
                Rcpp::Named( R_COLUMN_TIME, timeVec ),
                Rcpp::Named( R_COLUMN_TURBINE, turbineVec ),
                Rcpp::Named( R_COLUMN_CLUSTERING_SERIAL, clusSerialVec ),

                Rcpp::Named( R_COLUMN_LLH, llhVec ),
                Rcpp::Named( R_COLUMN_LLH_QUANT, llhQuantVec ),
                Rcpp::Named( R_COLUMN_LLH_OBJECTS_TOPO, llhObjsTopoVec ),
                Rcpp::Named( R_COLUMN_LLH_OBJECTS_EXP_DESIGN, llhObjsExpDesignVec ),
                Rcpp::Named( R_COLUMN_LLH_PROBES_TOPO, llhProbesTopoVec ),
                Rcpp::Named( R_COLUMN_LLH_PROBES_EXP_DESIGN, llhProbesExpDesignVec ),

                Rcpp::Named( R_COLUMN_LPP, lppVec ),
                Rcpp::Named( R_COLUMN_LPP_BLOCKS, lppBlocksVec ),
                Rcpp::Named( R_COLUMN_LPP_OBJECTS_PTN, lppObjsPtnVec ),
                Rcpp::Named( R_COLUMN_LPP_PROBES_PTN, lppProbesPtnVec ),

                Rcpp::Named( R_COLUMN_BASELINE_PEAK, baselinePeak ),
                Rcpp::Named( R_COLUMN_BASELINE_SHAPE, baselineShape ),
                Rcpp::Named( R_COLUMN_TRUE_MISS_RATE, trueMissRate ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
             );
    }
    {
        Rprintf( "Exporting the weights of probability components...\n" );
        Rcpp::CharacterVector weightNameVec( 4 );
        Rcpp::NumericVector   valueVec( 4 );
        weightNameVec[0] = R_COLUMN_LLH_OBJECTS_TOPO;
        valueVec[0] =  walk.probWeights().objects.topo;
        weightNameVec[1] = R_COLUMN_LLH_OBJECTS_EXP_DESIGN;
        valueVec[1] =  walk.probWeights().objects.conf;
        weightNameVec[2] = R_COLUMN_LLH_PROBES_TOPO;
        valueVec[2] =  walk.probWeights().probes.topo;
        weightNameVec[3] = R_COLUMN_LLH_PROBES_EXP_DESIGN;
        valueVec[3] =  walk.probWeights().probes.conf;
        rWalk.slot( R_SLOT_PROB_WEIGHTS ) = Rcpp::DataFrame::create(
            Rcpp::Named( "name", weightNameVec ),
            Rcpp::Named( "value", valueVec ),
            Rcpp::Named( R_STRINGS_AS_FACTORS, false )
        );
    }
    {
        Rprintf( "Exporting biclusterings structure...\n" );
        Rcpp::IntegerVector clusSerialVec( walk.indexing().size() );
        Rcpp::IntegerVector objSerialVec( walk.indexing().size() );
        Rcpp::IntegerVector probeSerialVec( walk.indexing().size() );
        Rcpp::NumericVector totalLnIcPdfVec( walk.indexing().size(), R_NaReal );
        Rcpp::NumericVector objectsLnIcPdfVec( walk.indexing().size(), R_NaReal );
        Rcpp::NumericVector probesLnIcPdfVec( walk.indexing().size(), R_NaReal );
        Rcpp::NumericVector totalLnPartsPdfVec( walk.indexing().size(), R_NaReal );
        Rcpp::NumericVector objectsLnPartsPdfVec( walk.indexing().size(), R_NaReal );
        Rcpp::NumericVector probesLnPartsPdfVec( walk.indexing().size(), R_NaReal );
        Rcpp::NumericVector cellsLnPdfVec( walk.indexing().size(), R_NaReal );

        int     elmIx = 0;
        for ( ChessboardBiclusteringsIndexing::const_value_iterator it = walk.indexing().valueMap().begin(); it != walk.indexing().valueMap().end(); ++it ) {
            const ChessboardBiclusteringScaffold& clus = (*it)->value();
            clusSerialVec[ elmIx ] = (*it)->serial();
            objSerialVec[ elmIx ] = clus.pObjectsPartition->serial();
            probeSerialVec[ elmIx ] = clus.pProbesPartition->serial();
            if ( pdfAdjust ) {
                log_prob_t lnPdf = pdfAdjust->lnIcPdf( clus );
                if ( !std::isnan( lnPdf ) ) totalLnIcPdfVec[ elmIx ] =  lnPdf;
                lnPdf = pdfAdjust->lnCellsMaskPdf( clus );
                if ( !std::isnan( lnPdf ) ) cellsLnPdfVec[ elmIx ] =  lnPdf;
                lnPdf = pdfAdjust->lnObjectsPartitionIcPdf( clus );
                if ( !std::isnan( lnPdf ) ) objectsLnIcPdfVec[ elmIx ] =  lnPdf;
                lnPdf = pdfAdjust->lnProbesPartitionIcPdf( clus );
                if ( !std::isnan( lnPdf ) ) probesLnIcPdfVec[ elmIx ] =  lnPdf;

                lnPdf = pdfAdjust->lnPartsPdf( clus );
                if ( !std::isnan( lnPdf ) ) totalLnPartsPdfVec[ elmIx ] =  lnPdf;
                lnPdf = pdfAdjust->lnObjectsPartitionPartsPdf( clus );
                if ( !std::isnan( lnPdf ) ) objectsLnPartsPdfVec[ elmIx ] =  lnPdf;
                lnPdf = pdfAdjust->lnProbesPartitionPartsPdf( clus );
                if ( !std::isnan( lnPdf ) ) probesLnPartsPdfVec[ elmIx ] =  lnPdf;
            }
            elmIx++;
        }
        Rprintf( "Creating biclusterings dataframe, nrow=%d\n", clusSerialVec.size() );
        rWalk.slot( R_SLOT_CLUSTERINGS ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_CLUSTERING_SERIAL, clusSerialVec ),
                Rcpp::Named( R_COLUMN_OBJECTS_PARTITION_SERIAL, objSerialVec ),
                Rcpp::Named( R_COLUMN_PROBES_PARTITION_SERIAL, probeSerialVec ),
                Rcpp::Named( R_COLUMN_TOTAL_LN_PDF, totalLnIcPdfVec ),
                Rcpp::Named( R_COLUMN_OBJECTS_LN_PDF, objectsLnIcPdfVec ),
                Rcpp::Named( R_COLUMN_PROBES_LN_PDF, probesLnIcPdfVec ),
                Rcpp::Named( R_COLUMN_CELLS_LN_PDF, cellsLnPdfVec ),
                Rcpp::Named( R_COLUMN_TOTAL_LN_PARTS_PDF, totalLnPartsPdfVec ),
                Rcpp::Named( R_COLUMN_OBJECTS_LN_PARTS_PDF, objectsLnPartsPdfVec ),
                Rcpp::Named( R_COLUMN_PROBES_LN_PARTS_PDF, probesLnPartsPdfVec ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
            );
    }
    {
        Rprintf( "Exporting objects partitions...\n" );
        std::vector<int> ptnSerial;
        std::vector<int> cluSerial;
        // iterate partitions
        for ( ChessboardBiclusteringsIndexing::object_partition_indexing::const_value_iterator it = walk.indexing().objectPartitionIndexing().valueMap().begin(); 
              it != walk.indexing().objectPartitionIndexing().valueMap().end(); ++it ) {
            typedef ChessboardBiclusteringsIndexing::object_partition_indexing::collection_type ptn_type;
            const ptn_type& ptn = (*it)->value();
            // iterate clusters of partitions
            for ( ptn_type::const_iterator elit = ptn.begin(); elit != ptn.end(); ++elit ) {
                ptnSerial.push_back( (*it)->serial() );
                cluSerial.push_back( (*elit)->serial() );
            }
        }
        Rprintf( "Creating objects.partitions dataframe, nrow=%d\n", cluSerial.size() );
        rWalk.slot( R_SLOT_OBJECTS_PARTITIONS ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_OBJECTS_PARTITION_SERIAL, ptnSerial ),
                Rcpp::Named( R_COLUMN_OBJECTS_CLUSTER_SERIAL, cluSerial ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
            );
    }
    if ( pdfAdjust ) {
        Rprintf( "Exporting objects independent components...\n" );
        std::vector<std::string> objectIds;
        std::vector<size_t> compIxs;

        for ( size_t i = 0; i < pdfAdjust->objectComponents().size(); i++ ) {
            const object_set_t& component = pdfAdjust->objectComponents()[ i ];
            for( object_set_t::const_iterator oit = component.begin(); oit != component.end(); ++oit ) {
                objectIds.push_back( objects.find( *oit )->second );
                compIxs.push_back( i );
            }
        }
        Rprintf( "Creating objects.components dataframe, nrow=%d\n", objectIds.size() );
        rWalk.slot( R_SLOT_OBJECTS_COMPONENTS ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_OBJECT, objectIds ),
                Rcpp::Named( R_COLUMN_COMPONENT_INDEX, compIxs ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
            );

        Rprintf( "Exporting objects subpartitions...\n" );
        Rcpp::IntegerVector ptnSerialVec( walk.indexing().objectPartitionIndexing().valueMap().size() * pdfAdjust->objectComponents().size() );
        Rcpp::IntegerVector compIxSerialVec( ptnSerialVec.size() );
        Rcpp::IntegerVector subptnSerialVec( ptnSerialVec.size() );

        size_t ix = 0;
        for ( ChessboardBiclusteringsIndexing::object_partition_indexing::const_value_iterator it = walk.indexing().objectPartitionIndexing().valueMap().begin(); 
            it != walk.indexing().objectPartitionIndexing().valueMap().end(); ++it
        ){
            size_t ptnSerial = (*it)->serial();
            const PartitionIndependentComponentsPDF::serial_vector_type& serials = pdfAdjust->objectsIcPdf().subpartitions( ptnSerial );
            if ( !serials.empty() && serials.size() != pdfAdjust->objectComponents().size() ) {
                THROW_RUNTIME_ERROR( "Incorrect size of components composition for objects partition #" << ptnSerial );
            }
            for ( size_t compIx = 0; compIx < pdfAdjust->objectComponents().size(); compIx++ ) {
                size_t longIx = ix * pdfAdjust->objectComponents().size() + compIx;
                ptnSerialVec[ longIx ] = ptnSerial;
                compIxSerialVec[ longIx ] = compIx;
                if ( !serials.empty() ) subptnSerialVec[ longIx ] = serials[ compIx ];
                else subptnSerialVec[ longIx ] = R_NaInt;
            }
            ix++;
        }
        Rprintf( "Creating objects.subpartitions dataframe, nrow=%d\n", ptnSerialVec.size() );
        rWalk.slot( R_SLOT_OBJECTS_SUBPARTITIONS ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_OBJECTS_PARTITION_SERIAL, ptnSerialVec ),
                Rcpp::Named( R_COLUMN_COMPONENT_INDEX, compIxSerialVec ),
                Rcpp::Named( R_COLUMN_SUBPARTITION_SERIAL, subptnSerialVec ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
            );
    }
    {
        Rprintf( "Exporting objects clusters...\n" );
        std::vector<int> cluSerial;
        std::vector<std::string> objectId;
        for ( ChessboardBiclusteringsIndexing::object_partition_indexing::element_indexing_type::const_value_iterator it = walk.indexing().objectPartitionIndexing().elementIndexing().valueMap().begin(); 
              it != walk.indexing().objectPartitionIndexing().elementIndexing().valueMap().end(); ++it ) {
            const object_set_t& clu = (*it)->value();
            for ( object_set_t::const_iterator elit = clu.begin(); elit != clu.end(); ++elit ) {
                cluSerial.push_back( (*it)->serial() );
                objectId.push_back( objects.find( *elit )->second );
            }
        }
        Rprintf( "Creating objects.clusters dataframe, nrow=%d\n", cluSerial.size() );
        Rcpp::IntegerVector cluSerialVec( cluSerial.size() );
        cluSerialVec = cluSerial;
        Rcpp::StringVector objectIdVec( objectId.size() );
        objectIdVec = objectId;
        rWalk.slot( R_SLOT_OBJECTS_CLUSTERS ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_OBJECTS_CLUSTER_SERIAL, cluSerialVec ),
                Rcpp::Named( R_COLUMN_OBJECT, objectIdVec ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
            );
    }
    {
        Rprintf( "Exporting probes partitions...\n" );
        std::vector<int> ptnSerial;
        std::vector<int> cluSerial;
        // iterate partitions
        for ( ChessboardBiclusteringsIndexing::probe_partition_indexing::const_value_iterator it = walk.indexing().probePartitionIndexing().valueMap().begin(); 
              it != walk.indexing().probePartitionIndexing().valueMap().end(); ++it ) {
            typedef ChessboardBiclusteringsIndexing::probe_partition_indexing::collection_type ptn_type;
            const ptn_type& ptn = (*it)->value();
            // iterate clusters of partition
            for ( ptn_type::const_iterator elit = ptn.begin(); elit != ptn.end(); ++elit ) {
                ptnSerial.push_back( (*it)->serial() );
                cluSerial.push_back( (*elit)->serial() );
            }
        }
        Rprintf( "Creating probes.partitions dataframe, nrow=%d\n", cluSerial.size() );
        rWalk.slot( R_SLOT_PROBES_PARTITIONS ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_PROBES_PARTITION_SERIAL, ptnSerial ),
                Rcpp::Named( R_COLUMN_PROBES_CLUSTER_SERIAL, cluSerial ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
            );
    }
    if ( pdfAdjust ) {
        Rprintf( "Exporting probes independent components...\n" );
        std::vector<std::string> probeIds;
        std::vector<size_t> compIxs;

        for ( size_t i = 0; i < pdfAdjust->probeComponents().size(); i++ ) {
            const probe_bitset_t& component = pdfAdjust->probeComponents()[ i ];
            foreach_bit( probe_index_t, probeIx, component ) {
                probeIds.push_back( probes.find( probeIx )->second );
                compIxs.push_back( i );
            }
        }
        Rprintf( "Creating probes.components dataframe, nrow=%d\n", probeIds.size() );
        rWalk.slot( R_SLOT_PROBES_COMPONENTS ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_PROBE, probeIds ),
                Rcpp::Named( R_COLUMN_COMPONENT_INDEX, compIxs ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
            );

        Rprintf( "Exporting probes subpartitions and adjusted PDF...\n" );
        Rcpp::IntegerVector ptnSerialVec( walk.indexing().probePartitionIndexing().valueMap().size() * pdfAdjust->probeComponents().size() );
        Rcpp::IntegerVector compIxSerialVec( ptnSerialVec.size() );
        Rcpp::IntegerVector subptnSerialVec( ptnSerialVec.size() );

        size_t ix = 0;
        for ( ChessboardBiclusteringsIndexing::probe_partition_indexing::const_value_iterator it = walk.indexing().probePartitionIndexing().valueMap().begin(); 
            it != walk.indexing().probePartitionIndexing().valueMap().end(); ++it
        ){
            size_t ptnSerial = (*it)->serial();
            const PartitionIndependentComponentsPDF::serial_vector_type& serials = pdfAdjust->probesIcPdf().subpartitions( ptnSerial );
            if ( !serials.empty() && serials.size() != pdfAdjust->probeComponents().size() ) {
                THROW_RUNTIME_ERROR( "Incorrect size of components composition for probes partition #" << ptnSerial );
            }
            for ( size_t compIx = 0; compIx < pdfAdjust->probeComponents().size(); compIx++ ) {
                size_t longIx = ix * pdfAdjust->probeComponents().size() + compIx;
                ptnSerialVec[ longIx ] = ptnSerial;
                compIxSerialVec[ longIx ] = compIx;
                if ( !serials.empty() ) subptnSerialVec[ longIx ] = serials[ compIx ];
                else subptnSerialVec[ longIx ] = R_NaInt;
            }
            ix++;
        }
        Rprintf( "Creating probes.subpartitions dataframe, nrow=%d\n", ptnSerialVec.size() );
        rWalk.slot( R_SLOT_PROBES_SUBPARTITIONS ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_PROBES_PARTITION_SERIAL, ptnSerialVec ),
                Rcpp::Named( R_COLUMN_COMPONENT_INDEX, compIxSerialVec ),
                Rcpp::Named( R_COLUMN_SUBPARTITION_SERIAL, subptnSerialVec ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
            );
    }
    {
        Rprintf( "Exporting probes clusters...\n" );
        std::vector<int> cluSerial;
        std::vector<std::string> probeId;
        for ( ChessboardBiclusteringsIndexing::probe_partition_indexing::element_indexing_type::const_value_iterator it = walk.indexing().probePartitionIndexing().elementIndexing().valueMap().begin(); 
              it != walk.indexing().probePartitionIndexing().elementIndexing().valueMap().end(); ++it ) {
            const probe_bitset_t& clu = (*it)->value();
            foreach_bit( probe_index_t, probeIx, clu ) {
                cluSerial.push_back( (*it)->serial() );
                probeId.push_back( probes.find( probeIx )->second );
            }
        }
        Rprintf( "Creating probe clusters dataframe, nrow=%d\n", cluSerial.size() );
        Rcpp::IntegerVector cluSerialVec( cluSerial.size() );
        cluSerialVec = cluSerial;
        Rcpp::StringVector probeIdVec( probeId.size() );
        probeIdVec = probeId;
        rWalk.slot( R_SLOT_PROBES_CLUSTERS ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_PROBES_CLUSTER_SERIAL, cluSerialVec ),
                Rcpp::Named( R_COLUMN_PROBE, probeIdVec ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
            );
    }
    {
        Rprintf( "Exporting chessboard biclustering blocks...\n" );
        std::vector<int> clusSerial;
        std::vector<int> objectSetSerial;
        std::vector<int> probeSetSerial;

        for ( ChessboardBiclusteringsIndexing::const_value_iterator it = walk.indexing().valueMap().begin(); it != walk.indexing().valueMap().end(); ++it ) {
            //Rprintf( "Getting scaffold %d...\n", (*it)->serial() );
            const ChessboardBiclusteringScaffold& clus = (*it)->value();
            for ( ChessboardBiclusteringIndexed::objects_cluster_collection_type::const_iterator objCluIt = clus.pObjectsPartition->value().begin(); 
                objCluIt != clus.pObjectsPartition->value().end(); ++objCluIt ) {
                //Rprintf( "Getting object's cluster serial...\n" );
                const ChessboardBiclusteringIndexed::objects_cluster_serial_type objCluSerial = (*objCluIt)->serial();
                //Rprintf( "Object's cluster serial is %d\n", objCluSerial );
                for ( ChessboardBiclusteringIndexed::probes_cluster_collection_type::const_iterator probeCluIt = clus.pProbesPartition->value().begin(); 
                probeCluIt != clus.pProbesPartition->value().end(); ++probeCluIt ) {
                    //Rprintf( "Getting probe's cluster serial...\n" );
                    const ChessboardBiclusteringIndexed::probes_cluster_serial_type probeCluSerial = (*probeCluIt)->serial();
                    //Rprintf( "Probes's cluster serial is %d\n", probeCluSerial );
                    if ( clus.isBlockEnabled( objCluSerial, probeCluSerial ) ) {
                        //Rprintf( "Block is enabled\n" );
                        clusSerial.push_back( (*it)->serial() );
                        objectSetSerial.push_back( objCluSerial );
                        probeSetSerial.push_back( probeCluSerial );
                    }
                }
            }
        }
        Rprintf( "Creating blocks dataframe, nrow=%d\n", clusSerial.size() );
        Rcpp::IntegerVector clusSerialVec( clusSerial.size() );
        clusSerialVec = clusSerial;
        Rcpp::IntegerVector objectSetVec( objectSetSerial.size() );
        objectSetVec = objectSetSerial;
        Rcpp::IntegerVector probeSetVec( probeSetSerial.size() );
        probeSetVec = probeSetSerial;
        rWalk.slot( R_SLOT_BLOCKS ) = Rcpp::DataFrame::create(
            Rcpp::Named( R_COLUMN_CLUSTERING_SERIAL, clusSerialVec ),
            Rcpp::Named( R_COLUMN_OBJECTS_CLUSTER_SERIAL, objectSetVec ),
            Rcpp::Named( R_COLUMN_PROBES_CLUSTER_SERIAL, probeSetVec ),
            Rcpp::Named( R_STRINGS_AS_FACTORS, false )
        );
    }
    {
        Rprintf( "Exporting objects data...\n" );
        std::vector<int> step;
        std::vector<std::string> object;
        std::vector<int> multiple;
        step.reserve( walk.stepsCount() * objects.size() );
        object.reserve( walk.stepsCount() * objects.size() );
        size_t stepIx = 0;
        for ( BIMAPWalk::const_step_iterator stepIt = walk.stepsBegin(); stepIt != walk.stepsEnd(); ++stepIt ){
            const ChessboardBiclusteringIndexed& clus = stepIt->clustering;
            for ( object_index_t objIx = 0; objIx < objects.size(); ++objIx ){
                step.push_back( stepIx );
                object.push_back( objects.find( objIx )->second );
                multiple.push_back( clus.objectsData()[ objIx ] );
            }
            stepIx++;
        }
        Rprintf( "Creating object data dataframe, nrow=%d\n", step.size() );
        Rcpp::IntegerVector stepVec( step.size() );
        stepVec = step;
        Rcpp::StringVector objectVec( object.size() );
        objectVec = object;
        Rcpp::IntegerVector multipleVec( multiple.size() );
        multipleVec = multiple;
        rWalk.slot( R_SLOT_OBJECTS_DATA ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_STEP, stepVec ),
                Rcpp::Named( R_COLUMN_OBJECT, objectVec ),
                Rcpp::Named( R_COLUMN_OBJECT_MULTIPLE, multipleVec ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
        );
    }
    if ( pdfAdjust ) {
        {
        Rprintf( "Exporting object clusters info...\n" );
        const ChessboardBiclusteringsPDFEval::obj_clu_stats_map& ocStatsMap = pdfAdjust->objectsClustersStatsMap();
        Rcpp::IntegerVector objectsSetSerialVec( ocStatsMap.size() );
        Rcpp::IntegerVector sizeVec( ocStatsMap.size() );
        Rcpp::IntegerVector nStepsVec( ocStatsMap.size() );
        Rcpp::IntegerVector nStepsIncludedVec( ocStatsMap.size() );
        Rcpp::NumericVector avgPairsCooccurVec( ocStatsMap.size() );

        size_t ix = 0;
        for ( ChessboardBiclusteringsPDFEval::obj_clu_stats_map::const_iterator
              it = ocStatsMap.begin(); it != ocStatsMap.end(); ++it
        ){
            objectsSetSerialVec[ ix ] = it->first;
            sizeVec[ ix ] = it->second.size;
            nStepsVec[ ix ] = it->second.nsteps;
            nStepsIncludedVec[ ix ] = it->second.nstepsIncluded;
            avgPairsCooccurVec[ ix ] = it->second.avgPairCoOccurrence;
            ix++;
        }
        Rprintf( "Creating objects clusters info dataframe, nrow=%d\n", objectsSetSerialVec.size() );
        Rcpp::DataFrame ocInfo = Rcpp::DataFrame::create(
            Rcpp::Named( R_COLUMN_OBJECTS_CLUSTER_SERIAL, objectsSetSerialVec ),
            Rcpp::Named( R_COLUMN_SIZE, sizeVec ),
            Rcpp::Named( R_COLUMN_NSTEPS, nStepsVec ),
            Rcpp::Named( R_COLUMN_NSTEPS_INCLUDED, nStepsIncludedVec ),
            Rcpp::Named( R_COLUMN_AVG_PAIRS_COOCCUR, avgPairsCooccurVec )
        );
        ocInfo.attr( "row.names" ) = objectsSetSerialVec;
        rWalk.slot( R_SLOT_OBJECTS_CLUSTERS_INFO ) = ocInfo;
        }

        {
        Rprintf( "Exporting probes clusters info...\n" );
        const ChessboardBiclusteringsPDFEval::probe_clu_stats_map& scStatsMap = pdfAdjust->probesClustersStatsMap();
        Rcpp::IntegerVector probesSetSerialVec( scStatsMap.size() );
        Rcpp::IntegerVector sizeVec( scStatsMap.size() );
        Rcpp::IntegerVector nStepsVec( scStatsMap.size() );
        Rcpp::IntegerVector nStepsIncludedVec( scStatsMap.size() );
        Rcpp::NumericVector avgPairsCooccurVec( scStatsMap.size() );

        size_t ix = 0;
        for ( ChessboardBiclusteringsPDFEval::probe_clu_stats_map::const_iterator
              it = scStatsMap.begin(); it != scStatsMap.end(); ++it
        ){
            probesSetSerialVec[ ix ] = it->first;
            sizeVec[ ix ] = it->second.size;
            nStepsVec[ ix ] = it->second.nsteps;
            nStepsIncludedVec[ ix ] = it->second.nstepsIncluded;
            avgPairsCooccurVec[ ix ] = it->second.avgPairCoOccurrence;
            ix++;
        }
        Rprintf( "Creating probes clusters info dataframe, nrow=%d\n", probesSetSerialVec.size() );
        Rcpp::DataFrame scInfo = Rcpp::DataFrame::create(
            Rcpp::Named( R_COLUMN_PROBES_CLUSTER_SERIAL, probesSetSerialVec ),
            Rcpp::Named( R_COLUMN_SIZE, sizeVec ),
            Rcpp::Named( R_COLUMN_NSTEPS, nStepsVec ),
            Rcpp::Named( R_COLUMN_NSTEPS_INCLUDED, nStepsIncludedVec ),
            Rcpp::Named( R_COLUMN_AVG_PAIRS_COOCCUR, avgPairsCooccurVec )
        );
        scInfo.attr( "row.names" ) = probesSetSerialVec;
        rWalk.slot( R_SLOT_PROBES_CLUSTERS_INFO ) = scInfo;
        }

        {
        Rprintf( "Exporting chessboard biclustering blocks \"on\" state frequency...\n" );
        const ChessboardBiclusteringsPDFEval::block_stats_map& blockStatsMap = pdfAdjust->blocksStatsMap();
        std::size_t blockStatsMapSize = 0;
        for ( ChessboardBiclusteringsPDFEval::block_stats_map::const_iterator mapit = blockStatsMap.begin(); mapit != blockStatsMap.end(); ++mapit ) {
            blockStatsMapSize += mapit->second->size();
        }
        Rcpp::IntegerVector objectsSetSerialVec( blockStatsMapSize );
        Rcpp::IntegerVector probesSetSerialVec( blockStatsMapSize );
        Rcpp::IntegerVector totalVec( blockStatsMapSize );
        Rcpp::IntegerVector enabledVec( blockStatsMapSize );

        size_t ix = 0;
        for ( ChessboardBiclusteringsPDFEval::block_stats_map::const_iterator mapit = blockStatsMap.begin(); mapit != blockStatsMap.end(); ++mapit ) {
            const ChessboardBiclusteringsPDFEval::object_block_stats_map& objMap = *mapit->second;
            for ( ChessboardBiclusteringsPDFEval::object_block_stats_map::const_iterator it = objMap.begin(); it != objMap.end(); ++it ) {
                objectsSetSerialVec[ ix ] = it->first;
                probesSetSerialVec[ ix ] = mapit->first;
                totalVec[ ix ] = it->second.count_total;
                enabledVec[ ix ] = it->second.count_on;
                ix++;
            }
        }
        Rprintf( "Creating blocks frequency dataframe, nrow=%d\n", objectsSetSerialVec.size() );
        rWalk.slot( R_SLOT_BLOCKS_FREQUENCY ) = Rcpp::DataFrame::create(
            Rcpp::Named( R_COLUMN_OBJECTS_CLUSTER_SERIAL, objectsSetSerialVec ),
            Rcpp::Named( R_COLUMN_PROBES_CLUSTER_SERIAL, probesSetSerialVec ),
            Rcpp::Named( R_COLUMN_TOTAL, totalVec ),
            Rcpp::Named( R_COLUMN_ENABLED, enabledVec )
        );
        }
    }
    {
        Rprintf( "Exporting signals...\n" );
        std::vector<int> step;
        std::vector<int> objectsSetSerial;
        std::vector<int> probesSetSerial;
        std::vector<signal_t> signal;

        const size_t estSize = walk.stepsCount() * objects.size() * probes.size() / 10;
        step.reserve( estSize );
        objectsSetSerial.reserve( estSize );
        probesSetSerial.reserve( estSize );
        signal.reserve( estSize );

        size_t stepIx = 0;
        for ( BIMAPWalk::const_step_iterator stepIt = walk.stepsBegin(); stepIt != walk.stepsEnd(); ++stepIt ){
            const ChessboardBiclusteringIndexed& clus = stepIt->clustering;
            clus.check();
            for ( ChessboardBiclusteringIndexed::block_data_map_type::const_iterator ccDataIt = clus.blocksData().begin(); 
                ccDataIt != clus.blocksData().end(); ++ccDataIt
            ){
                step.push_back( stepIx );
                objectsSetSerial.push_back( ccDataIt->first.first );
                probesSetSerial.push_back( ccDataIt->first.second );
                signal.push_back( ccDataIt->second );
            }
            stepIx++;
        }

        Rprintf( "Creating signals dataframe, nrow=%d\n", step.size() );
        Rcpp::IntegerVector stepVec( step.size() );
        stepVec = step;
        Rcpp::IntegerVector objectsSetSerialVec( objectsSetSerial.size() );
        objectsSetSerialVec = objectsSetSerial;
        Rcpp::IntegerVector probesSetSerialVec( probesSetSerial.size() );
        probesSetSerialVec = probesSetSerial;
        Rcpp::NumericVector signalVec( signal.size() );
        signalVec = signal;
        rWalk.slot( R_SLOT_SIGNALS ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_STEP, stepVec ),
                Rcpp::Named( R_COLUMN_OBJECTS_CLUSTER_SERIAL, objectsSetSerialVec ),
                Rcpp::Named( R_COLUMN_PROBES_CLUSTER_SERIAL, probesSetSerialVec ),
                Rcpp::Named( R_COLUMN_SIGNAL, signalVec ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
        );
    }

    Rprintf( "%d clustering steps exported\n", walk.stepsCount() );

    if ( stableBlocksScoring ) {
        LOG_INFO( "Exporting stable blocks scoring..." );
        Rcpp::IntegerVector objCluVec( stableBlocksScoring->blockStats().size() );
        Rcpp::IntegerVector probeCluVec( objCluVec.size() );
        Rcpp::NumericVector totalVec( objCluVec.size() );
        Rcpp::NumericVector enabledVec( objCluVec.size() );

        size_t i = 0;
        for ( BlocksScoring::block_stats_map::const_iterator it = stableBlocksScoring->blockStats().begin();
              it != stableBlocksScoring->blockStats().end(); ++it
        ){
            objCluVec[i] = it->first.first;
            probeCluVec[i] = it->first.second;
            totalVec[i] = it->second.count_total;
            enabledVec[i] = it->second.count_on;
        }
        rWalk.slot( R_SLOT_STABLE_BLOCKS_SCORING ) = Rcpp::DataFrame::create(
                Rcpp::Named( R_COLUMN_OBJECTS_CLUSTER_SERIAL, objCluVec ),
                Rcpp::Named( R_COLUMN_PROBES_CLUSTER_SERIAL, probeCluVec ),
                Rcpp::Named( R_COLUMN_TOTAL, totalVec ),
                Rcpp::Named( R_COLUMN_ENABLED, enabledVec ),
                Rcpp::Named( R_STRINGS_AS_FACTORS, false )
        );
    }

    // convert priors
    Rprintf( "Exporting priors...\n" );
    rWalk.slot( R_SLOT_PRIORS_WALK ) = ConvertPriorParamsWalkToRDataFrame( walk );
    Rprintf( "%d priors exported\n", walk.priorParamsStepsCount() );

    Rprintf( "BI-MAP result exporting done\n" );

    return ( rWalk );
}

class REESamplerExecutionMonitor: public ConsolePTCExecutionMonitor
{
protected:
    virtual void log_printf( const char* format, ... )
    {
        va_list ap;
        va_start( ap, format );
        Rvprintf( format, ap );
        va_end( ap );
    }
public:
    REESamplerExecutionMonitor(
        log_level_t verbosity = 0,
        const char* particleLog = NULL,
        const char* eeJumpLog = NULL
    ) : ConsolePTCExecutionMonitor( verbosity, particleLog, eeJumpLog )
    {}
    virtual ~REESamplerExecutionMonitor()
    {}
};

BlocksScoring* BIMAPBlocksScoring(
    const ChessboardBiclusteringsPDFEval& pdfAdjust,
    const IndexedPartitionsCollection<ObjectsCluster>& objPtnColl,
    const IndexedPartitionsCollection<ProbesCluster>& probePtnColl,
    prob_t stableObjectsClustersThreshold,
    prob_t stableProbesClustersThreshold
){
        LOG_INFO( "Finding stable clusters..." );
        std::vector<object_clundex_t> stableObjectsClusters;

        std::size_t minObjIncludedSteps = stableObjectsClustersThreshold * objPtnColl.walk().stepsCount();
        const part_stats_map_type& objectsClustersStats = pdfAdjust.objectsClustersStatsMap();
        for ( part_stats_map_type::const_iterator ocit = objectsClustersStats.begin();
            ocit != objectsClustersStats.end(); ++ocit
        ){
            if ( ocit->second.nstepsIncluded > minObjIncludedSteps ) {
                stableObjectsClusters.push_back( ocit->first );
            }
        }
        LOG_INFO( stableObjectsClusters << " stable object clusters found" );

        std::vector<probe_clundex_t> stableProbesClusters;
        std::size_t minProbeIncludedSteps = stableProbesClustersThreshold * probePtnColl.walk().stepsCount();
        const part_stats_map_type& probesClustersStats = pdfAdjust.probesClustersStatsMap();
        for ( part_stats_map_type::const_iterator pcit = probesClustersStats.begin();
            pcit != probesClustersStats.end(); ++pcit
        ){
            if ( pcit->second.nstepsIncluded > minProbeIncludedSteps ) {
                stableProbesClusters.push_back( pcit->first );
            }
        }
        LOG_INFO( stableObjectsClusters << " stable probes clusters found" );

        LOG_INFO( "Calculating stable blocks scores..." );
        return ( new BlocksScoring( pdfAdjust, objPtnColl, probePtnColl,
                                    stableObjectsClusters, stableProbesClusters ) );
}

/**
    Reads AP-MS data and runs BI-MAP.
 */
RcppExport SEXP BIMAPWalkEval(
    SEXP    proteinsDataExp,        /**< @param[in] proteins data frame: 1 - protein label, 2 - sequence length */
    SEXP    samplesDataExp,         /**< @param[in] samples data frame: 1 - sample label, 2 - bait protein label */
    SEXP    msRunsDataExp,          /**< @param[in] MS-runs data frame: 1 - MS run label, 2 - examined sample label */
    SEXP    measurementsDataExp,    /**< @param[in] MS measurements data frame: 1 - MS run label, 2 - protein label, 3 - TSC */
    SEXP    paramsExp               /**< @param[in] method params */
){
    BEGIN_RCPP

    OPAData        data = ReadMassSpecData( proteinsDataExp, samplesDataExp, msRunsDataExp, measurementsDataExp );

    PrecomputedDataParams       precomputedDataParams;
    GibbsSamplerParams          samplerParams;
    TurbineCascadeParams        cascadeParams;
    BIMAPSampleCollectorParams collectorParams;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    ChessboardBiclusteringPriors       priors;
    CellSignalParams            signalParams;
    ChessboardBiclustering iniClus( data.objectsCount(), data.probesCount() );
    std::string                 walkFilename;
    std::string                 particleLogFilename;
    std::string                 eeJumpLogFilename;
    bool                        createWalkRObject = true;
    double objectsComponentsThreshold = R_NaReal;
    double probesComponentsThreshold = R_NaReal;
    double objectsClotThreshold = R_NaReal;
    double stableObjectsClustersThreshold = R_NaReal;
    double stableProbesClustersThreshold = R_NaReal;

    {
        Rprintf( "Reading sampler and model parameters...\n" );
        RcppParamVector      params( paramsExp );

        MaybeReadRParam( walkFilename, params, BIMAP_PARAM_WALK_FILE );
        MaybeReadRParam( createWalkRObject, params, BIMAP_PARAM_WALK_CREATE_ROBJECT );
        MaybeReadRParam( collectorParams.walkSamples, params, BIMAP_PARAM_WALK_SAMPLES );

        MaybeReadRParam( cascadeParams.turbineParams.eeJumpEnergyTolerance, params, BIMAP_PARAM_EESAMPLER_ENERGY_TOLERANCE );
        MaybeReadRParam( cascadeParams.turbineParams.eeJumpRate, params, BIMAP_PARAM_EESAMPLER_JUMP_RATE );
        MaybeReadRParam( cascadeParams.turbineParams.generateRate, params, BIMAP_PARAM_EESAMPLER_GENERATE_RATE );
        MaybeReadRParam( cascadeParams.turbineParams.detentionIterations, params, BIMAP_PARAM_EESAMPLER_DETENTION_ITERATIONS );

        MaybeReadRParam( cascadeParams.levelsCount, params, BIMAP_PARAM_EESAMPLER_LEVELS_COUNT );
        MaybeReadRParam( cascadeParams.turbinesCount, params, BIMAP_PARAM_EESAMPLER_TURBINES_COUNT );
        MaybeReadRParam( cascadeParams.burnInIterations, params, BIMAP_PARAM_EESAMPLER_BURNIN_ITERATIONS );
        MaybeReadRParam( cascadeParams.ladderAdjustPeriod, params, BIMAP_PARAM_EESAMPLER_LADDER_ADJUST_PERIOD );
        MaybeReadRParam( cascadeParams.temperatureMultiplier, params, BIMAP_PARAM_EESAMPLER_TEMPERATURE_MULT );
        MaybeReadRParam( cascadeParams.turbineParams.particleSamplingPeriod, params, BIMAP_PARAM_SAMPLE_PERIOD_BICLUSTERING );

        MaybeReadRParam( precomputedDataParams.objectFreqThreshold, params, PRECOMP_PARAM_OBJECT_FREQ_THRESHOLD );
        MaybeReadRParam( precomputedDataParams.probeFreqThreshold, params, PRECOMP_PARAM_PROBE_FREQ_THRESHOLD );

        MaybeReadRParam( samplerParams.objectMembershipRate, params, BIMAP_PARAM_SAMPLE_RATE_OBJECT_MEMBERSHIP );
        MaybeReadRParam( samplerParams.objectsSplitMergeRate, params, BIMAP_PARAM_SAMPLE_RATE_OBJECTS_SPLIT_MERGE );
        MaybeReadRParam( samplerParams.probeMembershipRate, params, BIMAP_PARAM_SAMPLE_RATE_PROBE_MEMBERSHIP );
        MaybeReadRParam( samplerParams.probesSplitMergeRate, params, BIMAP_PARAM_SAMPLE_RATE_PROBES_SPLIT_MERGE );
        MaybeReadRParam( samplerParams.blockFlipRate, params, BIMAP_PARAM_SAMPLE_RATE_BLOCK_FLIP );
        MaybeReadRParam( samplerParams.signalRate, params, BIMAP_PARAM_SAMPLE_RATE_SIGNAL );
        MaybeReadRParam( samplerParams.blockResamples, params, BIMAP_PARAM_BLOCK_RESAMPLES );
        MaybeReadRParam( collectorParams.priorsStoragePeriod, params, BIMAP_PARAM_SAMPLE_PERIOD_PRIORS );

        MaybeReadRParam( signalParams.sequenceLengthFactor, params, BIMAP_PARAM_SIGNAL_SEQUENCE_LENGTH_FACTOR );
        MaybeReadRParam( signalParams.scShape, params, BIMAP_PARAM_SIGNAL_SHAPE );

        MaybeReadRParam( priors.objectClustering.concentration, params, BIMAP_PARAM_PRIOR_OBJECTS_CLUSTERING_CONCENTRATION );
        MaybeReadRParam( priors.objectClustering.discount, params, BIMAP_PARAM_PRIOR_OBJECTS_CLUSTERING_DISCOUNT );

        MaybeReadRParam( priors.probeClustering.concentration, params, BIMAP_PARAM_PRIOR_PROBES_CLUSTERING_CONCENTRATION );
        MaybeReadRParam( priors.probeClustering.discount, params, BIMAP_PARAM_PRIOR_PROBES_CLUSTERING_DISCOUNT );
        MaybeReadRParam( priors.cellEnablementProb, params, BIMAP_PARAM_PRIOR_BLOCK_ENABLED );

        MaybeReadRParam( priors.noise.successes, params, BIMAP_PARAM_PRIOR_TRUE_MISSES );
        MaybeReadRParam( priors.noise.failures, params, BIMAP_PARAM_PRIOR_FALSE_HITS );

        MaybeReadRParam( hyperpriors.signalHyperprior.meanMean, params, 
                            BIMAP_PARAM_HPRIOR_BASELINE_SIGNAL );
        MaybeReadRParam( hyperpriors.signalHyperprior.meanVarScale, params, 
                            BIMAP_PARAM_HPRIOR_BASELINE_SCALE );
        MaybeReadRParam( hyperpriors.signalHyperprior.varDistrib.shape, params, 
                            BIMAP_PARAM_HPRIOR_SIGNAL_VAR_SHAPE );
        MaybeReadRParam( hyperpriors.signalHyperprior.varDistrib.scale, params, 
                            BIMAP_PARAM_HPRIOR_SIGNAL_VAR_SCALE );

        MaybeReadRParam( particleLogFilename, params, BIMAP_PARAM_LOG_PARTICLES_FILE );
        MaybeReadRParam( eeJumpLogFilename, params, BIMAP_PARAM_LOG_EEJUMPS_FILE );

        MaybeReadRParam( objectsComponentsThreshold, params, BIMAP_PARAM_OBJECTS_COMPONENTS_THRESHOLD );
        MaybeReadRParam( probesComponentsThreshold, params, BIMAP_PARAM_PROBES_COMPONENTS_THRESHOLD );

        MaybeReadRParam( objectsClotThreshold, params, BIMAP_PARAM_OBJECTS_CLOT_THRESHOLD );

        MaybeReadRParam( stableObjectsClustersThreshold, params, BIMAP_PARAM_STABLE_OBJECTS_CLUSTERS_THRESHOLD );
        MaybeReadRParam( stableProbesClustersThreshold, params, BIMAP_PARAM_STABLE_PROBES_CLUSTERS_THRESHOLD );

        priors.updateCachedDistributions();

        objects_clu_code_map_t  objsCluCode;
        try {
            SEXP  objsClusExp = PROTECT( params[ BIMAP_PARAM_INITIAL_OBJECTS_PARTITION ] );
            if ( !Rf_isNull( objsClusExp ) ) {
                objsCluCode = ReadObjectsClusters( iniClus, objsClusExp, data );
            }
            UNPROTECT( 1 );
        }
        catch ( Rcpp::index_out_of_bounds& e ) {
        }
        if ( objsCluCode.empty() ) {
            // put each object into separate cluster
            Rprintf( "Putting each object to separate cluster...\n" );
            for ( object_index_t objIx = 0; objIx < iniClus.objectsCount(); objIx++ ) {
                iniClus.addObjectCluster( objIx );
            }
            iniClus.cleanupClusters();
        }
        probes_clu_code_map_t  probesCluCode;
        try {
            SEXP  probesClusExp = PROTECT( params[ BIMAP_PARAM_INITIAL_PROBES_PARTITION ] );
            if ( !Rf_isNull( probesClusExp ) ) {
                probesCluCode = ReadProbesClusters( iniClus, probesClusExp, data );
            }
            UNPROTECT( 1 );
        }
        catch ( Rcpp::index_out_of_bounds& e ) {
        }
        if ( probesCluCode.empty() ) {
            // put each probe into separate cluster
            Rprintf( "Putting each probe to separate cluster...\n" );
            for ( probe_index_t probeIx = 0; probeIx < iniClus.probesCount(); probeIx++ ) {
                iniClus.addProbeCluster( probeIx );
            }
            iniClus.cleanupClusters();
        }
        try {
            SEXP  blocksExp = PROTECT( params[ BIMAP_PARAM_INITIAL_BLOCKS ] );
            if ( !Rf_isNull( blocksExp ) ) {
                if ( objsCluCode.empty() ) {
                    throw std::runtime_error( "Objects partition not specified, cannot read blocks" );
                }
                if ( probesCluCode.empty() ) {
                    throw std::runtime_error( "Probes partition not specified, cannot read blocks" );
                }
                ReadChessboardBiclusteringBlocks( iniClus, objsCluCode, probesCluCode, blocksExp );
            }
            UNPROTECT( 1 );
        }
        catch ( Rcpp::index_out_of_bounds& e ) {
        }
    }

    Rprintf( "Initializing sampler...\n" );
    ChessboardBiclusteringsIndexing  chessboardBiclusteringsIndexing( data.probesCount() );
    REESamplerExecutionMonitor execMonitor( 1, particleLogFilename.empty() ? NULL : particleLogFilename.c_str(), 
                                               eeJumpLogFilename.empty() ? NULL : eeJumpLogFilename.c_str() );
    Rprintf( "Evaluating preliminary signals...\n" );
    PrecomputedData precomputed( data, precomputedDataParams, signalParams );

    BIMAPSamplerHelper  helper( precomputed, hyperpriors, priors,
                                 samplerParams, 0, &execMonitor );

    Rprintf( "Bayesian Biclustering iterations in progress...\n" );
    BIMAPWalk res = BIMAPSampler_run( helper, chessboardBiclusteringsIndexing,
                                        samplerParams, cascadeParams,
                                        collectorParams,
                                        iniClus, &execMonitor );

    Rprintf( "Bayesian Biclustering finished\n" );
    boost::scoped_ptr<ChessboardBiclusteringsPDFEval> pdfAdjust;
    boost::scoped_ptr<IndexedPartitionsCollection<ObjectsCluster> > objPtnColl;
    boost::scoped_ptr<IndexedPartitionsCollection<ProbesCluster> > probePtnColl;
    if ( !R_IsNA( objectsComponentsThreshold ) || !R_IsNA( probesComponentsThreshold ) ) {
        if ( R_IsNA( objectsComponentsThreshold ) ) objectsComponentsThreshold = 0;
        if ( R_IsNA( probesComponentsThreshold ) ) probesComponentsThreshold = 0;

        LOG_INFO( "Detecting independent objects components...\n" );
        objPtnColl.reset( new IndexedPartitionsCollection<ObjectsCluster>( res ) );
        objPtnColl->init( helper.rndNumGen, res.stepsCount() * objectsComponentsThreshold,
                          res.stepsCount() * objectsClotThreshold );

        LOG_INFO( "Detecting independent probes components...\n" );
        probePtnColl.reset( new IndexedPartitionsCollection<ProbesCluster>( res ) );
        probePtnColl->init( nullptr, res.stepsCount() * probesComponentsThreshold,
                           std::numeric_limits<prob_t>::quiet_NaN() );

        LOG_INFO( "Calculating blocks statistics...\n" );
        pdfAdjust.reset( new ChessboardBiclusteringsPDFEval( res,
                                           *objPtnColl, *probePtnColl ) );
    }
    objects_label_map_type  objects = data.objectsLabelMap();
    probes_label_map_type   probes = data.probesLabelMap();
    if ( !walkFilename.empty() ) {
        Rprintf( "Saving BIMAP walk to '%s'...\n", walkFilename.c_str() );
        BIMAPResultsSave( walkFilename.c_str(), objects, probes, res, pdfAdjust.get() );
        Rprintf( "Saving done\n" );
    }
    boost::scoped_ptr<BlocksScoring> blocksScoring;
    if ( !R_IsNA( stableObjectsClustersThreshold )
         && !R_IsNA( stableProbesClustersThreshold )
         && pdfAdjust
    ){
        blocksScoring.reset( BIMAPBlocksScoring( *pdfAdjust, *objPtnColl, *probePtnColl,
                                                 stableObjectsClustersThreshold,
                                                 stableProbesClustersThreshold ) );
    }
    if ( createWalkRObject ) {
        return ( ConvertBIMAPWalkToRObject( objects, probes, res, pdfAdjust.get(), blocksScoring.get() ) );
    }
    else {
        return ( R_NilValue );
    }

    END_RCPP
}

/**
 *  Loads BI-MAP MCMC walk file and calculates independent components.
 */
RcppExport SEXP BIMAPWalkLoad(
    SEXP filenameExp,                   /** @param[in] filenameExp walk filename, character string */
    SEXP objectsComponentsThresholdExp, /** @param[in] objectsComponentsThresholdExp threshold for defining independent objects' components, numeric */
    SEXP probesComponentsThresholdExp,  /** @param[in] probesComponentsThresholdExp  threshold for defining independent probes' components, numeric */
    SEXP objectsClotThresholdExp,       /** @param[in] objectsClotThresholdExp   threshold for defining "clots" of objects (highly co-occurring), numeric */
    SEXP stableObjectsClustersThresholdExp, /** @param[in] stableObjectsClustersThresholdExp threshold for defining stable (frequently included) objects' clusters, numeric */
    SEXP stableProbesClustersThresholdExp   /** @param[in] stableProbesClustersThresholdExp threshold for defining stable (frequently included) probes' components, numeric */
){
    BEGIN_RCPP

    std::string filename;
    {
        Rcpp::StringVector filenameVec( filenameExp );
        filename = filenameVec[0];
    }

    Rprintf( "Reading BIMAP walk file %s...\n", filename.c_str() );
    ChessboardBiclusteringsIndexing            indexing;
    objects_label_map_type  objects;
    probes_label_map_type   probes;
    BIMAPWalk* walk = nullptr;
    ChessboardBiclusteringsPDFEval* pdfAdjust = nullptr;
    BIMAPResultsLoad( filename.c_str(), indexing, objects, probes, &walk, &pdfAdjust );

    if ( !walk ) {
        Rprintf( "No BIMAP walk loaded...\n" );
        return ( R_NilValue );
    }
    Rprintf( "Checking...\n" );
    if ( !walk->check() ) {
        throw std::runtime_error( "Error detected in BIMAP walk, loaded from file" );
    }

    double objectsComponentsThreshold = Rcpp::as<double>( objectsComponentsThresholdExp ); 
    double probesComponentsThreshold = Rcpp::as<double>( probesComponentsThresholdExp ); 
    double stableObjectsClustersThreshold = Rcpp::as<double>( stableObjectsClustersThresholdExp );
    double stableProbesClustersThreshold = Rcpp::as<double>( stableProbesClustersThresholdExp );
    double objectsClotThreshold = Rcpp::as<double>( objectsClotThresholdExp ); 
    boost::scoped_ptr<IndexedPartitionsCollection<ObjectsCluster> > objPtnColl;
    boost::scoped_ptr<IndexedPartitionsCollection<ProbesCluster> > probePtnColl;
    boost::scoped_ptr<BlocksScoring> blocksScoring;
    if ( !R_IsNA( objectsComponentsThreshold ) || !R_IsNA( probesComponentsThreshold ) ) {
        if ( pdfAdjust ) {
            Rprintf( "Deleting old PDF adjustment...\n" );
            delete pdfAdjust;
        }
        if ( R_IsNA( objectsComponentsThreshold ) ) objectsComponentsThreshold = 0;
        if ( R_IsNA( probesComponentsThreshold ) ) probesComponentsThreshold = 0;
        if ( R_IsNA( objectsClotThreshold ) ) objectsClotThreshold = std::numeric_limits<prob_t>::quiet_NaN();
        cemm::math::RngPtr rndNumGen( gsl_rng_alloc( gsl_rng_default ) );

        if ( R_IsNA( objectsComponentsThreshold ) ) objectsComponentsThreshold = 0;
        if ( R_IsNA( probesComponentsThreshold ) ) probesComponentsThreshold = 0;

        LOG_INFO( "Detecting independent objects components...\n" );
        objPtnColl.reset( new IndexedPartitionsCollection<ObjectsCluster>( *walk ) );
        objPtnColl->init( rndNumGen, walk->stepsCount() * objectsComponentsThreshold,
                          walk->stepsCount() * objectsClotThreshold );

        LOG_INFO( "Detecting independent probes components...\n" );
        probePtnColl.reset( new IndexedPartitionsCollection<ProbesCluster>( *walk ) );
        probePtnColl->init( nullptr, walk->stepsCount() * probesComponentsThreshold,
                           std::numeric_limits<prob_t>::quiet_NaN() );

        Rprintf( "Calculating PDF adjustment with independent components threshold o=%f, s=%f", 
                 objectsComponentsThreshold, probesComponentsThreshold );
        if ( std::isfinite( objectsClotThreshold ) ) {
            Rprintf( "\nand generating object clusters intersections with clot threshold o=%f", 
                    objectsClotThreshold );
        }
        Rprintf( "...\n" );
        pdfAdjust = new ChessboardBiclusteringsPDFEval( *walk,
                                           *objPtnColl, *probePtnColl );
    }
    if ( objPtnColl && probePtnColl &&
         !R_IsNA( objectsComponentsThreshold ) && !R_IsNA( probesComponentsThreshold )
    ){
        blocksScoring.reset( BIMAPBlocksScoring( *pdfAdjust, *objPtnColl, *probePtnColl,
                                                 stableObjectsClustersThreshold, stableProbesClustersThreshold ) );
    }

    //Rprintf( "Exporting walk to R...\n" );
    SEXP res = ConvertBIMAPWalkToRObject( objects, probes, *walk,
                                          pdfAdjust, blocksScoring.get() );
    delete walk;
    if ( pdfAdjust ) delete pdfAdjust;

    return ( res );

    END_RCPP
}

#if 0
RcppExport SEXP OPADataLoad(
    SEXP filenameExp
){
    BEGIN_RCPP

    std::string filename;
    {
        Rcpp::StringVector filenameVec( filenameExp );
        filename = filenameVec[0];
    }

    Rprintf( "Reading OPAData file %s...\n", filename.c_str() );
    OPAData data = OPAData::load( filename.c_str() );

    Rprintf( "Exporting OPAData to R...\n" );
    return ( ConvertBIMAPWalkToRObject( walk, objects, probes ) );

    END_RCPP;
}
#endif

/**
    Saves OPA data to disk in internal format.
 */
RcppExport SEXP OPADataSave(
    SEXP    filenameExp,            /**< @param[in] name of file to save OPA data to */
    SEXP    proteinsDataExp,        /**< @param[in] proteins data frame: 1 - protein label, 2 - sequence length */
    SEXP    samplesDataExp,         /**< @param[in] samples data frame: 1 - sample label, 2 - bait protein label */
    SEXP    msRunsDataExp,          /**< @param[in] MS-runs data frame: 1 - MS run label, 2 - examined sample label */
    SEXP    measurementsDataExp,    /**< @param[in] MS measurements data frame: 1 - MS run label, 2 - protein label, 3 - TSC */
    SEXP    paramsExp               /**< @param[in] method params */
){
    BEGIN_RCPP

    OPAData  data = ReadMassSpecData( proteinsDataExp, samplesDataExp, msRunsDataExp, measurementsDataExp );

    std::string filename;
    {
        Rcpp::StringVector filenameVec( filenameExp );
        filename = filenameVec[0];
    }
    Rprintf( "Saving OPA data to '%s'...\n", filename.c_str() );
    data.save( filename.c_str() );
    Rprintf( "Data saved\n" );

    END_RCPP
}

Rcpp::NumericMatrix CreateDistanceMatrix(
    const symmetric_array2d<double>&    matrix,
    const Rcpp::CharacterVector&        elmNames
){
    if ( (int)matrix.size() != elmNames.size() ) {
        THROW_EXCEPTION( std::invalid_argument, "Names length and square matrix dimensions do not match" );
    }
    Rcpp::NumericMatrix rmtx( matrix.size(), matrix.size() );
    for ( size_t i = 0; i < matrix.size(); i++ ) {
        for ( size_t j = 0; j < matrix.size(); j++ ) {
            rmtx( i, j ) = matrix( i, j );
        }
    }
    rmtx.attr("dimnames") = Rcpp::List::create( elmNames, elmNames );
    return ( rmtx );
}

/**
    Calculates distances between objects and probes.
 */
RcppExport SEXP CalcDistances(
    SEXP    proteinsDataExp,        /**< @param[in] proteins data frame: 1 - protein label, 2 - sequence length */
    SEXP    samplesDataExp,         /**< @param[in] samples data frame: 1 - sample label, 2 - bait protein label */
    SEXP    msRunsDataExp,          /**< @param[in] MS-runs data frame: 1 - MS run label, 2 - examined sample label */
    SEXP    measurementsDataExp,    /**< @param[in] MS measurements data frame: 1 - MS run label, 2 - protein label, 3 - TSC */
    SEXP    paramsExp               /**< @param[in] method params */
){
    BEGIN_RCPP

    OPAData  data = ReadMassSpecData( proteinsDataExp, samplesDataExp, msRunsDataExp, measurementsDataExp );

    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;

    {
        Rprintf( "Reading parameters...\n" );
        RcppParamVector      params( paramsExp );

        MaybeReadRParam( signalParams.sequenceLengthFactor, params, BIMAP_PARAM_SIGNAL_SEQUENCE_LENGTH_FACTOR );
        MaybeReadRParam( signalParams.scShape, params, BIMAP_PARAM_SIGNAL_SHAPE );

        MaybeReadRParam( precomputedDataParams.objectFreqThreshold, params, PRECOMP_PARAM_OBJECT_FREQ_THRESHOLD );
        MaybeReadRParam( precomputedDataParams.probeFreqThreshold, params, PRECOMP_PARAM_PROBE_FREQ_THRESHOLD );
    }

    Rprintf( "Preparing...\n" );
    PrecomputedData::sc_cache_type scDistribCache( signalParams, 0, data.maxMeasurement().sc, signalParams.scRateStep );
    CellSignalLLHMaximizer signalLLH( data, signalParams, scDistribCache );
    Rprintf( "Calculating objects signal LLH...\n" );
    symmetric_array2d<double> objectCoSignals = signalLLH.evalObjectCoSignalDistances();
    Rprintf( "Calculating probes signal LLH...\n" );
    symmetric_array2d<double> probeCoSignals = signalLLH.evalProbeCoSignalDistances();
    Rprintf( "Calculating objects co-occurrence...\n" );
    std::vector<PrecomputedData::observations_mask_t> objObs = PrecomputedData::ObjectsObservations( data );
    std::vector<PrecomputedData::observations_mask_t> probeObs = PrecomputedData::ProbesObservations( data );
    PrecomputedData::observations_mask_t spProbes =
        CoOccurenceDistanceMatrix<log_prob_t>::SpecificHits( objObs, precomputedDataParams.probeFreqThreshold );
    CoOccurenceDistanceMatrix<log_prob_t> objectCoOccurrence(
        PrecomputedData::ObjectsObservations( data ),
        precomputedDataParams.probeFreqThreshold, true );
    double probesBias = CoOccurenceDistanceMatrix<log_prob_t>::ObservationsBias( objObs, spProbes );
    symmetric_array2d<double> objCohits = CoOccurenceDistanceMatrix<log_prob_t>::CoOccurenceDistances( objObs, spProbes, probesBias );
    Rprintf( "Calculating probes co-occurrence...\n" );
    PrecomputedData::observations_mask_t spObjs =
        CoOccurenceDistanceMatrix<log_prob_t>::SpecificHits( probeObs, precomputedDataParams.objectFreqThreshold );
    CoOccurenceDistanceMatrix<log_prob_t> probeCoOccurrence(
        PrecomputedData::ProbesObservations( data ),
        precomputedDataParams.objectFreqThreshold, true );
    double objsBias = CoOccurenceDistanceMatrix<log_prob_t>::ObservationsBias( probeObs, spObjs );
    symmetric_array2d<double> probesCohits = CoOccurenceDistanceMatrix<log_prob_t>::CoOccurenceDistances( probeObs, spObjs, objsBias );
    Rprintf( "Calculation finished\n" );

    objects_label_map_type  objects = data.objectsLabelMap();
    probes_label_map_type   probes = data.probesLabelMap();

    Rcpp::CharacterVector objVec( objects.size() );
    for ( objects_label_map_type::const_iterator oit = objects.begin(); oit != objects.end(); ++oit ) {
        objVec[ oit->first ] = oit->second;
    }
    Rcpp::CharacterVector probeVec( probes.size() );
    for ( probes_label_map_type::const_iterator sit = probes.begin(); sit != probes.end(); ++sit ) {
        probeVec[ sit->first ] = sit->second;
    }
    Rcpp::CharacterVector spObjVec( spObjs.count() );
    size_t i = 0;
    foreach_bit( object_index_t, objIx, spObjs ) {
        spObjVec[ i++ ] = objects.find( objIx )->second;
    }
    Rcpp::CharacterVector spProbeVec( spProbes.count() );
    i = 0;
    foreach_bit( probe_index_t, probeIx, spProbes ) {
        spProbeVec[ i++ ] = probes.find( probeIx )->second;
    }
    return ( Rcpp::List::create(
        Rcpp::Named( "objects.signal" ) = CreateDistanceMatrix( objectCoSignals, objVec ),
        Rcpp::Named( "probes.signal" ) = CreateDistanceMatrix( probeCoSignals, probeVec ),

        Rcpp::Named( "objects.bias" ) = objsBias,
        Rcpp::Named( "objects.specific" ) = spObjVec,
        Rcpp::Named( "probes.cohits" ) = CreateDistanceMatrix( probesCohits, probeVec ),

        Rcpp::Named( "objects.cohits" ) = CreateDistanceMatrix( objCohits, objVec ),
        Rcpp::Named( "probes.specific" ) = spProbeVec,
        Rcpp::Named( "probes.bias" ) = objectCoOccurrence.observationsBias()
    ) );

    END_RCPP
}

} }
