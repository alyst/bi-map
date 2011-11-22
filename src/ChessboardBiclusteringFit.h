#pragma once

#include "BasicTypedefs.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

#include "ChessboardBiclustering.h"
#include "PrecomputedData.h"
#include "DataSignalNoiseCache.h"
#include "ChessboardBiclusteringMetrics.h"

#include "dynamic_bitset_utils.h"

#include "mcmc/MetropolisHastingsStep.h"
#include "array2d.h"
#include "symmetric_array2d.h"
#include "statically_tracked.h"

#define LN_PROB_ERROR_TOL 1E-5

#define REPORT_PROB_ERROR( msg ) LOG_WARN( msg )
#define REPORT_PROB_ERROR_IF( cond, msg ) LOG_WARN_IF( cond, msg )
//#define REPORT_PROB_ERROR( msg ) THROW_RUNTIME_ERROR( msg )
//#define REPORT_PROB_ERROR_IF( cond, msg ) if ( cond ) { THROW_RUNTIME_ERROR( msg ) }

class StaticChessboardBiclustering;

/**
    Calculation with caching for chessboard biclustering.
    Evaluates likelihoods, priors etc.
 */
class ChessboardBiclusteringFit: public ChessboardBiclustering {
public:
    typedef unsigned char count_type;
    typedef array2d<count_type> block_counts_matrix_type;
    typedef array2d<log_prob_t> lnprob_matrix_type;

protected:
    typedef array2d<log_prob_t> block_lnprob_matrix_type;
    typedef symmetric_array2d<log_prob_t> lnprob_symmatrix_type;
    typedef std::vector<log_prob_t> lnprob_vector_type;

    const PrecomputedData&          _precomputed;
    const ChessboardBiclusteringPriors&    _priors;

    mutable boost::shared_ptr<DataSignalNoiseCache>  _signalNoiseCache;

    typedef std::vector<LLHMetrics> llhmetrics_vector_type;

    mutable StatsMetrics        _metrics;   /** statistical metrics */

    mutable llhmetrics_vector_type  _objClusterLLH;
    mutable llhmetrics_vector_type  _probeClusterLLH;

    mutable lnprob_vector_type  _objClusterLPP;

    mutable block_lnprob_matrix_type _blockLLH;
    mutable block_lnprob_matrix_type _blockLPP;

    mutable bool _blockIsSignalLLHValid;
    mutable block_lnprob_matrix_type _blockIsSignalLLH;
    mutable bool _blockIsNoiseLLHValid;
    mutable block_lnprob_matrix_type _blockIsNoiseLLH;

    block_counts_matrix_type _blockToSample; /** how many each block have to be sampled */

    void updateLPP() const;
    void updateBlocksProbesLLH() const;
    void updateClustersLLH() const;
    void updateMetrics() const;

    void resetAllCaches() const;
    void resetObjectsClusterCache( object_clundex_t objCluIx, bool contentsChanged ) const;
    void resetProbesClusterCache( probe_clundex_t probeCluIx ) const;
    void resetCachesAfterLoading();

    friend class boost::serialization::access;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> boost::serialization::make_nvp( "chessboardBiclustering", boost::serialization::base_object<ChessboardBiclustering>( *this ) );
        ar >> boost::serialization::make_nvp( "blockToSample", _blockToSample );
        ar >> boost::serialization::make_nvp( "blockLLH", _blockLLH );
        ar >> boost::serialization::make_nvp( "blockLPP", _blockLPP );
        ar >> boost::serialization::make_nvp( "blockNoiseLLHValid", _blockIsNoiseLLHValid );
        if ( _blockIsNoiseLLHValid ) {
            ar >> boost::serialization::make_nvp( "blockNoiseLLH", _blockIsNoiseLLH );
        } else {
            _blockIsNoiseLLH.reset( objectsClusters().size(), probesClusters().size(), unset() );
        }
        ar >> boost::serialization::make_nvp( "blockSignalLLHValid", _blockIsSignalLLHValid );
        if ( _blockIsSignalLLHValid ) {
            ar >> boost::serialization::make_nvp( "blockSignalLLH", _blockIsSignalLLH );
        } else {
            _blockIsSignalLLH.reset( objectsClusters().size(), probesClusters().size(), unset() );
        }
        ar >> boost::serialization::make_nvp( "metrics", _metrics );
        ar >> boost::serialization::make_nvp( "signalNoiseCache", _signalNoiseCache );
#ifdef _DEBUG
        // debugging version: check that serializing of llh/lpp works by recalculating after cache reset
        StatsMetrics prevMetrics = _metrics;
        _metrics.unset( true, true );
        StatsMetrics curMetric = _metrics;
        REPORT_PROB_ERROR_IF( std::abs( curMetric.lpp() - prevMetrics.lpp() ) > LN_PROB_ERROR_TOL,
                              "Incorrect LPP serialization: real=" << curMetric.lpp()
                              << " serialized=" << prevMetrics.lpp() );
        REPORT_PROB_ERROR_IF( std::abs( curMetric.llh( LLHWeights() ) - prevMetrics.llh( LLHWeights() ) ) > LN_PROB_ERROR_TOL,
                              "Incorrect LLH serialization: real=" << curMetric.llh( LLHWeights() )
                              << " serialized=" << prevMetrics.llh( LLHWeights() ) );
#else
        resetCachesAfterLoading();
#endif
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << boost::serialization::make_nvp( "chessboardBiclustering", boost::serialization::base_object<ChessboardBiclustering>( *this ) );
        ar << boost::serialization::make_nvp( "blockToSample", _blockToSample );
        ar << boost::serialization::make_nvp( "blockLLH", _blockLLH );
        ar << boost::serialization::make_nvp( "blockLPP", _blockLPP );
        ar << boost::serialization::make_nvp( "blockNoiseLLHValid", _blockIsNoiseLLHValid );
        if ( _blockIsNoiseLLHValid ) {
            ar << boost::serialization::make_nvp( "blockNoiseLLH", _blockIsNoiseLLH );
        }
        ar << boost::serialization::make_nvp( "blockSignalLLHValid", _blockIsSignalLLHValid );
        if ( _blockIsSignalLLHValid ) {
            ar << boost::serialization::make_nvp( "blockSignalLLH", _blockIsSignalLLH );
        }
        ar << boost::serialization::make_nvp( "metrics", _metrics );
        ar << boost::serialization::make_nvp( "signalNoiseCache", _signalNoiseCache );
    }

public:
    typedef ChessboardBiclustering::const_block_iterator const_block_iterator;

    ChessboardBiclusteringFit( const PrecomputedData& precomputed,
                        const ChessboardBiclusteringPriors& priors );
    ChessboardBiclusteringFit( const PrecomputedData& precomputed,
                        const ChessboardBiclusteringPriors& priors,
                        const ChessboardBiclustering& clus );

    ChessboardBiclusteringFit& operator=( const ChessboardBiclustering& that );
    ChessboardBiclusteringFit& operator=( const ChessboardBiclusteringFit& that );
    ChessboardBiclusteringFit& operator=( const StaticChessboardBiclustering& that );

    const PrecomputedData& precomputed() const {
        return ( _precomputed );
    }

    const OPAData& data() const {
        return ( _precomputed.data() );
    }

    std::set<object_clundex_t> boundObjectsClusters( const probe_bitset_t& probes ) const;
    std::set<object_clundex_t> boundObjectsClusters( probe_clundex_t clusterIx ) const {
        return ( boundObjectsClusters( probesCluster( clusterIx ).items() ) );
    }

    std::set<probe_clundex_t> boundProbesClusters( const object_set_t& objects ) const;
    std::set<probe_clundex_t> boundProbesClusters( object_clundex_t clusterIx ) const {
        return ( boundProbesClusters( objectsCluster( clusterIx ).items() ) );
    }

    const CellSignalParams& signalParams() const {
        return ( _precomputed.signalParams() );
    }

    const ChessboardBiclusteringPriors& priors() const {
        return ( _priors );
    }

    const StatsMetrics metrics() const {
        if ( _metrics.is_unset() ) {
            updateMetrics();
        }
        return ( _metrics );
    }
    const log_prob_t llh( const LLHWeights& weights ) const {
        return ( metrics().llh( weights ) );
    }
    log_prob_t lpp() const {
        return ( metrics().lpp() );
    }
    log_prob_t totalLnP( const LLHWeights& weights ) const {
#ifdef _DEBUG
        // debugging version: check that caching of llh/lpp works by recalculating after cache reset
        StatsMetrics prevMetrics = metrics();
        resetAllCaches();
        StatsMetrics curMetrics = metrics();
        REPORT_PROB_ERROR_IF( std::abs( curMetrics.lpp() - prevMetrics.lpp() ) > LN_PROB_ERROR_TOL,
                              "Incorrect LPP caching: new=" << curMetrics.lpp()
                              << " old=" << prevMetrics.lpp() );
        REPORT_PROB_ERROR_IF( std::abs( curMetrics.llh( weights ) - prevMetrics.llh( weights ) ) > LN_PROB_ERROR_TOL,
                              "Incorrect LLH caching: new=" << curMetrics.llh( weights )
                              << " old=" << prevMetrics.llh( weights ) );
        return ( curMetrics.totalLnP( weights ) );
#else
        return ( metrics().totalLnP( weights ) );
#endif
    }

    log_prob_t blockLLH( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const {
        log_prob_t res = _blockLLH( objCluIx, probeCluIx );
        if ( is_unset( res ) ) {
            updateClustersLLH();
            res = _blockLLH( objCluIx, probeCluIx );
        }
        return ( res );
    }

    LLHMetrics objectsClusterLLH( object_clundex_t cluIx ) const {
        LLHMetrics& res = _objClusterLLH[ cluIx ];
        if ( res.is_unset() ) {
            updateClustersLLH();
        }
        return ( res );
    }

    LLHMetrics probesClusterLLH( probe_clundex_t cluIx ) const {
        LLHMetrics& res = _probeClusterLLH[ cluIx ];
        if ( res.is_unset() ) {
            updateClustersLLH();
        }
        return ( res );
    }

    log_prob_t blockIsNoiseLLH( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const {
        if ( !_blockIsNoiseLLHValid ) {
            updateBlocksProbesLLH();
        }
        return ( _blockIsNoiseLLH( objCluIx, probeCluIx ) );
    }
    const block_lnprob_matrix_type::section_type objectsClusterIsNoiseLLH( object_clundex_t cluIx ) const {
        if ( !_blockIsNoiseLLHValid ) {
            updateBlocksProbesLLH(); // check if whole section1 is up-to-date
        }
        return ( _blockIsNoiseLLH.section1( cluIx ) );
    }
    const block_lnprob_matrix_type::section_type probesClusterIsNoiseLLH( probe_clundex_t cluIx ) const {
        if ( !_blockIsNoiseLLHValid ) {
            updateBlocksProbesLLH(); // check if whole section2 is up-to-date
        }
        return ( _blockIsNoiseLLH.section2( cluIx ) );
    }

    log_prob_t blockIsSignalLLH( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const {
        if ( !_blockIsSignalLLHValid ) {
            updateBlocksProbesLLH();
        }
        return ( _blockIsSignalLLH( objCluIx, probeCluIx ) );
    }
    const block_lnprob_matrix_type::section_type objectsClusterIsSignalLLH( object_clundex_t cluIx ) const {
        if ( !_blockIsSignalLLHValid ) {
            updateBlocksProbesLLH(); // check if whole section1 is up-to-date
        }
        return ( _blockIsSignalLLH.section1( cluIx ) );
    }
    const block_lnprob_matrix_type::section_type probesClusterIsSignalLLH( probe_clundex_t cluIx ) const {
        if ( !_blockIsSignalLLHValid ) {
            updateBlocksProbesLLH(); // check if whole section2 is up-to-date
        }
        return ( _blockIsSignalLLH.section2( cluIx ) );
    }

    DataSignalNoiseCache& signalNoiseCache() const {
        if ( !_signalNoiseCache ) {
            _signalNoiseCache.reset( new DataSignalNoiseCache( _precomputed, _priors, *this ) );
        }
        _signalNoiseCache->update();
        return ( *_signalNoiseCache );
    }

    size_t blockToSample( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const {
        return ( (size_t)_blockToSample( objCluIx, probeCluIx ) );
    }
    const block_counts_matrix_type& blocksToSample() const {
        return ( _blockToSample );
    }
    void setBlockSamples( object_clundex_t objCluIx, probe_clundex_t probeCluIx, size_t counts ) {
        count_type& cnts = _blockToSample( objCluIx, probeCluIx );
        cnts = std::max( cnts, (count_type)counts );
    }
    void decBlockSamples( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) {
        count_type& cnts = _blockToSample( objCluIx, probeCluIx );
        if ( cnts > 0 ) cnts--;
    }
    void setObjectsClusterSamples( object_clundex_t objCluIx, size_t counts );
    void setProbesClusterSamples( probe_clundex_t probeCluIx, size_t counts );
    void setAllBlockSamples( size_t counts );

protected:
    virtual void beforeObjectsClusterRemoved( object_clundex_t cluIx ) const;
    virtual void beforeProbesClusterRemoved( probe_clundex_t cluIx ) const;
    virtual void afterObjectsClusterInserted( object_clundex_t cluIx ) const;
    virtual void afterProbesClusterInserted( probe_clundex_t cluIx ) const;
    virtual void afterObjectsClusterChanged( object_clundex_t cluIx ) const;
    virtual void afterProbesClusterChanged( probe_clundex_t cluIx ) const;
    virtual void afterObjectMultipleChanged( object_index_t objIx ) const;
    virtual void afterSignalChanged( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const;
    virtual void afterBlockFlipped( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const;
    virtual void afterSignalPriorChanged() const;
    virtual void afterNoiseParamsChanged() const;
};

BOOST_CLASS_IMPLEMENTATION( ChessboardBiclusteringFit, object_serializable )

/**
 *  Chessboard biclustering + LLH + LPP.
 *  Realization of static_particle_type concept.
 */
struct StaticChessboardBiclustering {
    typedef StatsMetrics pre_energy_type;
    typedef ChessboardBiclusteringFit::block_counts_matrix_type block_counts_matrix_type;

    boost::shared_ptr<ChessboardBiclustering>       clustering;
    boost::shared_ptr<block_counts_matrix_type>     blocksToSample; /** how many each block have to be sampled */
    StatsMetrics                                    metrics;

    StaticChessboardBiclustering()
    {}

    StaticChessboardBiclustering( const ChessboardBiclustering& clustering )
    : clustering( new ChessboardBiclustering( clustering ) )
    {
    }

    StaticChessboardBiclustering( const ChessboardBiclusteringFit& clustering )
    : clustering( new ChessboardBiclustering( clustering ) )
    , blocksToSample( new block_counts_matrix_type( clustering.blocksToSample() ) )
    , metrics( clustering.metrics() )
    {
    }

    const pre_energy_type& preEnergy() const {
        return ( metrics );
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( metrics );
        ar & BOOST_SERIALIZATION_NVP( clustering );
        ar & BOOST_SERIALIZATION_NVP( blocksToSample );
    }

    friend std::ostream& operator<<( std::ostream& out, const StaticChessboardBiclustering& cc )
    {
        if ( cc.clustering ) {
            return ( out << ( *cc.clustering ) );
        } else {
            return ( out << "(null)" );
        }
    }
};

BOOST_CLASS_IMPLEMENTATION( StaticChessboardBiclustering, object_serializable )
BOOST_CLASS_TRACKING( StaticChessboardBiclustering, track_never )

struct ChessboardBiclusteringEnergyEval {
    friend class boost::serialization::access;

    typedef StatsMetrics pre_energy_type;
    typedef log_prob_t energy_type;
    typedef StaticChessboardBiclustering particle_type;

    LLHWeights  weights;

    energy_type operator()( const pre_energy_type& preEnergy ) const {
        return ( -preEnergy.totalLnP( weights ) );
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( weights );
    }
};

BOOST_CLASS_IMPLEMENTATION( ChessboardBiclusteringEnergyEval, object_serializable )
BOOST_CLASS_TRACKING( ChessboardBiclusteringEnergyEval, track_never )

namespace boost { namespace serialization {

template<class Archive>
inline void save_construct_data(
    Archive & ar, const ChessboardBiclusteringFit* t, const unsigned int file_version
){
    // save data required to construct instance
    const PrecomputedData* pPrecomputed = &t->precomputed();
    statically_tracked<PrecomputedData> tPrecomputed( "precomputed", pPrecomputed );
    ar << boost::serialization::make_nvp( "precomputed", tPrecomputed );

    const ChessboardBiclusteringPriors* pPriors = &t->priors();
    statically_tracked<ChessboardBiclusteringPriors> tPriors( "priors", pPriors );
    ar << boost::serialization::make_nvp( "priors", tPriors );
}

template<class Archive>
inline void load_construct_data(
    Archive & ar, ChessboardBiclusteringFit* t, const unsigned int file_version
){
    const PrecomputedData* pPrecomputed;
    statically_tracked<PrecomputedData> tPrecomputed( "precomputed", pPrecomputed );
    ar >> boost::serialization::make_nvp( "precomputed", tPrecomputed );

    const ChessboardBiclusteringPriors* pPriors;
    statically_tracked<ChessboardBiclusteringPriors> tPriors( "priors", pPriors );
    ar >> boost::serialization::make_nvp( "priors", tPriors );
    // invoke inplace constructor to initialize instance of my_class
    ::new(t)ChessboardBiclusteringFit( *pPrecomputed, *pPriors );
}

}
}
