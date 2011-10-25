#pragma once

#include "BasicTypedefs.h"

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/indexed_by.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_array.hpp>

#include <boost/serialization/version.hpp>
#include <boost/serialization/nvp.hpp>

#include "ChessboardBiclustering.h"
#include "ChessboardBiclusteringsIndexing.h"

/**
    Sampling walk in the space of chessboard biclusterings
    and their prior parameters.
 */
class BIMAPWalk {
private:
    friend class boost::serialization::access;

    typedef double  step_time_t;

    struct BIMAPStep {
        step_time_t                 time;
        size_t                      turbineIx;
        log_prob_t                  llh;
        log_prob_t                  lpp;
        ChessboardBiclusteringIndexed      clustering;

        BIMAPStep( step_time_t time, size_t turbineIx, const ChessboardBiclustering& clustering, log_prob_t llh, log_prob_t lpp, ChessboardBiclusteringsIndexing& indexing )
        : time( time ), turbineIx( turbineIx ), llh( llh ), lpp( lpp ), clustering( indexing.index( clustering ) )
        {}

        BIMAPStep( step_time_t time, size_t turbineIx, const ChessboardBiclusteringIndexed& clustering, log_prob_t llh, log_prob_t lpp )
        : time( time ), turbineIx( turbineIx ), llh( llh ), lpp( lpp ), clustering( clustering )
        {}
    };

    typedef boost::multi_index::multi_index_container<
        BIMAPStep,
        boost::multi_index::indexed_by<
            // by order of appending
            boost::multi_index::sequenced<>,
            // by sample time
            boost::multi_index::ordered_non_unique<boost::multi_index::member<BIMAPStep, const step_time_t, &BIMAPStep::time> >
        >
    > StepsIndexing;

    struct BIMAPPriorParamsStep {
        friend class boost::serialization::access;

        step_time_t                     time;
        size_t                          turbineIx;
        ChessboardBiclusteringDerivedPriors    priors;

        BIMAPPriorParamsStep()
        : time( 0 ), turbineIx( 0 )
        {}

        BIMAPPriorParamsStep( step_time_t time, size_t turbineIx, const ChessboardBiclusteringDerivedPriors& priors )
        : time( time ), turbineIx( turbineIx ), priors( priors )
        {}

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP( time );
            ar & BOOST_SERIALIZATION_NVP( turbineIx );
            ar & BOOST_SERIALIZATION_NVP( priors );
        }
    };

    typedef boost::multi_index::multi_index_container<
        BIMAPPriorParamsStep,
        boost::multi_index::indexed_by<
            // by order of appending
            boost::multi_index::sequenced<>,
            // by time
            boost::multi_index::ordered_unique<boost::multi_index::member<BIMAPPriorParamsStep, const step_time_t, &BIMAPPriorParamsStep::time> >
        >
    > PriorParamsStepsIndexing;

    ChessboardBiclusteringsIndexing&   _crossClusteringsIndexing;        /** walker's indexer of biclusterings */
    StepsIndexing               _steps;
    PriorParamsStepsIndexing    _priorParamsSteps;

    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        boost::serialization::collection_size_type stepsCount( _steps.get<0>().size() );
        ar << BOOST_SERIALIZATION_NVP( stepsCount );
        for ( StepsIndexing::const_iterator stepit = _steps.get<0>().begin(); stepit != _steps.get<0>().end(); ++stepit ) {
            ChessboardBiclusteringIndexed& cclus = const_cast<ChessboardBiclusteringIndexed&>( stepit->clustering );
            ar << boost::serialization::make_nvp( "time", stepit->time );
            ar << boost::serialization::make_nvp( "turbineIx", stepit->turbineIx );
            ar << boost::serialization::make_nvp( "llh", stepit->llh );
            ar << boost::serialization::make_nvp( "lpp", stepit->lpp );
            default_serial_t serial = cclus.serial();
            ar << boost::serialization::make_nvp( "crossClusteringSerial", serial );
            ar << boost::serialization::make_nvp( "clusteringData", cclus );
            LOG_DEBUG2( stepit->iteration << ": " << cclus._objectClusterData.size() );
            ar << boost::serialization::make_nvp( "crossClustersData", cclus._crossClustersData );
            ar << boost::serialization::make_nvp( "objectsData", cclus._objectsData );
        }
        ar << boost::serialization::make_nvp( "priorParamSteps", _priorParamsSteps );
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        boost::serialization::collection_size_type stepsCount;
        ar >> BOOST_SERIALIZATION_NVP( stepsCount );

        _steps.clear();
        while ( stepsCount-- ) {
            step_time_t time;
            ar >> BOOST_SERIALIZATION_NVP( time );
            size_t      turbineIx;
            ar >> BOOST_SERIALIZATION_NVP( turbineIx );
            log_prob_t      llh;
            ar >> BOOST_SERIALIZATION_NVP( llh );
            log_prob_t      lpp;
            ar >> BOOST_SERIALIZATION_NVP( lpp );
            default_serial_t crossClusteringSerial;
            ar >> BOOST_SERIALIZATION_NVP( crossClusteringSerial );
            ChessboardBiclusteringsIndexing::const_serial_iterator pScaffoldIt = _crossClusteringsIndexing.iterator_to( crossClusteringSerial );
            if ( pScaffoldIt == _crossClusteringsIndexing.serial_not_found() ) {
                THROW_RUNTIME_ERROR( "Chessboard biclustering #" << crossClusteringSerial << " not found the index" );
            }
            ChessboardBiclusteringData clusteringData;
            ar >> BOOST_SERIALIZATION_NVP( clusteringData );
            ChessboardBiclusteringIndexed::cross_cluster_data_map_type crossClustersData;
            ar >> BOOST_SERIALIZATION_NVP( crossClustersData );
            multiple_map_t      objectsData;
            ar >> BOOST_SERIALIZATION_NVP( objectsData );

            _steps.get<0>().push_back( BIMAPStep( time, turbineIx, ChessboardBiclusteringIndexed( *pScaffoldIt, crossClustersData, objectsData, clusteringData ), llh, lpp ) );
        }
        ar >> boost::serialization::make_nvp( "priorParamSteps", _priorParamsSteps );
    }

protected:
    friend class BIMAPSampler;

public:
    BIMAPWalk( ChessboardBiclusteringsIndexing& crossClusteringsIndexing );

    typedef StepsIndexing::const_iterator const_step_iterator;
    typedef PriorParamsStepsIndexing::const_iterator const_priors_step_iterator;

    void step( step_time_t time, size_t turbineIx, const ChessboardBiclustering& clu, log_prob_t llh, log_prob_t lpp ) {
        _steps.get<0>().push_back( BIMAPStep( time, turbineIx, clu, llh, lpp, _crossClusteringsIndexing ) );
    }
    void step( step_time_t time, size_t turbineIx, const ChessboardBiclusteringDerivedPriors& priors ) {
        _priorParamsSteps.get<0>().push_back( BIMAPPriorParamsStep( time, turbineIx, priors ) );
    }

    const ChessboardBiclusteringsIndexing& indexing() const {
        return ( _crossClusteringsIndexing );
    }
    ChessboardBiclusteringsIndexing& indexing() {
        return ( _crossClusteringsIndexing );
    }

    size_t stepsCount() const {
        return ( _steps.size() );
    }
    const_step_iterator stepsBegin() const {
        return ( _steps.begin() );
    }
    const_step_iterator stepsEnd() const {
        return ( _steps.end() );
    }

    size_t priorParamsStepsCount() const {
        return ( _priorParamsSteps.size() );
    }
    const_priors_step_iterator priorParamsStepsBegin() const {
        return ( _priorParamsSteps.begin() );
    }
    const_priors_step_iterator priorParamsStepsEnd() const {
        return ( _priorParamsSteps.end() );
    }

    bool check() const;

    size_t filterSteps( size_t minScaffoldCounts = 1, size_t minObjectsPtnCounts = 2, size_t minProbesPtnCounts = 2 );
};

BOOST_CLASS_VERSION(BIMAPWalk, 1);
BOOST_CLASS_VERSION(BIMAPWalk::BIMAPPriorParamsStep, 1);
