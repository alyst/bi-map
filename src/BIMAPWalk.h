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
#include "ChessboardBiclusteringMetrics.h"

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
        StatsMetrics                metrics;
        ChessboardBiclusteringIndexed      clustering;

        BIMAPStep( step_time_t time, size_t turbineIx, const ChessboardBiclustering& clustering,
                   const StatsMetrics& metrics, ChessboardBiclusteringsIndexing& indexing )
        : time( time ), turbineIx( turbineIx ), metrics( metrics ), clustering( indexing.index( clustering ) )
        {}

        BIMAPStep( step_time_t time, size_t turbineIx, const ChessboardBiclusteringIndexed& clustering,
                   const StatsMetrics& metrics )
        : time( time ), turbineIx( turbineIx ), metrics( metrics ), clustering( clustering )
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

    ChessboardBiclusteringsIndexing&   _chessboardBiclusteringsIndexing;        /** walker's indexer of biclusterings */
    StepsIndexing               _steps;
    PriorParamsStepsIndexing    _priorParamsSteps;
    LLHWeights                  _probWeights;

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
            ar << boost::serialization::make_nvp( "metrics", stepit->metrics );
            default_serial_t serial = cclus.serial();
            ar << boost::serialization::make_nvp( "chessboardBiclusteringSerial", serial );
            ar << boost::serialization::make_nvp( "clusteringData", cclus );
            LOG_DEBUG2( stepit->iteration << ": " << cclus._objectClusterData.size() );
            ar << boost::serialization::make_nvp( "blocksData", cclus._blocksData );
            ar << boost::serialization::make_nvp( "objectsData", cclus._objectsData );
        }
        ar << boost::serialization::make_nvp( "priorParamSteps", _priorParamsSteps );
        ar << boost::serialization::make_nvp( "probWeights", _probWeights );
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
            StatsMetrics metrics;
            ar >> BOOST_SERIALIZATION_NVP( metrics );
            default_serial_t chessboardBiclusteringSerial;
            ar >> BOOST_SERIALIZATION_NVP( chessboardBiclusteringSerial );
            ChessboardBiclusteringsIndexing::const_serial_iterator pScaffoldIt = _chessboardBiclusteringsIndexing.iterator_to( chessboardBiclusteringSerial );
            if ( pScaffoldIt == _chessboardBiclusteringsIndexing.serial_not_found() ) {
                THROW_RUNTIME_ERROR( "Chessboard biclustering #" << chessboardBiclusteringSerial << " not found the index" );
            }
            ChessboardBiclusteringData clusteringData;
            ar >> BOOST_SERIALIZATION_NVP( clusteringData );
            ChessboardBiclusteringIndexed::block_data_map_type blocksData;
            ar >> BOOST_SERIALIZATION_NVP( blocksData );
            multiple_map_t      objectsData;
            ar >> BOOST_SERIALIZATION_NVP( objectsData );

            _steps.get<0>().push_back( BIMAPStep( time, turbineIx, ChessboardBiclusteringIndexed( *pScaffoldIt, blocksData, objectsData, clusteringData ), metrics ) );
        }
        ar >> boost::serialization::make_nvp( "priorParamSteps", _priorParamsSteps );
        ar >> boost::serialization::make_nvp( "probWeights", _probWeights );
    }

protected:
    friend class BIMAPSampler;

public:
    BIMAPWalk( ChessboardBiclusteringsIndexing& chessboardBiclusteringsIndexing,
               const LLHWeights& probWeights = LLHWeights() );

    typedef StepsIndexing::const_iterator const_step_iterator;
    typedef PriorParamsStepsIndexing::const_iterator const_priors_step_iterator;

    void step( step_time_t time, size_t turbineIx, const ChessboardBiclustering& clu, const StatsMetrics& metrics ) {
        _steps.get<0>().push_back( BIMAPStep( time, turbineIx, clu, metrics, _chessboardBiclusteringsIndexing ) );
    }
    void step( step_time_t time, size_t turbineIx, const ChessboardBiclusteringDerivedPriors& priors ) {
        _priorParamsSteps.get<0>().push_back( BIMAPPriorParamsStep( time, turbineIx, priors ) );
    }

    const ChessboardBiclusteringsIndexing& indexing() const {
        return ( _chessboardBiclusteringsIndexing );
    }
    ChessboardBiclusteringsIndexing& indexing() {
        return ( _chessboardBiclusteringsIndexing );
    }

    const LLHWeights& probWeights() const {
        return ( _probWeights );
    }
    void setProbWeights( const LLHWeights& weights ) {
        _probWeights = weights;
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
