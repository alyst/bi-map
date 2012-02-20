#pragma once

#include "BasicTypedefs.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>

#include <cemm/containers/array2d.h>
#include <cemm/statically_tracked.h>

#include "PrecomputedData.h"
#include "ChessboardBiclustering.h"

namespace cemm { namespace bimap {

/**
 *  Caches probabilities of OPAData cells to be from signal/noise.
 */
class DataSignalNoiseCache
{
public:
    typedef array2d<log_prob_t> lnprob_matrix_type;

protected:
    const PrecomputedData&          _precomputed;
    const ChessboardBiclusteringPriors&    _priors;
    ChessboardBiclusteringData             _clusData;

    bool                _noiseLnPdfValid;
    lnprob_matrix_type  _noiseLnPdf;
    bool                _signalLnPdfValid;
    lnprob_matrix_type  _signalLnPdf;

    void evalMatrices();

private:
    friend class boost::serialization::access;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> boost::serialization::make_nvp( "clusData", _clusData );
        ar >> boost::serialization::make_nvp( "noiseLnPdfValid", _noiseLnPdfValid );
        if ( _noiseLnPdfValid ) {
            ar >> boost::serialization::make_nvp( "noiseLnPdf", _noiseLnPdf );
        }
        else {
            _noiseLnPdf.reset( data().objectsCount(), data().probesCount(), unset<log_prob_t>() );
        }
        ar >> boost::serialization::make_nvp( "signalLnPdfValid", _signalLnPdfValid );
        if ( _signalLnPdfValid ) {
            ar >> boost::serialization::make_nvp( "signalLnPdf", _signalLnPdf );
        }
        else {
            _signalLnPdf.reset( data().objectsCount(), data().probesCount(), unset<log_prob_t>() );
        }
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << boost::serialization::make_nvp( "clusData", _clusData );
        ar << boost::serialization::make_nvp( "noiseLnPdfValid", _noiseLnPdfValid );
        if ( _noiseLnPdfValid ) {
            ar << boost::serialization::make_nvp( "noiseLnPdf", _noiseLnPdf );
        }
        ar << boost::serialization::make_nvp( "signalLnPdfValid", _signalLnPdfValid );
        if ( _signalLnPdfValid ) {
            ar << boost::serialization::make_nvp( "signalLnPdf", _signalLnPdf );
        }
    }

public:
    DataSignalNoiseCache( const PrecomputedData& precomputed, const ChessboardBiclusteringPriors& priors );
    DataSignalNoiseCache( const PrecomputedData& precomputed, const ChessboardBiclusteringPriors& priors, const ChessboardBiclustering& clus );

    void setChessboardBiclusteringData( const ChessboardBiclusteringData& data );

    const PrecomputedData& precomputed() const {
        return ( _precomputed );
    }

    const OPAData& data() const {
        return ( _precomputed.data() );
    }

    const CellSignalParams& signalParams() const {
        return ( _precomputed.signalParams() );
    }
    const ChessboardBiclusteringPriors& priors() const {
        return ( _priors );
    }

    void update() {
        if ( !_noiseLnPdfValid || !_signalLnPdfValid ) evalMatrices();
    }

    const lnprob_matrix_type& signalLnPdf() const {
        return ( _signalLnPdf );
    }

    log_prob_t signalLLH( const object_set_t& objects, const probe_bitset_t& probes ) const;
    log_prob_t noiseLLH( const object_set_t& objects, const probe_bitset_t& probes ) const;
};

} }

namespace boost { namespace serialization {

template<class Archive>
inline void save_construct_data(
    Archive & ar, const cemm::bimap::DataSignalNoiseCache* t, const unsigned int file_version
){
    // save data required to construct instance
    const cemm::bimap::PrecomputedData* pPrecomputed = &t->precomputed();
    cemm::statically_tracked<cemm::bimap::PrecomputedData> tPrecomputed( "precomputed", pPrecomputed );
    ar << boost::serialization::make_nvp( "precomputed", tPrecomputed );

    const cemm::bimap::ChessboardBiclusteringPriors* pPriors = &t->priors();
    cemm::statically_tracked<cemm::bimap::ChessboardBiclusteringPriors> tPriors( "priors", pPriors );
    ar << boost::serialization::make_nvp( "priors", tPriors );
}

template<class Archive>
inline void load_construct_data(
    Archive & ar, cemm::bimap::DataSignalNoiseCache* t, const unsigned int file_version
){
    const cemm::bimap::PrecomputedData* pPrecomputed;
    cemm::statically_tracked<cemm::bimap::PrecomputedData> tPrecomputed( "precomputed", pPrecomputed );
    ar >> boost::serialization::make_nvp( "precomputed", tPrecomputed );

    const cemm::bimap::ChessboardBiclusteringPriors* pPriors;
    cemm::statically_tracked<cemm::bimap::ChessboardBiclusteringPriors> tPriors( "priors", pPriors );
    ar >> boost::serialization::make_nvp( "priors", tPriors );
    // invoke inplace constructor to initialize instance of my_class
    ::new(t)cemm::bimap::DataSignalNoiseCache( *pPrecomputed, *pPriors );
}

} }

BOOST_CLASS_IMPLEMENTATION( cemm::bimap::DataSignalNoiseCache, object_serializable )
