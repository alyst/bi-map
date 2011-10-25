#include "dynamic_bitset_utils.h"
#include "CellsLLHEval.h"
#include "DataSignalNoiseCache.h"

BOOST_CLASS_EXPORT( DataSignalNoiseCache )

/**
 *  Serialization-only ctor, called by load_construct_data().
 *  @todo make it private
 */
DataSignalNoiseCache::DataSignalNoiseCache(
    const PrecomputedData&     precomputed, 
    const ChessboardBiclusteringPriors&        priors
) : _precomputed( precomputed )
  , _priors( priors )
  , _noiseLnPdfValid( false )
  , _signalLnPdfValid( false )
{
}

DataSignalNoiseCache::DataSignalNoiseCache(
    const PrecomputedData& precomputed, 
    const ChessboardBiclusteringPriors&    priors,
    const ChessboardBiclustering&          clus
) : _precomputed( precomputed )
  , _priors( priors )
  , _clusData( clus.clusteringData() )
  , _noiseLnPdfValid( false )
  , _signalLnPdfValid( false )
{
}

void DataSignalNoiseCache::setChessboardBiclusteringData(
    const ChessboardBiclusteringData& data
){
    if ( _clusData._baselineSignalParams != data._baselineSignalParams
      || _clusData._derivedPriors.signalPrior != data._derivedPriors.signalPrior
    ){
        _signalLnPdfValid = false;
    }
    if ( _clusData._noiseParams != data._noiseParams ) {
        _noiseLnPdfValid = false;
    }
    _clusData = data;
}

void DataSignalNoiseCache::evalMatrices()
{
    CellsLLHEval eval( _precomputed, _clusData );
    if ( _noiseLnPdf.size1() != data().objectsCount() || _noiseLnPdf.size2() != data().probesCount() ) {
        _noiseLnPdfValid = false;
        _noiseLnPdf.reset( data().objectsCount(), data().probesCount(), unset() );
    }
    if ( _signalLnPdf.size1() != data().objectsCount() || _signalLnPdf.size2() != data().probesCount() ) {
        _signalLnPdfValid = false;
        _signalLnPdf.reset( data().objectsCount(), data().probesCount(), unset() );
    }
    for ( probe_index_t probeIx = 0; probeIx < data().probesCount(); probeIx++ ) {
        for ( object_index_t objIx = 0; objIx < data().objectsCount(); objIx++ ) {
            if ( !_noiseLnPdfValid ) {
                _noiseLnPdf( objIx, probeIx ) = eval.cellNoiseLnPdf( objIx, probeIx );
            }
            if ( !_signalLnPdfValid ) {
                _signalLnPdf( objIx, probeIx ) = eval.cellSignalLnPdf( objIx, probeIx );
            }
        }
    }
    _noiseLnPdfValid = true;
    _signalLnPdfValid = true;
}

log_prob_t DataSignalNoiseCache::noiseLLH(
    const object_set_t& objects, 
    const probe_bitset_t& probes
) const {
    double res = 0;
    for ( object_set_t::const_iterator objIt = objects.begin(); objIt != objects.end(); ++objIt ) {
        const object_index_t objIx = *objIt;
        foreach_bit( probe_index_t, probeIx, probes ) {
            res += _noiseLnPdf( objIx, probeIx );
        }
    }
    BOOST_ASSERT( std::isfinite( res ) );
    return ( res );
}
