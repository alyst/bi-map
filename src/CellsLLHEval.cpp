#include "dynamic_bitset_utils.h"
#include "math/Distributions.h"

#include "CellsLLHEval.h"

CellsLLHEval::CellsLLHEval(
    const PrecomputedData&      precomputed,
    const ChessboardBiclusteringData&  clusData
) : _precomputed( precomputed )
  , _clusData( clusData )
  , _noiseModel( clusData._noiseParams )
{
}

log_prob_t CellsLLHEval::cellNoiseLnPdf(
    object_index_t      objectIx,
    probe_index_t       probeIx
) const {
    const OPAProbe& probe = data().probe( probeIx );
    const OPAAssay* pLastAssay = &data().assay( probe.assayIndexes().back() );
    const assay_index_t firstAssayIx = probe.assayIndexes().front();
    const OPAData::celldata_t*  cellDataVec = &data().measurement( objectIx, firstAssayIx );

    double res = 0;
    for ( const OPAAssay* pAssay = &data().assay( firstAssayIx ); pAssay <= pLastAssay; ++pAssay ) {
        const OPAData::celldata_t&  cellData = *(cellDataVec++);
        res += _noiseModel.lnPdf( cellData.sc );
    }
    LOG_DEBUG2( "Cell[N](" << objIx << "," << probeIx << ") pdf=" << res );
    return ( res );
}

log_prob_t CellsLLHEval::cellSignalLnPdf(
    object_index_t      objectIx,
    probe_index_t       probeIx
) const {
    const OPAProbe& probe = data().probe( probeIx );
    const OPAAssay* pLastAssay = &data().assay( probe.assayIndexes().back() );
    const assay_index_t firstAssayIx = probe.assayIndexes().front();
    const OPAData::celldata_t*  cellDataVec = &data().measurement( objectIx, firstAssayIx );

    double res = 0;

    ObjectsClusterSignal baseSignalModel( _precomputed.signalParams(),
                                          data().object( objectIx ), 1 );
    signal_t    signal = _precomputed.preliminarySignals()( objectIx, probeIx );
    for ( const OPAAssay* pAssay = &data().assay( firstAssayIx ); pAssay <= pLastAssay; ++pAssay ) {
        const OPAData::celldata_t&  cellData = *(cellDataVec++);
        res += ObjectsClusterSignal( baseSignalModel, *pAssay, signal )
               .distribTable( _precomputed.scDistribCache() ).lnPdf( cellData.sc );
    }
    LOG_DEBUG2( "Cell[S](" << objIx << "," << probeIx << ") pdf=" << res );
    return ( res );
}

/**
 *  Likelihood of block being in on or off state.
 */
log_prob_t CellsLLHEval::BlockLLH(
    const object_set_t&         objects,    /** objects of block */
    const probe_bitset_t&       probes,     /** probes of block */
    const lnprob_matrix_type&   lnPdf       /** matrix of log probabilities
                                                for data cells (either noise or min signal) */
){
    double lnPdfSum = 0.0;

    for ( object_set_t::const_iterator objIt = objects.begin(); objIt != objects.end(); ++objIt ) {
        const object_index_t objIx = *objIt;
        foreach_bit( probe_index_t, probeIx, probes ) {
            lnPdfSum += lnPdf( objIx, probeIx );
        }
    }
    return ( lnPdfSum );
}
