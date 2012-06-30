//#include "cemm/bimap/dynamic_bitset_utils.h"
#include <cemm/math/Distributions.h>

#include "cemm/bimap/CellsLLHEval.h"

namespace cemm { namespace bimap {

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
    LOG_DEBUG2( "Cell[N](" << objectIx << "," << probeIx << ") pdf=" << res );
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
    LOG_DEBUG2( "Cell[S](" << objectIx << "," << probeIx << ") pdf=" << res );
    return ( res );
}

} }
