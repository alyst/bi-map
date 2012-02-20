#include "cemm/bimap/BIMAPWalk.h"

namespace cemm { namespace bimap {

BIMAPWalk::BIMAPWalk(
    ChessboardBiclusteringsIndexing&  chessboardBiclusteringsIndexing,
    const LLHWeights&                 probWeights
) : _chessboardBiclusteringsIndexing( chessboardBiclusteringsIndexing )
  , _probWeights( probWeights )
{
}

bool BIMAPWalk::check() const
{
    for ( const_step_iterator stepIt = stepsBegin(); stepIt != stepsEnd(); ++stepIt ) {
        LOG_DEBUG2( "Checking iteration #" << stepIt->iteration );
        stepIt->clustering.check();
    }
    return ( true );
}

size_t BIMAPWalk::filterSteps(
    size_t minScaffoldCounts,
    size_t minObjsPtnCounts,
    size_t minProbesPtnCounts
){
    size_t  filteredCount = 0;
    for ( StepsIndexing::nth_index_iterator<0>::type stepIt = _steps.get<0>().begin(); stepIt != _steps.get<0>().end(); ) {
        const BIMAPStep& step = *stepIt;
        if ( ( step.clustering.scaffoldRefCount() < minScaffoldCounts + 1 )
          || ( step.clustering.objectsPartitionRefCount() < minObjsPtnCounts + 1 )
          || ( step.clustering.probesPartitionRefCount() < minProbesPtnCounts + 1 )
        ){
            stepIt = _steps.get<0>().erase( stepIt );
            filteredCount++;
        }
        else {
            ++stepIt;
        }
    }
    return ( filteredCount );
}

} }