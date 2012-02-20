#pragma once

#include "cemm/bimap/BasicTypedefs.h"

#include "cemm/bimap/ChessboardBiclusteringFit.h"
#include "cemm/bimap/ChessboardBiclusteringPriorEval.h"
#include "cemm/bimap/ChessboardBiclusteringLLHEval.h"
#include "cemm/bimap/ChessboardBiclusteringStructureLLHEval.h"

namespace cemm { namespace bimap {

inline ChessboardBiclusteringPriorEval PriorEval( const ChessboardBiclusteringFit& fit )
{
    return ( ChessboardBiclusteringPriorEval( fit.data(), fit.priors(), (const ChessboardBiclustering&)fit ) );
}

inline ChessboardBiclusteringLLHEval LLHEval( const ChessboardBiclusteringFit& fit )
{
    return ( ChessboardBiclusteringLLHEval( fit.signalNoiseCache(), (const ChessboardBiclustering&)fit ) );
}

inline ChessboardBiclusteringStructureLLHEval StructureLLHEval( const ChessboardBiclusteringFit& fit )
{
    return ( ChessboardBiclusteringStructureLLHEval( fit.precomputed(), fit.priors().cellEnablementProb,
                                              fit.priors().probesCluOffObjectsCluRate, fit.priors().objectsCluOffProbesCluRate,
                                              (const ChessboardBiclustering&)fit ) );
}

} }