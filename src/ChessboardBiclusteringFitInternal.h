#pragma once

#include "BasicTypedefs.h"

#include "ChessboardBiclusteringFit.h"
#include "ChessboardBiclusteringPriorEval.h"
#include "ChessboardBiclusteringLLHEval.h"
#include "ChessboardBiclusteringStructureLLHEval.h"

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
