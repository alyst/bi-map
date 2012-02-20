#pragma once

#include "BasicTypedefs.h"

#include "ChessboardBiclusteringFit.h"
#include "ChessboardBiclusteringGibbsSampler.h"

#include "cemm/eesampler/ParticleCache.h"

namespace cemm { namespace bimap {

using cemm::eesampler::ParticleCache;

class ChessboardBiclusteringCrossoverGenerator {
public:
    typedef StaticChessboardBiclustering particle_type;
    typedef std::vector<particle_type> particle_container_type;

private:
    const PrecomputedData&          _precomputed;
    const ChessboardBiclusteringPriors&    _priors;
    mutable ChessboardBiclusteringGibbsSampler _sampler;
    typedef ParticleCache<particle_type>::energies_proxy_type particle_cache_type;

    void push_to_result( particle_container_type& res, const ChessboardBiclusteringFit& ptn ) const;

public:
    ChessboardBiclusteringCrossoverGenerator( 
                   const gsl_rng* rndNumGen,
                   const PrecomputedData& precomputed,
                   const ChessboardBiclusteringPriors& priors,
                   const ChessboardBiclusteringHyperPriors& hyperpriors,
                   const ChessboardBiclusteringEnergyEval& energyEval );

    virtual particle_container_type operator()( const gsl_rng* rng, 
                                                const particle_cache_type& cache ) const;
    virtual ~ChessboardBiclusteringCrossoverGenerator()
    {
    }
};

} }