#pragma once

#include "BasicTypedefs.h"

#include "ChessboardBiclusteringFit.h"
#include "ChessboardBiclusteringGibbsSampler.h"

#include "eesampler/EnergyDisk.h"

class ChessboardBiclusteringCrossoverGenerator {
public:
    typedef StaticChessboardBiclustering particle_type;
    typedef std::vector<particle_type> particle_container_type;

private:
    const PrecomputedData&          _precomputed;
    const ChessboardBiclusteringPriors&    _priors;
    mutable ChessboardBiclusteringGibbsSampler _sampler;
    typedef EnergyDisk<particle_type>::energies_proxy_type energy_disk_type;

    void push_to_result( particle_container_type& res, const ChessboardBiclusteringFit& ptn ) const;

public:
    ChessboardBiclusteringCrossoverGenerator( 
                   const gsl_rng* rndNumGen,
                   const PrecomputedData& precomputed,
                   const ChessboardBiclusteringPriors& priors,
                   const ChessboardBiclusteringHyperPriors& hyperpriors );

    virtual particle_container_type operator()( const gsl_rng* rng, 
                                                const energy_disk_type& disk ) const;
    virtual ~ChessboardBiclusteringCrossoverGenerator()
    {
    }
};
