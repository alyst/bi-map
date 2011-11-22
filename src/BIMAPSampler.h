#pragma once

#include "BasicTypedefs.h"

#include "ChessboardBiclustering.h"
#include "ChessboardBiclusteringCrossover.h"
#include "ChessboardBiclusteringGibbsSampler.h"
#include "BIMAPWalk.h"
#include "math/PitmanYorProcess.h"
#include "eesampler/TurbineCascadeUnit.h"

/**
    Implementation of DynamicParticle concept for ChessboardBiclusteringFit.
 */
class DynamicChessboardBiclustering {
public:
    typedef StaticChessboardBiclustering static_particle_type;
    typedef ChessboardBiclusteringEnergyEval static_particle_energy_eval_type;

private:
    const PrecomputedData&          _precomputed;
    const ChessboardBiclusteringPriors&    _priors;
    const ChessboardBiclusteringHyperPriors&   _hyperpriors;
    const GibbsSamplerParams&       _params;
    ChessboardBiclusteringGibbsSampler     _sampler;
    ChessboardBiclusteringFit              _pClus;

public:
    DynamicChessboardBiclustering( const gsl_rng* rndNumGen,
                            const PrecomputedData& precomputed,
                            const ChessboardBiclusteringPriors& priors,
                            const ChessboardBiclusteringHyperPriors& hyperpriors,
                            const GibbsSamplerParams& params,
                            double minEnergy, double temperature )
    : _precomputed( precomputed )
    , _priors( priors )
    , _hyperpriors( hyperpriors )
    , _params( params )
    , _sampler( rndNumGen, precomputed, hyperpriors, params,
                SamplingTransform( -minEnergy, temperature ) )
    , _pClus( _precomputed, _priors )
    {
    }

    DynamicChessboardBiclustering& operator=( const ChessboardBiclusteringFit& clus )
    {
        _pClus = clus;
        return ( *this );
    }

    DynamicChessboardBiclustering& operator=( const static_particle_type& particle )
    {
        _pClus = particle;
        return ( *this );
    }

    void setParams( double minEnergy, double temperature )
    {
        _sampler.setTransform( SamplingTransform( -minEnergy, temperature ) );
    }

    operator const ChessboardBiclusteringFit&() const
    {
        return ( _pClus );
    }

    operator static_particle_type() const
    {
        return ( static_particle_type( _pClus ) );
    }

    double energy() const
    {
        return ( -_pClus.totalLnP() );
    }

    void iterate()
    {
        _sampler.doSamplingStep( _pClus );
    }

    const static_particle_energy_eval_type& energyEval() const {
        return ( ChessboardBiclusteringEnergyEval() );
    }
    void setEnergyEval( const static_particle_energy_eval_type& energyEval ) {
    }
};

/**
    Implementation of DynamicParticleFactory concept for ChessboardBiclusteringFit.
 */
struct DynamicChessboardBiclusteringFactory {
    typedef DynamicChessboardBiclustering dynamic_particle_type;
    typedef StaticChessboardBiclustering static_particle_type;
    typedef ChessboardBiclusteringEnergyEval static_particle_energy_eval_type;

    const gsl_rng*                      rndNumGen;
    const PrecomputedData&              precomputed;
    const ChessboardBiclusteringPriors&        priors;
    const ChessboardBiclusteringHyperPriors&   hyperpriors;
    const GibbsSamplerParams&           params;

    DynamicChessboardBiclusteringFactory( 
            const gsl_rng*                      rndNumGen,
            const PrecomputedData&              precomputed, 
            const ChessboardBiclusteringPriors&        priors, 
            const ChessboardBiclusteringHyperPriors&   hyperpriors,
            const GibbsSamplerParams&           params
    ) : rndNumGen( rndNumGen )
    , precomputed( precomputed )
    , priors( priors )
    , hyperpriors( hyperpriors )
    , params( params )
    {}

    DynamicChessboardBiclustering* operator()( double minEnergy, double temperature ) const
    {
        return ( new DynamicChessboardBiclustering( rndNumGen, precomputed,
                                             priors, hyperpriors, params, 
                                             minEnergy, temperature ) );
    }
    ChessboardBiclusteringEnergyEval updateEnergyEval( const ChessboardBiclusteringEnergyEval& energyEval,
        std::vector<StatsMetrics>& energyLandscape ) const;
};

/**
    Protein Subcomplex Inference By Bayesian Biclustering (BIMAP) evaluator.
 */
class BIMAPSamplerHelper {
public:
    typedef StaticChessboardBiclustering particle_type;

    const PrecomputedData&              precomputed;
    const ChessboardBiclusteringHyperPriors&   hyperpriors;
    const ChessboardBiclusteringPriors&        priors;
    const GibbsSamplerParams&           params;

    static const gsl_rng_type*          gslRngType;
    gsl_rng*                            rndNumGen;

    DynamicChessboardBiclusteringFactory       dynamicCrossCluFactory;
    ChessboardBiclusteringCrossoverGenerator   crossoverGenerator;

protected:
    GeometricDistribution randomNoiseParams( const ChessboardBiclustering& clus ) const;
    double randomProbeSignal( const ChessboardBiclusteringDerivedPriors& priors ) const;

    void generateMissingProbeSignals( ChessboardBiclustering& clus, const object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const;

    static gsl_rng* CreateRndNumGen( size_t rndNumGenSeed = 0 );

public:
    BIMAPSamplerHelper( const PrecomputedData&             precomputed,
                         const ChessboardBiclusteringHyperPriors&  hyperpriors, 
                         const ChessboardBiclusteringPriors&       priors, 
                         const GibbsSamplerParams&          params,
                         size_t                             rndNumGenSeed = 0,
                         TurbineCascadeExecutionMonitor*    pMonitor = NULL );
    ~BIMAPSamplerHelper();

    const OPAData& data() const {
        return ( precomputed.data() );
    }
    ChessboardBiclusteringDerivedPriors randomPriors() const;
    void randomizeMissingData( ChessboardBiclustering& clus ) const;
    ChessboardBiclustering randomClustering() const;
    ChessboardBiclustering trivialClustering() const;
};

struct BIMAPSampleCollectorParams {
    size_t    walkSamples;          /** number of chessboard biclustering samples to collect */
    size_t    priorsStoragePeriod;  /** period (# of clustering samples) at which prior parameters are sampled
                                        @note i.e. prior samples are taken every priorsStoragePeriod * clusteringStoragePeriod */
    size_t    reportingPeriod;      /** period (# of samples) as which collector logs how many it collected */
    double    maxReportingDelay;    /** max time (in seconds) without logging # of collected samples */

    BIMAPSampleCollectorParams();
};

class BIMAPSampleCollector {
private:
    BIMAPWalk                          _walk;
    const BIMAPSampleCollectorParams&  _params;
    double                             _lastSampleReportTime;

public:
    typedef StaticChessboardBiclustering particle_type;
    typedef ChessboardBiclusteringEnergyEval particle_energy_eval_type;

    BIMAPSampleCollector( ChessboardBiclusteringsIndexing& chessboardBiclusteringsIndexing,
                          const BIMAPSampleCollectorParams& params );

    const BIMAPWalk& walk() const { return ( _walk ); }
    bool storeSample( double time, turbine_ix_t originIx,
                      const particle_type& particle, const particle_energy_eval_type& energyEval );
};

BIMAPWalk BIMAPSampler_run(
    const BIMAPSamplerHelper&  helper,
    ChessboardBiclusteringsIndexing&   chessboardBiclusteringsIndexing,
    const GibbsSamplerParams&   gibbsSamplerParams,
    const TurbineCascadeParams& eeCascadeParams,           /// parameters of equi-energy sample
    const BIMAPSampleCollectorParams&  collectorParams,   /// parameters of sample collector
    const ChessboardBiclustering&      iniClus,                   /** initial clustering */
    const TurbineCascadeExecutionMonitor*   pMon = NULL
);
