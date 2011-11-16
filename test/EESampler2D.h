#include <BasicTypedefs.h>

#include <set>

#include <math/Distributions.h>

#include <eesampler/TurbineCascadeUnit.h>
#include <mcmc/MetropolisHastingsStep.h>
#include <mcmc/TransitionDistributions.h>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>

struct Particle2d: public std::pair<double, double> {
    typedef std::pair<double, double> super_type;

    Particle2d( double x = 0, double y = 0 ) : super_type( x, y ) {};
    Particle2d( const super_type& a ) : super_type( a ) {};

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp( "base", boost::serialization::base_object<super_type>( *this ) );
    }

    friend std::ostream& operator<<( std::ostream& out, const Particle2d& p )
    {
        out << "(" << p.first << "," << p.second << ")";
        return ( out );
    }
};

BOOST_CLASS_EXPORT( Particle2d );

struct Particle2dDistrParam {
    typedef Particle2d value_type;

    const value_type  center;
    const value_type  sigma;
    const double ratio;

    Particle2dDistrParam(
        const value_type& center,
        const value_type& sigma,
        double ratio )
        : center( center ), sigma( sigma ), ratio( ratio )
    {}
};

struct StaticParticle2d : public Particle2d {
    double      _energy;

    StaticParticle2d()
    : Particle2d(), _energy( gsl_nan() )
    {}

    StaticParticle2d(
        const Particle2d&   particle,
        double              energy
    ) : Particle2d( particle ), _energy( energy )
    {}

    double energy() const {
        return ( _energy );
    }

    StaticParticle2d& operator=( const Particle2d& a )
    {
        Particle2d::operator=( a );
        _energy = gsl_nan();
        return ( *this );
    }

    friend std::ostream& operator<<( std::ostream& out, const StaticParticle2d& p )
    {
        return ( out << (const Particle2d&)p );
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp( "base", boost::serialization::base_object<Particle2d>( *this ) );
        ar & boost::serialization::make_nvp( "energy", _energy );
    }
};

struct Particle2dEval {
    typedef Particle2d value_type;

    const Particle2dDistrParam param;

    Particle2dEval(
        const Particle2dDistrParam& param
    ) : param( param )
    {}

    double operator()( const value_type& val ) const {
        return ( log( gsl_ran_bivariate_gaussian_pdf( val.first - param.center.first, val.second - param.center.second,
                                                      param.sigma.first, param.sigma.second, param.ratio ) ) /*
              + log( gsl_ran_bivariate_gaussian_pdf( val.first - 3, val.second - 5, 0.1, 0.2, 0.5 ) ) */ );
    }
};

struct Particle2dEnergyEval {
    typedef StaticParticle2d particle_type;

    double operator()( const StaticParticle2d& val ) const {
        return ( val.energy() );
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    }
};

struct DynamicParticle2d {
    typedef Particle2d particle_type;
    typedef Particle2dEnergyEval static_particle_energy_eval_type;

    const Particle2dEval eval;
    const gsl_rng*  rng;

    double          minEnergy;
    double          temperature;
    mutable StaticParticle2d    _pcl;

    DynamicParticle2d( const gsl_rng* rng, double minEnergy, double temperature,
                       const Particle2dDistrParam& param )
    : eval( param )
    , rng( rng )
    , minEnergy( minEnergy )
    , temperature( temperature )
    {}

    DynamicParticle2d& operator=(const StaticParticle2d& pcl)
    {
        _pcl = pcl;
        return ( *this );
    }

    void setParams( double minEnergy, double temperature )
    {
        DynamicParticle2d::temperature = temperature;
        DynamicParticle2d::minEnergy = minEnergy;
    }

    operator const StaticParticle2d&() const
    {
        return ( _pcl );
    }

    void iterate()
    {
        energy();
        _pcl = MetropolisHastringsPosteriorSample<Particle2d>( 
                rng, BivariateGaussianTransitionDistribution( 0.1 * sqrt(temperature), 0.1 * sqrt(temperature), 0 ),
                eval, DegeneratedDistribution<Particle2d>(),
                _pcl, -_pcl._energy, SamplingTransform( -minEnergy, temperature ) ).value;
        _pcl._energy = gsl_nan();
    }

    double energy() const {
        if ( gsl_isnan( _pcl._energy ) ) {
            _pcl._energy = -eval( _pcl );
        }
        return ( -_pcl._energy );
    }

    static_particle_energy_eval_type staticParticleEnergyEval() const {
        return ( static_particle_energy_eval_type() );
    }
};

struct DynamicParticle2dFactory {
    typedef DynamicParticle2d dynamic_particle_type;
    typedef Particle2dEnergyEval static_particle_energy_eval_type;
    typedef Particle2dEnergyEval::particle_type static_particle_type;

    const gsl_rng* rng;
    const Particle2dDistrParam param;

    DynamicParticle2dFactory( const gsl_rng* rng, const Particle2dDistrParam param )
    : rng( rng )
    , param( param )
    {}

    dynamic_particle_type* operator()( double minEnergy, double temperature ) const {
        return ( new dynamic_particle_type( rng, minEnergy, temperature, param ) );
    }
};

struct Particle2dEquiEnergySamplerTestParam
{
    TurbineCascadeParams    eesParams;
    Particle2dDistrParam    distrParams;
    size_t                  samplesCount;

    Particle2dEquiEnergySamplerTestParam(
        const Particle2d& mean,
        const Particle2d& sigma,
        double ratio,
        size_t samplesCount,
        size_t levelsCount,
        size_t turbinesCount,
        double eeJumpRate = 0.3,
        double generateRate = 0.0
    ) : distrParams( mean, sigma, ratio )
      , samplesCount( samplesCount )
    {
        eesParams.levelsCount = levelsCount;
        eesParams.turbinesCount = turbinesCount;
        eesParams.turbineParams.eeJumpRate = eeJumpRate;
        eesParams.turbineParams.generateRate = generateRate;
    }
};

struct Particle2dCollector {
    const size_t samplesToCollect;

    std::vector<StaticParticle2d> samples;

    Particle2dCollector( size_t samplesToCollect = 0 )
    : samplesToCollect( samplesToCollect )
    {}

    bool storeSample( double time, turbine_ix_t originIx, const StaticParticle2d& particle )
    {
        samples.push_back( particle );
        LOG_DEBUG1_IF( ( samples.size() % 100 == 0 ), samples.size() << " of " << samplesToCollect << " samples collected" );
        return ( samples.size() >= samplesToCollect );
    }
};

class Particle2dInterpolationGenerator {
private:
    typedef ParticleCache<StaticParticle2d>::energies_proxy_type particle_cache_type;

public:
    typedef std::vector<StaticParticle2d> particle_container_type;

    particle_container_type operator()( const gsl_rng* rng,
                                        const particle_cache_type& cache ) const
    {
        particle_container_type res;
        if ( cache.size() < 2 )                   return ( res );
        size_t ix1 = gsl_rng_uniform_int( rng, cache.size()-1 );
        size_t ix2 = gsl_rng_uniform_int( rng, cache.size()-2 );
        if ( ix2 == ix1 ) ix2++;
        const Particle2d& pt1 = cache.nthParticle( ix1 )->particle;
        const Particle2d& pt2 = cache.nthParticle( ix2 )->particle;
        const double k = gsl_rng_uniform( rng );

        Particle2d lerp1;
        Particle2d lerp2;
        lerp1.first = pt1.first * k + pt2.first * ( 1 - k );
        lerp1.second = pt1.second * k + pt2.second * ( 1 - k );
        //lerp1.val = BivariateGaussianTransitionDistribution( 0.1, 0.1, 0 ).generate( rng, lerp1.val );
        //lerp1.first = pit1->particle.first;
        //lerp1.second = pit2->particle.second;
        //lerp2.first = pit2->particle.first;
        //lerp2.second = pit1->particle.second;

        res.push_back( StaticParticle2d( lerp1, gsl_nan() ) );
        //res.push_back( lerp2 );
        return ( res );
    }
};
