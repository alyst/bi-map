#pragma once

#include "../BasicTypedefs.h"

#include <vector>
#include <set>
#include <boost/numeric/ublas/storage.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>

/**
    Sample generated by Pitman-Yor process.
 */
class PitmanYorSample {
public:
    typedef std::size_t size_type;
    typedef std::size_t cluster_index_type;
    typedef std::size_t sample_index_type;
    typedef std::set<sample_index_type> sample_set_type;

private:
    typedef std::vector<cluster_index_type> sample_to_cluster_map;
    typedef std::vector<sample_set_type> cluster_to_samples_map;

    sample_to_cluster_map   _sampleToCluster;
    cluster_to_samples_map  _clusterToSamples;

protected:
    friend class PitmanYorProcess;

    PitmanYorSample()
    {
    }

    void push_next_sample( cluster_index_type clusterIx );

public:
    size_type samplesCount() const {
        return ( _sampleToCluster.size() );
    }

    size_type clustersCount() const {
        return ( _clusterToSamples.size() );
    }

    const sample_set_type& operator[]( cluster_index_type clusterIx ) const {
        return ( cluster( clusterIx ) );
    }

    const sample_set_type& cluster( cluster_index_type clusterIx ) const {
        return ( _clusterToSamples[ clusterIx ] );
    }

    cluster_index_type clusterOfSample( sample_index_type sampleIx ) const {
        return ( _sampleToCluster[ sampleIx ] );
    }
};

/**
    Pitman-Yor process.
    @see Combinatorial Stochastic Processes, Jim Pitman, Two-parameter model.
 */
struct PitmanYorProcess {
    prob_t  concentration;
    prob_t  discount;

    typedef boost::numeric::ublas::unbounded_array<log_prob_t> log_prob_cache_type;

    mutable log_prob_cache_type  __lnPoch_concentration_cache;
    mutable log_prob_cache_type  __lnPoch_discount_cache;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( concentration );
        ar & BOOST_SERIALIZATION_NVP( discount );
    }

    static log_prob_t lnPoch( log_prob_cache_type& lnPoch_cache, double a, size_t n )
    {
        if ( n >= lnPoch_cache.size() ) {
            lnPoch_cache.resize( n + 1, unset() );
        }
        log_prob_t& val = lnPoch_cache[ n ];
        if ( is_unset( val ) ) {
            val = gsl_sf_lnpoch( a, n );
        }
        return ( val );
    }

    typedef std::size_t size_type;
    typedef std::vector<size_type> size_vector_type;
    typedef std::vector<prob_t> weight_vector_type;

    weight_vector_type sampleLogStickBreaks( const gsl_rng* r, const size_vector_type& groupSizes, int curClusterIx = -1, int newClusterIx = -1 ) const;
    weight_vector_type sampleLogWeights( const gsl_rng* r, const size_vector_type& groupSizes, int curClusterIx = -1, int newClusterIx = -1 ) const;

    prob_t clusterAssignmentPrior( size_type clusterSize, size_type clustersCount, size_type elementsCount ) const;

    prob_t expectedClustersCount( size_type samplesCount ) const;

    PitmanYorProcess( prob_t concentration, prob_t discount );
    PitmanYorSample random( const gsl_rng* r, size_type samplesCount ) const;

    /** Base accumulator for prior probability of partitioning calculation */
    struct PriorLnPAccumBase {
        const PitmanYorProcess& py;
        size_t  groupsCount;
        size_t  samplesCount;
        log_prob_t  accum;

        PriorLnPAccumBase( const PitmanYorProcess& py )
        : py( py )
        , groupsCount( 0 ), samplesCount( 0 )
        , accum( 0 )
        {}

        PriorLnPAccumBase& operator<<( size_t groupSize )
        {
            samplesCount += groupSize;
            groupsCount++;
            return ( *this );
        }
    };

    /** see 3.6 */
    struct PYLnPAccum: public PriorLnPAccumBase {
        PYLnPAccum(const PitmanYorProcess& py);

        PYLnPAccum& operator<<( size_t groupSize );
        log_prob_t operator()() const;
    };

    /** see 2.19 */
    struct DirichletLnPAccum: public PriorLnPAccumBase {
        DirichletLnPAccum(const PitmanYorProcess& py);

        DirichletLnPAccum& operator<<( size_t groupSize );
        log_prob_t operator()() const;
    };

    template<class LnPAccum, class GroupSizes>
    log_prob_t lnP(
        const GroupSizes&     groupSizes
    ) const {
        LnPAccum accum( *this );
        for ( size_t i = 0; i < groupSizes.size(); i++ ) {
            size_t groupSize = groupSizes[i];
            if ( groupSize > 0 ) accum << groupSize;
        }
        log_prob_t res = accum();
        BOOST_ASSERT( res <= 0 );
        return ( res );
    }

    /**
        Logarithm of given partition probability.

        From 'Combinatorical Stohastic processes'.
    */
    template<class GroupSizes>
    log_prob_t lnP(
        const GroupSizes&     groupSizes
    ) const {
        return ( discount > 0 
                 ? lnP<PYLnPAccum>( groupSizes ) 
                 : lnP<DirichletLnPAccum>( groupSizes ) );
    }
};
