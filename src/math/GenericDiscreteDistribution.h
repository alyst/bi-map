#pragma once

#include "../BasicTypedefs.h"

#include <gsl/gsl_randist.h>

#include <vector>

typedef std::vector<prob_t> probability_vector_t;

template<class UnimodalDistribution>
std::pair<probability_vector_t, int> ProbabilitiesTable(
    const UnimodalDistribution& distr, prob_t min_factor = 1E-3 );

/**
    C++ wrapper for \a gsl_ran_discrete_t API.
 */
struct GenericDiscreteDistribution {
private:
    gsl_ran_discrete_t*     _data;
    int                     _offset;

public:
    GenericDiscreteDistribution( const GenericDiscreteDistribution& that ); /// non-copyable, no API to copy gsl_ran_discrete_t

    GenericDiscreteDistribution( const probability_vector_t& pdf, int offset = 0 )
        : _data( gsl_ran_discrete_preproc( pdf.size(), pdf.data() ) )
        , _offset( offset )
    {
    }

    template<class UnimodalDistribution>
    GenericDiscreteDistribution( const UnimodalDistribution& distr, prob_t min_factor = 1E-3 )
    : _data( NULL ), _offset( 0 )
    {
        std::pair<probability_vector_t, int> table = ProbabilitiesTable( distr, min_factor );
        _data = gsl_ran_discrete_preproc( table.first.size(), table.first.data() );
        _offset = table.second;
    }

    ~GenericDiscreteDistribution()
    {
        if ( _data ) {
            gsl_ran_discrete_free( _data );
        }
    }

    GenericDiscreteDistribution& operator=( const GenericDiscreteDistribution& that ); /// non-copyable, no API to copy gsl_ran_discrete_t

    GenericDiscreteDistribution& operator=( const probability_vector_t& pdf )
    {
        if ( !isEmpty() ) {
            gsl_ran_discrete_free( _data );
        }
        _data = gsl_ran_discrete_preproc( pdf.size(), pdf.data() );
        _offset = 0;
        return ( *this );
    }

    bool isEmpty() const {
        return ( _data == NULL );
    }

    size_t random( const gsl_rng* rnd_gen ) const
    {
        if ( isEmpty() )   throw std::runtime_error( "Discrete probability not initialized" );
        return ( gsl_ran_discrete( rnd_gen, _data ) + _offset );
    }
};


template<class UnimodalDistribution>
std::pair<probability_vector_t, int> ProbabilitiesTable(
    const UnimodalDistribution& distr,
    prob_t                      min_factor = 1E-3
){
    probability_vector_t pdf;
    const int mode_ = distr.mode();
    prob_t    p_max = 0;

    for ( int i = mode_; i == mode_ || pdf.back() > min_factor * p_max; i++ ) {
        prob_t p = distr.pdf( i );
        if ( p > p_max ) p_max = p;
        pdf.push_back( p );
        //std::cout << "diff[" << i << "]=" << pdf.back() << std::endl;
    }
    probability_vector_t revPdf;
    for ( int i = mode_ - 1; i >= 0 && ( revPdf.empty() || revPdf.back() > min_factor * p_max ); i-- ) {
        prob_t p = distr.pdf( i );
        if ( p > p_max ) p_max = p;
        revPdf.push_back( p );
        //std::cout << "diff[" << i << "]=" << revPdf.back() << std::endl;
    }
    pdf.insert( pdf.begin(), revPdf.rbegin(), revPdf.rend() );
    #if 0
    for ( int i = 0; i < pdf.size(); i++ ) {
        std::cout << "diff[" << tableOffset + i << "]=" << pdf[ i ] << std::endl;
    }
    #endif
    return ( std::make_pair( pdf, mode_ - revPdf.size() ) );
}
