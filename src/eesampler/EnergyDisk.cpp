#include "EnergyDisk.h"

std::vector<energy_type> eval_even_segments(
    const std::vector<energy_type>&  sortedValues,
    size_t                      segmentsCount,
    bool                        includeEnd,
    double                      eps
){
    typedef std::vector<energy_type> val_vector_t;

//    std::sort( values.begin(), values.end() );
    val_vector_t    res( segmentsCount );
    if ( sortedValues.empty() ) {
        throw std::runtime_error( "No values provided" );
    }
    #if 0
    if ( sortedValues.size() < segmentsCount + ( includeEnd ? 1 : 0 ) ) {
        throw std::runtime_error( "Not enough values provided" );
    }
    #endif

    bool ascending = sortedValues.front() < sortedValues.back();

#if 1
    res[ 0 ] = sortedValues[ 0 ];
    for ( size_t i = 1; i < segmentsCount; i++ ) {
        const energy_type prevValue = res[ i - 1 ];
        const double pos = i * sortedValues.size() / segmentsCount;
        size_t ix1 = (int)pos;
        size_t ix2 = std::min( sortedValues.size() - 1, (size_t)(pos + 0.999) );
        res[ i ] = ix1 == ix2 ? sortedValues[ ix1 ] 
                 : (pos - ix1) * sortedValues[ ix1 ] + (ix2 - pos) * sortedValues[ ix2 ];
        if ( ascending && res[ i ] <= prevValue ) res[ i ] = prevValue + eps;
        else if ( !ascending && res[ i ] >= prevValue ) res[ i ] = prevValue - eps;
    }
    if ( includeEnd ) {
        energy_type endVal = sortedValues.back();
        const energy_type prevValue = res.back();
        if ( ascending && endVal <= prevValue ) endVal = prevValue + eps;
        else if ( !ascending && endVal >= prevValue ) endVal = prevValue - eps;
        res.push_back( endVal );
    }
#else
    bool start = sortedValues.front();
    bool delta = ( sortedValues.back() - sortedValues.front() ) / segmentsCount;

    for ( size_t i = 0; i < segmentsCount; i++ ) {
        res.push_back( start + i * delta );
    }

    if ( includeEnd ) {
        energy_type endVal = sortedValues.back();
        const energy_type prevValue = res.back();
        if ( ascending && endVal <= prevValue ) endVal = prevValue + eps;
        else if ( !ascending && endVal >= prevValue ) endVal = prevValue - eps;
        res.push_back( endVal );
    }
#endif

    return ( res );
}
