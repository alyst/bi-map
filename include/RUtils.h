#pragma once

#include <Rcpp/Vector.h>
#include <Rcpp/DataFrame.h>
#include <boost/lexical_cast.hpp>

typedef Rcpp::Vector<VECSXP> RcppParamVector;

template<class ValueType>
bool MaybeReadRParam( ValueType& dest, RcppParamVector& params, const char* paramName )
{
    bool res = false;
    try {
        Rcpp::RObject destExp = params[ std::string( paramName ) ];
        if ( !destExp.isNULL() ) {
            dest = Rcpp::as<ValueType>( destExp );
            res = true;
        }
    }
    catch ( Rcpp::index_out_of_bounds& e ) {
    }
    catch ( Rcpp::not_compatible& e ) {
    }
    std::ostringstream str;
    str << paramName << " = " << dest;
    if ( !res ) str << " (using default)";
    str << "\n";
    Rprintf( str.str().c_str() );
    return ( res );
}

inline void Rprintf_columns( const Rcpp::DataFrame& dataFrame, const char* dataFrameName )
{
    Rcpp::StringVector colnames = dataFrame.names();
    Rprintf( "Data frame '%s' %d columns:\n", dataFrameName, colnames.size() );
    for ( int i = 0; i < colnames.size(); i++ ) {
        Rprintf( "    %d: %s\n", i, ((std::string)( colnames[i] )).c_str() );
    }
}

/**
    R DataFrame row proxy.
 */
class RDataFrameRowProxy {
private:
    Rcpp::DataFrame&    _frame;
    size_t              _rowIx;

public:
    RDataFrameRowProxy( Rcpp::DataFrame& frame, size_t rowIx )
    : _frame( frame ), _rowIx( rowIx )
    {}

    Rcpp::GenericVector::Proxy operator[]( size_t colIx )
    {
        return ( Rcpp::as<Rcpp::GenericVector>( _frame( colIx ) )( _rowIx ) );
    }
};
