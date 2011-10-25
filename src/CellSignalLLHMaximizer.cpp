#include "CellSignalLLHMaximizer.h"

#include "dynamic_bitset_utils.h"
#include "math/gsl_exception.h"

/**
 *  Serialization-only version of ctor.
 */
CellSignalLLHMaximizer::CellSignalLLHMaximizer(
    const OPAData&                  data,
    const CellSignalParams&         signalParams,
    ObjectsClusterSignal::distrib_cache_type& scDistribCache,
    bool  evalPreliminarySignals
) : _data( data )
  , _signalParams( signalParams )
  , _scDistribCache( &scDistribCache )
  , _multCache( data.objectsCount(), 1 )
  , _llhMaximizer( gsl_min_fminimizer_alloc( gsl_min_fminimizer_brent ) )
{
    if ( evalPreliminarySignals ) {
        _preliminarySignals = PreliminarySignals( data, signalParams );
        _preliminaryAssaySignals = PreliminaryAssaySignals( data, signalParams );
    }
}

CellSignalLLHMaximizer::~CellSignalLLHMaximizer()
{
    gsl_min_fminimizer_free( _llhMaximizer );
}

/**
 *  Evaluate signal for each object x probe pair,
 *  assuming no clustering and Poisson model with Gamma prior.
 */
CellSignalLLHMaximizer::signals_matrix_type CellSignalLLHMaximizer::PreliminarySignals(
    const OPAData&              data,
    const CellSignalParams&     signalParams
){
    //double signalShapeDelta = log( 1.0 - signalParams.scShape );
    signals_matrix_type res( data.objectsCount(), data.probesCount(), unset() );
    for ( object_index_t objIx = 0; objIx < data.objectsCount(); objIx++ ) {
        const OPAObject& obj = data.object( objIx );
        for ( probe_index_t probeIx = 0; probeIx < data.probesCount(); probeIx++ ) {
            const OPAProbe& probe = data.probe( probeIx );
            // collect all measurements
            size_t scSum = 0;
            for ( size_t assayPos = 0; assayPos < probe.assayIndexes().size(); assayPos++ ) {
                assay_index_t assayIx = probe.assayIndexes()[ assayPos ];
                scSum += ceil( data.measurement( objIx, assayIx ).sc / data.assay( assayIx ).multiplier() );
            }
            // generate posterior preliminary rate of Poisson distribution
            double rate = signalParams.preliminarySignalPrior.poissonPosterior( scSum, probe.assayIndexes().size() ).mean();
            // correct signal with respect to object's sequence length (and other factors)
            res( objIx, probeIx ) = log( rate )
                                       - signalParams.sequenceLengthFactor * log( obj.sequenceLength() )
                                       /*+ signalShapeDelta*/;
        }
    }
#if 0
    for ( object_index_t objIx = 0; objIx < _preliminarySignals.size1(); ++objIx ) {
        for ( probe_index_t probeIx = 0; probeIx < _preliminarySignals.size2(); ++probeIx ) {
            std::cout << boost::format("%.3f ") % _preliminarySignals( objIx, probeIx );
        }
        std::cout << "\n";
    }
#endif
    return ( res );
}

/**
 *  Evaluate signal for each object x probe pair,
 *  assuming no clustering and Poisson model with Gamma prior.
 */
CellSignalLLHMaximizer::signals_matrix_type CellSignalLLHMaximizer::PreliminaryAssaySignals(
    const OPAData&              data,
    const CellSignalParams&     signalParams
){
    //double signalShapeDelta = log( 1.0 - signalParams.scShape );
    signals_matrix_type probeSignals( data.objectsCount(), data.probesCount(), unset() );
    signals_matrix_type res( data.objectsCount(), data.assaysCount(), unset() );
    for ( object_index_t objIx = 0; objIx < data.objectsCount(); objIx++ ) {
        const OPAObject& obj = data.object( objIx );
        const OPAData::celldata_t* pM = &data.measurement( objIx, 0 );
        for ( size_t assayIx = 0; assayIx < data.assaysCount(); assayIx++ ) {
            // collect all measurements
            size_t sc = ceil( pM[ assayIx ].sc / data.assay( assayIx ).multiplier() );
            // generate posterior preliminary rate of Poisson distribution
            double rate = signalParams.preliminarySignalPrior.poissonPosterior( sc, 1 ).mean();
            // correct signal with respect to object's sequence length (and other factors)
            res( objIx, assayIx ) = log( rate )
                                  - signalParams.sequenceLengthFactor * log( obj.sequenceLength() )
                                /*+ signalShapeDelta*/;
        }
    }
    return ( res );
}

/**
 *  Calculate preliminary signal by taking into account all objects and probes.
 */
signal_t CellSignalLLHMaximizer::clusterSignalA(
    const object_set_t&     objects,
    const assay_bitset_t&   assays,
    const multiple_map_t&   multiples
) const {
    // handle simple case
    if ( objects.size() == 0 || assays.empty() ) return ( unset() );
    else if ( objects.size() == 1 && assays.count() == 1 ) {
        return ( _preliminaryAssaySignals( *objects.begin(), assays.find_first() ) );
    }
    // handle multi-object/multi-assay case
    size_t  minLength = std::numeric_limits<size_t>::max();
    for ( object_set_t::const_iterator oit = objects.begin(); oit != objects.end(); ++oit ) {
        const OPAObject& obj = _data.object( *oit );
        if ( minLength > obj.sequenceLength() ) minLength = obj.sequenceLength();
    }
    size_t scSum = 0;
    size_t dataPoints = 0;
    for ( object_set_t::const_iterator oit = objects.begin(); oit != objects.end(); ++oit ) {
        const OPAObject& obj = _data.object( *oit );
        size_t  multiple = multiples[ *oit ];
        const OPAData::celldata_t* pM = &_data.measurement( obj.index(), 0 );
        foreach_bit( assay_index_t, assayIx, assays ) {
            size_t counts = pM[ assayIx ].sc;
            scSum += ceil( counts / pow( (double)obj.sequenceLength() / minLength, signalParams().sequenceLengthFactor ) 
                            / multiple / _data.assay( assayIx ).multiplier() );
            dataPoints++;
        }
    }
    double rate = _signalParams.preliminarySignalPrior.poissonPosterior( scSum, dataPoints ).mean();
    if ( is_unset( rate ) ) THROW_RUNTIME_ERROR( "Bad preliminary signal generated" );
    // correct signal with respect to object's sequence length (and other factors)
    return ( log( rate )
             - _signalParams.sequenceLengthFactor * log( minLength )
             /*+ log( 1 - _signalParams.scShape )*/ );
}


struct CellBlockLLHEvalData {
    const CellSignalLLHMaximizer& data;
    const assay_bitset_t&                   _assays;

    std::vector<ObjectsClusterSignal>       _objBaseSignals;
    std::vector<const OPAData::celldata_t*> _objMeasurements;

    CellBlockLLHEvalData(
        const CellSignalLLHMaximizer&  data,
        const object_set_t&     objects,
        const assay_bitset_t&   assays,
        const CellSignalLLHMaximizer::multiple_map_t&   multiples
    )
    : data( data )
    , _assays( assays )
    , _objBaseSignals( objects.size() )
    , _objMeasurements( objects.size() )
    {
        size_t i = 0;
        for ( object_set_t::const_iterator oit = objects.begin(); oit != objects.end(); ++oit ) {
            _objMeasurements[ i ] = &data.data().measurement( *oit, 0 );
            _objBaseSignals[ i ] = ObjectsClusterSignal( data.signalParams(),
                                                         data.data().object( *oit ),
                                                         multiples[ *oit ] );
            i++;
        }
    }

    double eval( signal_t signal ) const
    {
        double   signalLLH = 0;
        for ( size_t i = 0; i < _objBaseSignals.size(); i++ ) {
            foreach_bit( assay_index_t, assayIx, _assays ) {
                signalLLH += ObjectsClusterSignal( _objBaseSignals[i], data.data().assay( assayIx ), signal )
                             .distribTable( data.scDistribCache() )
                             .lnNormPdf( _objMeasurements[i][assayIx].sc );
            }
        }
        LOG_DEBUG3( "CellBlockLLH(" << signal << ")=" << signalLLH );
        return ( signalLLH );
    }
};

double _CellBlockLLHEval(
    signal_t     signal,
    void*        p
){
    return ( -((const CellBlockLLHEvalData*)p)->eval( signal ) );
}

/**
 *  Finds the signal that maximizes LLH for given cells block.
 *  @FIXME function is not reentrant as _llhMaximizer is reused.
 * 
 *  @return maximal LLH, -infinity if maximization failed
 */
log_prob_t CellSignalLLHMaximizer::maximizeSignalLLHA(
    const object_set_t&     objects,
    const assay_bitset_t&   assays,
    const multiple_map_t&   multiples,
    signal_t&               signal      /** initial approximation of the signal,
                                            on exit contains the solution */
) const {
    double   signalLLH = 0;
    const size_t MaxIter = 20;
    const double SignalError = 5;
    const double SignalAccuracy = 1E-2;
    const double MaxSignal = 5;

    CellBlockLLHEvalData evalData( *this, objects, assays, multiples );
    gsl_function    cellBlockLLHEval;
    cellBlockLLHEval.function = &_CellBlockLLHEval;
    cellBlockLLHEval.params = &evalData;

    gsl_error_handler_t* prev_handler = gsl_set_error_handler( &raise_gsl_exception );

    try {
        gsl_min_fminimizer_set( _llhMaximizer, &cellBlockLLHEval, std::min( signal, MaxSignal ),
                                std::min( signal, MaxSignal ) - SignalError, std::min( signal, MaxSignal ) + SignalError );

        int status;
        size_t iter = 0;
        do {
            iter++;
            status = gsl_min_fminimizer_iterate( _llhMaximizer );

            signal_t signalMin = gsl_min_fminimizer_x_lower( _llhMaximizer );
            signal_t signalMax = gsl_min_fminimizer_x_upper( _llhMaximizer );

            signalLLH = -gsl_min_fminimizer_f_minimum( _llhMaximizer );
            status = gsl_min_test_interval( signalMin, signalMax, SignalAccuracy, 0.0 );
            if ( status == GSL_SUCCESS ) break;
        } while ( status == GSL_CONTINUE && iter < MaxIter );
        signal = gsl_min_fminimizer_x_minimum( _llhMaximizer );
    }
    catch ( const gsl_exception& e ) {
        signalLLH = -std::numeric_limits<log_prob_t>::infinity();
    }
    gsl_set_error_handler( prev_handler );
    return ( signalLLH );
}

/**
 *   Distances between vectors of measurements
 *   based on LLH of cluster containing these 2 elements.
 */
log_prob_t CellSignalLLHMaximizer::evalObjectsCoSignalLLH(
    object_index_t  obj1ix,
    object_index_t  obj2ix
) const {
    typedef GeometricDistribution noise_params_type;
    static noise_params_type noiseModel = noise_params_type::ByFailureRate( 1E-4, 0 );

    assay_bitset_t assays( data().assaysCount() );
    object_set_t objs;
    objs.insert( obj1ix );
    objs.insert( obj2ix );
    const size_t MaxMultiple = 3; // maximum object multiple to test
    log_prob_t maxLLH = unset();
    for ( size_t mult1 = 1; mult1 <= MaxMultiple; mult1++ ) {
        for ( size_t mult2 = 1; mult2 <= MaxMultiple; mult2++ ) {
            if ( mult1 == mult2 && mult1 > 1 ) continue;
            _multCache[ obj1ix ] = mult1;
            _multCache[ obj2ix ] = mult2;

            const OPAData::celldata_t*  cell1DataVec = &data().measurement( obj1ix, 0 );
            const OPAData::celldata_t*  cell2DataVec = &data().measurement( obj2ix, 0 );
            double llh = 0;
            // calculate signal for each assay individually
            for ( assay_index_t assayIx = 0; assayIx < data().assaysCount(); assayIx++ ) {
                assays.set( assayIx );

                double noiseLLH = noiseModel.lnPdf( (cell1DataVec++)->sc )
                                + noiseModel.lnPdf( (cell2DataVec++)->sc ); // q-value for geom.distr

                // evaluate signalLLH using signal that maximizes it
                signal_t signal = clusterSignalA( objs, assays, _multCache );
                log_prob_t  signalLLH = maximizeSignalLLHA( objs, assays,
                                                           _multCache, signal );

                llh += std::max( signalLLH, noiseLLH );
                assays.set( assayIx, false );
            }
            if ( is_unset( maxLLH ) || maxLLH < llh ) maxLLH = llh;
        }
    }
    BOOST_ASSERT( !is_unset( maxLLH ) );
    return ( maxLLH );
}

/**
 *   Distances between vectors of measurements
 *   based on LLH of cluster containing these 2 elements.
 */
symmetric_array2d<log_prob_t> CellSignalLLHMaximizer::evalObjectCoSignalDistances() const
{
    // evaluate matrix
    LOG_INFO( "Calculating object distances..." );
    symmetric_array2d<log_prob_t>   distances( data().objectsCount() );
    for ( object_index_t elm1Ix = 0; elm1Ix < distances.size(); ++elm1Ix ) {
        for ( object_index_t elm2Ix = elm1Ix; elm2Ix < distances.size(); ++elm2Ix ) {
            distances( elm1Ix, elm2Ix ) = -evalObjectsCoSignalLLH( elm1Ix, elm2Ix );
        }
    }
    return ( distances );
}

log_prob_t CellSignalLLHMaximizer::evalProbesCoSignalLLH(
    probe_index_t  probe1ix,
    probe_index_t  probe2ix
) const {
    typedef GeometricDistribution noise_params_type;
    static noise_params_type noiseModel = noise_params_type::ByFailureRate( 1E-4, 0 );

    assay_bitset_t assays( _data.assaysCount() );

    double maxLLH = unset();
    // _multCache.assign( _multCache.size(), 1 ); not required -- we don't care about signal
    const OPAProbe& probe1 = _data.probe( probe1ix );
    const OPAProbe& probe2 = _data.probe( probe2ix );
    for ( assay_container_t::const_iterator a1it = probe1.assayIndexes().begin(); a1it != probe1.assayIndexes().end(); ++a1it ) {
        assays.set( *a1it );
        for ( assay_container_t::const_iterator a2it = probe2.assayIndexes().begin(); a2it != probe2.assayIndexes().end(); ++a2it ) {
            assays.set( *a2it );
            double llh = 0;
            for ( object_index_t objIx = 0; objIx < data().objectsCount(); objIx++ ) {
                object_set_t objs;
                objs.insert( objIx );
                const OPAData::celldata_t* cellDataVec = &data().measurement( objIx, 0 );

                double noiseLLH = 0;
                foreach_bit( assay_index_t, assayIx, assays ) {
                    noiseLLH += noiseModel.lnPdf( cellDataVec[assayIx].sc );
                }

                // evaluate signal LLH using signal that maximizes it
                signal_t signal = clusterSignalA( objs, assays, _multCache );
                log_prob_t  signalLLH = maximizeSignalLLHA( objs, assays,
                                                            _multCache, signal );
                llh += std::max( signalLLH, noiseLLH );
            }
            if ( *a1it != *a2it ) assays.set( *a2it, false );
            maxLLH = is_unset( maxLLH ) ? llh : std::max( maxLLH, llh );
        }
        assays.set( *a1it, false );
    }
    BOOST_ASSERT( !is_unset( maxLLH ) );
    return ( maxLLH );
}

/**
 *   Distances between vectors of measurements
 *   based on LLH of cluster containing these 2 elements.
 */
symmetric_array2d<log_prob_t> CellSignalLLHMaximizer::evalProbeCoSignalDistances() const
{
    // evaluate matrix
    LOG_INFO( "Calculating probes distances..." );
    symmetric_array2d<log_prob_t>   distances( data().probesCount() );
    for ( probe_index_t elm1Ix = 0; elm1Ix < distances.size(); ++elm1Ix ) {
        for ( probe_index_t elm2Ix = elm1Ix; elm2Ix < distances.size(); ++elm2Ix ) {
            distances( elm1Ix, elm2Ix ) = -evalProbesCoSignalLLH( elm1Ix, elm2Ix );
        }
    }
    return ( distances );
}
