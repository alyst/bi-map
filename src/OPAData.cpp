#include <sstream>

#include "dynamic_bitset_utils.h"
#include "math/logmath.h"

#include <fstream>
#include <boost/filesystem/path.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include "statically_tracked.h"
#include "misc_utils.h"

#include "OPAData.h"

ENABLE_STATIC_TRACKING( OPAData )

std::string entity_error::compose_error_msg(
    const std::string& objectClass,
    const std::string& objectName,
    const std::string& message
){
    std::ostringstream  out;
    out << objectClass << " <" << objectName << "> " << message;
    return ( out.str() );
}

OPAData::OPAData( bool mapBaitsToObjects,
                  size_t objectsUniverseSize )
    : _matrixModified( true )
    , _mapBaitsToObjects( mapBaitsToObjects )
    , _objectsUniverseSize( objectsUniverseSize )
{
}

/**
 *  Gets baits of given set of probes.
 */
object_set_t OPAData::probesToBaits(
    const probe_bitset_t& probes
) const {
    object_set_t res;
    foreach_bit( probe_index_t, probeIx, probes ) {
        object_index_t baitIndex = _probes[ probeIx ].baitIndex();
        if ( baitIndex != OBJECT_NA ) res.insert( baitIndex );
    }
    return ( res );
}

/**
 *  Gets set of assays for given sets of preys.
 */
assay_bitset_t OPAData::probesToAssays(
    const probe_bitset_t& probes
) const {
    assay_bitset_t res( assaysCount() );
    foreach_bit( probe_index_t, probeIx, probes ) {
        const OPAProbe& st = probe( probeIx );
        for ( assay_container_t::const_iterator ait = st.assayIndexes().begin(); ait != st.assayIndexes().end(); ++ait ) {
            res.set( *ait );
        }
    }
    return ( res );
}

size_t OPAData::assaysCount( size_t probeIx ) const
{
  return ( probe( probeIx ).assayIndexes().size() );
}

size_t OPAData::assaysCount(const probe_bitset_t& probes ) const
{
    size_t res = 0;
    foreach_bit( probe_index_t, probeIx, probes ) {
        res += assaysCount( probeIx );
    }
    return ( res );
}

OPAData::const_object_ptr_t OPAData::addObject(
    const object_label_t&   label, 
    size_t                  seqLength
){
    if ( _matrix.size1() > 0 ) {
        throw std::runtime_error( "Cannot add object, measurements matrix already initialized" );
    }
    const_object_ptr_t  pPrevObject = object( label );
    if ( pPrevObject ) {
        throw entity_already_exists( "Object", label );
    }
    object_ptr_t pObject( new OPAObject( _objects.size(), label, seqLength ) );
    _objects.push_back( pObject );
    _objectLabelMap[ label ] = pObject;

    return ( pObject );
}

OPAData::const_probe_ptr_t OPAData::addProbe(
    const probe_label_t&    label, 
    const object_label_t&   baitLabel
){
    if ( _matrix.size2() > 0 ) {
        throw std::runtime_error( "Cannot add probe, measurements matrix already initialized" );
    }
    const_object_ptr_t  pBait = NULL;
    if ( isMapBaitsToObjects() ) {
        if ( !baitLabel.empty() ) {
            pBait = object( baitLabel );
            if ( !pBait ) {
                throw entity_not_found( "Bait", baitLabel );
            }
        } else {
            LOG_WARN( "No bait label for probe '" << label << "', assuming no bait" );
        }
    }
    const_probe_ptr_t  pPrevProbe = probe( label );
    if ( pPrevProbe ) {
        throw entity_already_exists( "Probe", label );
    }
    probe_ptr_t pProbe( new OPAProbe( _probes.size(), label, pBait ? pBait->_index : OBJECT_NA ) );
    _probes.push_back( pProbe );
    _probeLabelMap[ label ] = pProbe;
    if ( pBait ) _baitToProbes.insert( std::make_pair( pProbe->baitIndex(), pProbe->index() ) );

    return ( pProbe );
}

assay_index_t OPAData::addAssay(
    const assay_label_t&    label,
    const probe_label_t&    probeLabel,
    prob_t                  multiplier
){
    if ( _matrix.size2() > 0 ) {
        throw std::runtime_error( "Cannot add assay, measurements matrix already initialized" );
    }
    probe_ptr_t     pProbe = probe( probeLabel );
    if ( !pProbe ) {
        throw entity_not_found( "Probe", probeLabel );
    }
    if ( _assaysLabelMap.find( label ) != _assaysLabelMap.end() ) {
        throw entity_already_exists( "Assay", label );
    }
    assay_index_t assayIx = _assays.size();
    _assays.push_back( OPAAssay( pProbe->_index, assayIx, label, multiplier ) );
    _assaysLabelMap[ label ] = assayIx;
    pProbe->addAssay( assayIx );

    return ( assayIx );
}

void OPAData::addMeasurement(
    const object_label_t&   objectLabel,
    const assay_label_t&    assayLabel,
    celldata_t    measurement
){
    if ( _matrix.size1() != _objects.size() || _matrix.size2() != _assays.size() ) {
        initMatrixDataContainers();
    }
    const const_object_ptr_t  pObj = object( objectLabel );
    if ( !pObj ) {
        throw entity_not_found( "Object", objectLabel );
    }
    assay_index_t assayIx = assayIndex( assayLabel );
    if ( assayIx == ASSAY_NA ) {
        throw entity_not_found( "Assay", assayLabel );
    }
    _matrix( pObj->_index, assayIx ) = measurement;
    // ensure factorial is precalculated for this measurement
    ln_factorial_cache_fill( measurement.sc );
    _matrixModified = true;
}

const OPAData::celldata_t& OPAData::maxMeasurement() const
{
    const celldata_t* maxM = NULL;
    for ( object_index_t objIx = 0; objIx < _matrix.size1(); objIx++ ) {
        for ( assay_index_t aIx = 0; aIx < _matrix.size2(); aIx++ ) {
            if ( maxM == NULL || maxM->sc < _matrix( objIx, aIx ).sc ) {
                maxM = _matrix.data( objIx, aIx );
            }
        }
    }
    return ( *maxM );
}

void OPAData::initMatrixDataContainers()
{
    // initialize measurements matrix (could be done only once, otherwise current data would be lost)
    _matrix.reset( _objects.size(), _assays.size() );
    _matrixModified = true;
}

OPAData::hit_counts_type OPAData::getObjectsHitCounts( const object_set_t& objects ) const
{
    updateMatrixDependent();
    hit_counts_type res( _objectHits[ 0 ].size(), 0 );
    for ( object_set_t::const_iterator objIt = objects.begin(); objIt != objects.end(); ++objIt ) {
        const hit_mask_type& assaysMask = _objectHits[ *objIt ];
        foreach_bit( assay_index_t, assayIx, assaysMask ) {
            res[ assayIx ]++;
        }
    }
    return ( res );
}

OPAData::hit_counts_type OPAData::getProbesHitCounts( const probe_bitset_t& probes ) const
{
    updateMatrixDependent();
    hit_counts_type res( _probeHits[ 0 ].size(), 0 );
    foreach_bit( probe_index_t, probeIx, probes ) {
        const hit_mask_type& objMask = _probeHits[ probeIx ];
        foreach_bit( object_index_t, objIx, objMask ) {
            res[ objIx ]++;
        }
    }
    return ( res );
}

/**
 *  Recalculate object and probe hits caches.
 */
void OPAData::resetHits() const
{
    // initialize mask containers
    _objectHits.resize( objectsCount(), hit_mask_type( assaysCount() ) );
    _probeHits.resize( probesCount(), hit_mask_type( objectsCount() ) );
    for ( object_index_t objIx = 0; objIx < objectsCount(); objIx++ ) {
        for ( probe_index_t probeIx = 0; probeIx < probesCount(); probeIx++ ) {
            const OPAProbe& probe = _probes[ probeIx ];
            for ( size_t assayIx = probe.assayIndexes().front(); assayIx <= probe.assayIndexes().back(); ++assayIx ) {
                const celldata_t measurement = _matrix( objIx, assayIx );
                // ensure factorial is precalculated for this measurement
                ln_factorial_cache_fill( measurement.sc );
                if ( measurement.sc > 0 ) {
                    // if ever object encountered, it's set
                    _objectHits[ objIx ].set( assayIx );
                    _probeHits[ probeIx ].set( objIx );
                }
            }
        }
    }
}

void OPAData::resetIndexes()
{
    // restore label maps
    _objectLabelMap.clear();
    for ( object_index_t objIx = 0; objIx < _objects.size(); ++objIx ) {
        _objectLabelMap[ _objects[ objIx ].label() ] = &_objects[ objIx ];
    }
    // restore probes label and baits-to-probe map
    _probeLabelMap.clear();
    _baitToProbes.clear();
    for ( probe_index_t probeIx = 0; probeIx < _probes.size(); ++probeIx ) {
        OPAProbe& probe = _probes[ probeIx ];
        probe._assayIndices.clear(); 
        _probeLabelMap[ probe.label() ] = &probe;
        if ( probe.baitIndex() != OBJECT_NA ) {
            _baitToProbes.insert( std::make_pair( probe.baitIndex(), probe.index() ) );
        }
    }
    for ( assay_index_t assayIx = 0; assayIx < _assays.size(); ++assayIx ) {
        OPAAssay& assay = _assays[ assayIx ];
        _probes[ assay.probeIndex() ].addAssay( assayIx );
        _assaysLabelMap[ assay.label() ] = assayIx;
    }
}

void OPAData::resetSumMatrix() const
{
    _sumMatrix.reset( objectsCount(), probesCount(), 0 );
    for ( probe_index_t probeIx = 0; probeIx < _probes.size(); ++probeIx ) {
        const OPAProbe& probe = _probes[ probeIx ];
        const OPAAssay* pLastAssay = &assay( probe.assayIndexes().back() );
        const assay_index_t firstAssayIx = probe.assayIndexes().front();

        for ( object_index_t objIx = 0; objIx < _objects.size(); ++objIx ) {
            MassSpectraData sum = 0;
            const celldata_t*  cellDataVec = &measurement( objIx, firstAssayIx );
            for ( const OPAAssay* pAssay = &assay( firstAssayIx ); pAssay <= pLastAssay; ++pAssay ) {
                const OPAData::celldata_t&  cellData = *(cellDataVec++);
                sum.sc += cellData.sc;
                sum.pc += cellData.pc;
            }
            _sumMatrix( objIx, probeIx ) = sum;
        }
    }
}

/**
 *  Map from object index to object label.
 */
objects_label_map_type OPAData::objectsLabelMap() const
{
    objects_label_map_type  objects;
    for ( object_index_t i = 0; i < objectsCount(); i++ ) {
        objects[ i ] = object( i ).label();
    }
    return ( objects );
}

probes_label_map_type OPAData::probesLabelMap() const
{
    probes_label_map_type   probes;
    for ( probe_index_t i = 0; i < probesCount(); i++ ) {
        probes[ i ] = probe( i ).label();
    }
    return ( probes );
}

template<class Archive, class OutputStream>
void SerializeOPAData(
    OutputStream&   out,
    const OPAData&  data
){
    if ( !out.good() ) THROW_RUNTIME_ERROR( "Error writing to stream" );

    Archive dataArchive( out );
    //Rprintf( "Writing OPA data...\n" );
    dataArchive << BOOST_SERIALIZATION_NVP( data );
}

void OPAData::save(
    const char* filename
) const {
    if ( endsWith( filename, ".xml" ) ) {
        std::ofstream out( filename, std::ios_base::out );
        SerializeOPAData<boost::archive::xml_oarchive>( out, *this );
    }
    else if ( endsWith( filename, ".xml.z" ) ) {
        std::ofstream compressedOut( filename, std::ios_base::out | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::output> out;
        out.push( boost::iostreams::zlib_compressor() );
        out.push( compressedOut );
        SerializeOPAData<boost::archive::xml_oarchive>( out, *this );
    }
    else if ( endsWith( filename, ".xml.gz" ) ) {
        std::ofstream compressedOut( filename, std::ios_base::out | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::output> out;
        out.push( boost::iostreams::gzip_compressor() );
        out.push( compressedOut );
        SerializeOPAData<boost::archive::xml_oarchive>( out, *this );
    }
    else if ( endsWith( filename, ".xml.bz2" ) ) {
        std::ofstream compressedOut( filename, std::ios_base::out | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::output> out;
        out.push( boost::iostreams::bzip2_compressor() );
        out.push( compressedOut );
        SerializeOPAData<boost::archive::xml_oarchive>( out, *this );
    }
    else if ( endsWith( filename, ".bar" ) ) {
        std::ofstream out( filename, std::ios_base::out | std::ios_base::binary );
        SerializeOPAData<boost::archive::binary_oarchive>( out, *this );
    }
    else if ( endsWith( filename, ".bar.z" ) ) {
        std::ofstream compressedOut( filename, std::ios_base::out | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::output> out;
        out.push( boost::iostreams::zlib_compressor() );
        out.push( compressedOut );
        SerializeOPAData<boost::archive::binary_oarchive>( out, *this );
    }
    else if ( endsWith( filename, ".bar.gz" ) ) {
        std::ofstream compressedOut( filename, std::ios_base::out | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::output> out;
        out.push( boost::iostreams::gzip_compressor() );
        out.push( compressedOut );
        SerializeOPAData<boost::archive::binary_oarchive>( out, *this );
    }
    else if ( endsWith( filename, ".bar.bz2" ) ) {
        std::ofstream compressedOut( filename, std::ios_base::out | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::output> out;
        out.push( boost::iostreams::bzip2_compressor() );
        out.push( compressedOut );
        SerializeOPAData<boost::archive::binary_oarchive>( out, *this );
    }
    else {
        THROW_RUNTIME_ERROR( "Unsupported extension of walk file" );
    }
}

template<class Archive, class InputStream>
OPAData DeserializeOPAData(
    InputStream&                            dataStream
){
    if ( !dataStream.good() ) THROW_RUNTIME_ERROR( "Error reading OPAData stream" );

    Archive dataArchive( dataStream );

    LOG_INFO( "Reading OPAData...\n" );
    OPAData  data;
    dataArchive >> boost::serialization::make_nvp( "data", data );

    return ( data );
}

OPAData OPAData::load(
    const char*     filename
){
    LOG_INFO( "Reading OPAData file " << filename << "..." );
    if ( endsWith( filename, ".xml" ) ) {
        std::ifstream dataStream( filename, std::ios_base::in );
        return ( DeserializeOPAData<boost::archive::xml_iarchive>( dataStream ) );
    }
    else if ( endsWith( filename, ".xml.z" ) ) {
        std::ifstream walkCompressedStream( filename, std::ios_base::in | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::input> dataStream;
        dataStream.push( boost::iostreams::zlib_decompressor() );
        dataStream.push( walkCompressedStream );
        return ( DeserializeOPAData<boost::archive::xml_iarchive>( dataStream ) );
    }
    else if ( endsWith( filename, ".xml.gz" ) ) {
        std::ifstream walkCompressedStream( filename, std::ios_base::in | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::input> dataStream;
        dataStream.push( boost::iostreams::gzip_decompressor() );
        dataStream.push( walkCompressedStream );
        return ( DeserializeOPAData<boost::archive::xml_iarchive>( dataStream ) );
    }
    else if ( endsWith( filename, ".xml.bz2" ) ) {
        std::ifstream walkCompressedStream( filename, std::ios_base::in | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::input> dataStream;
        dataStream.push( boost::iostreams::bzip2_decompressor() );
        dataStream.push( walkCompressedStream );
        return ( DeserializeOPAData<boost::archive::xml_iarchive>( dataStream ) );
    }
    else if ( endsWith( filename, ".bar" ) ) {
        std::ifstream dataStream( filename, std::ios_base::in | std::ios_base::binary );
        return ( DeserializeOPAData<boost::archive::binary_iarchive>( dataStream ) );
    }
    else if ( endsWith( filename, ".bar.z" ) ) {
        std::ifstream walkCompressedStream( filename, std::ios_base::in | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::input> dataStream;
        dataStream.push( boost::iostreams::zlib_decompressor() );
        dataStream.push( walkCompressedStream );
        return ( DeserializeOPAData<boost::archive::binary_iarchive>( dataStream ) );
    }
    else if ( endsWith( filename, ".bar.gz" ) ) {
        std::ifstream walkCompressedStream( filename, std::ios_base::in | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::input> dataStream;
        dataStream.push( boost::iostreams::gzip_decompressor() );
        dataStream.push( walkCompressedStream );
        return ( DeserializeOPAData<boost::archive::binary_iarchive>( dataStream ) );
    }
    else if ( endsWith( filename, ".bar.bz2" ) ) {
        std::ifstream walkCompressedStream( filename, std::ios_base::in | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::input> dataStream;
        dataStream.push( boost::iostreams::bzip2_decompressor() );
        dataStream.push( walkCompressedStream );
        return ( DeserializeOPAData<boost::archive::binary_iarchive>( dataStream ) );
    }
    else {
        THROW_EXCEPTION( std::invalid_argument, "Unsupported extension of walk file: " << filename );
    }
}
