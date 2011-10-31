#include "BIMAPResultsSerialize.h"

#include <fstream>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include "unordered_serialization.h"
#include "misc_utils.h"

template<class Archive, class InputStream>
void DeserializeResults(
    InputStream&                resultsStream,
    ChessboardBiclusteringsIndexing&   indexing,
    objects_label_map_type&     objects,
    probes_label_map_type&      probes,
    BIMAPWalk**                pWalk,
    ChessboardBiclusteringsPDFEval** pPDFEval
){
    if ( pWalk ) {
        *pWalk = NULL;
    } else {
        throw std::invalid_argument( "No walk pointer" );
    }
    if ( !resultsStream.good() ) {
        throw std::runtime_error( "Error reading walk stream" );
    }
    Archive resultsArchive( resultsStream );

    LOG_INFO( "Reading objects and probes...\n" );
    resultsArchive >> boost::serialization::make_nvp( "objects", objects );
    resultsArchive >> boost::serialization::make_nvp( "probes", probes );

    LOG_INFO( "Reading chessboard biclusterings index...\n" );
    resultsArchive >> boost::serialization::make_nvp( "chessboardBiclusteringsIndexing", indexing );

    LOG_INFO( "Reading walk...\n" );
    *pWalk = new BIMAPWalk( indexing );
    resultsArchive >> boost::serialization::make_nvp( "walk", **pWalk );

    if ( pPDFEval ) {
      LOG_INFO( "Reading PDF adjustment...\n" );
      resultsArchive >> boost::serialization::make_nvp( "pdfAdjust", *pPDFEval );
    } else {
      LOG_INFO( "Not loading PDF adjustment...\n" );
    }
}

template<class Archive, class OutputStream>
void SerializeResults(
    OutputStream&                       resultsStream, 
    const objects_label_map_type&       objects,
    const probes_label_map_type&        probes,
    const BIMAPWalk&                   walk,
    const ChessboardBiclusteringsPDFEval*    pdfAdjust
){
    if ( !resultsStream.good() ) {
        THROW_RUNTIME_ERROR( "Error writing to walk stream" );
    }

    Archive resultsArchive( resultsStream );

    //Rprintf( "Writing objects and probes...\n" );
    resultsArchive << BOOST_SERIALIZATION_NVP( objects );
    resultsArchive << BOOST_SERIALIZATION_NVP( probes );

    //Rprintf( "Writing chessboard biclusterings index...\n" );
    resultsArchive << boost::serialization::make_nvp( "chessboardBiclusteringsIndexing", walk.indexing() );

    //Rprintf( "Writing walk...\n" );
    resultsArchive << BOOST_SERIALIZATION_NVP( walk );

    LOG_INFO( "Saving PDF adjustment...\n" );
    resultsArchive << boost::serialization::make_nvp( "pdfAdjust", pdfAdjust );
}

void BIMAPResultsLoad(
    const char*                 filename,
    ChessboardBiclusteringsIndexing&   indexing,
    objects_label_map_type&     objects,
    probes_label_map_type&      probes,
    BIMAPWalk**                pWalk,
    ChessboardBiclusteringsPDFEval** pPdfAdjust
){
    LOG_INFO( "Reading BIMAP results file " << filename << "..." );
    if ( !boost::filesystem::exists( filename ) ) {
        THROW_EXCEPTION( std::invalid_argument, "File '" << filename << "' not found" );
    }
    if ( endsWith( filename, ".xml" ) ) {
        std::ifstream resultsStream( filename, std::ios_base::in );
        DeserializeResults<boost::archive::xml_iarchive>( resultsStream, indexing,
                                                          objects, probes, pWalk, pPdfAdjust );
    }
    else if ( endsWith( filename, ".xml.z" ) ) {
        std::ifstream resultsCompressedStream( filename, std::ios_base::in | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::input> resultsStream;
        resultsStream.push( boost::iostreams::zlib_decompressor() );
        resultsStream.push( resultsCompressedStream );
        DeserializeResults<boost::archive::xml_iarchive>( resultsStream, indexing, 
                                                          objects, probes, pWalk, pPdfAdjust );
    }
    else if ( endsWith( filename, ".xml.gz" ) ) {
        std::ifstream resultsCompressedStream( filename, std::ios_base::in | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::input> resultsStream;
        resultsStream.push( boost::iostreams::gzip_decompressor() );
        resultsStream.push( resultsCompressedStream );
        DeserializeResults<boost::archive::xml_iarchive>( resultsStream, indexing,
                                                          objects, probes, pWalk, pPdfAdjust );
    }
    else if ( endsWith( filename, ".xml.bz2" ) ) {
        std::ifstream resultsCompressedStream( filename, std::ios_base::in | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::input> resultsStream;
        resultsStream.push( boost::iostreams::bzip2_decompressor() );
        resultsStream.push( resultsCompressedStream );
        DeserializeResults<boost::archive::xml_iarchive>( resultsStream, indexing,
                                                          objects, probes, pWalk, pPdfAdjust );
    }
    else if ( endsWith( filename, ".bar" ) ) {
        std::ifstream resultsStream( filename, std::ios_base::in | std::ios_base::binary );
        DeserializeResults<boost::archive::binary_iarchive>( resultsStream, indexing,
                                                             objects, probes, pWalk, pPdfAdjust );
    }
    else if ( endsWith( filename, ".bar.z" ) ) {
        std::ifstream resultsCompressedStream( filename, std::ios_base::in | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::input> resultsStream;
        resultsStream.push( boost::iostreams::zlib_decompressor() );
        resultsStream.push( resultsCompressedStream );
        DeserializeResults<boost::archive::binary_iarchive>( resultsStream, indexing,
                                                             objects, probes, pWalk, pPdfAdjust );
    }
    else if ( endsWith( filename, ".bar.gz" ) ) {
        std::ifstream resultsCompressedStream( filename, std::ios_base::in | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::input> resultsStream;
        resultsStream.push( boost::iostreams::gzip_decompressor() );
        resultsStream.push( resultsCompressedStream );
        DeserializeResults<boost::archive::binary_iarchive>( resultsStream, indexing,
                                                             objects, probes, pWalk, pPdfAdjust );
    }
    else if ( endsWith( filename, ".bar.bz2" ) ) {
        std::ifstream resultsCompressedStream( filename, std::ios_base::in | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::input> resultsStream;
        resultsStream.push( boost::iostreams::bzip2_decompressor() );
        resultsStream.push( resultsCompressedStream );
        DeserializeResults<boost::archive::binary_iarchive>( resultsStream, indexing,
                                                             objects, probes, pWalk, pPdfAdjust );
    }
    else {
        THROW_EXCEPTION( std::invalid_argument, "Unsupported extension of results file: " << filename );
    }
}

void BIMAPResultsSave(
    const char*                         filename,
    const objects_label_map_type&       objects,
    const probes_label_map_type&        probes,
    const BIMAPWalk&                   walk,
    const ChessboardBiclusteringsPDFEval*    pdfAdjust
){
    if ( endsWith( filename, ".xml" ) ) {
        std::ofstream resultsStream( filename, std::ios_base::out );
        SerializeResults<boost::archive::xml_oarchive>( resultsStream, objects, probes, walk, pdfAdjust );
    }
    else if ( endsWith( filename, ".xml.z" ) ) {
        std::ofstream resultsCompressedStream( filename, std::ios_base::out | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::output> resultsStream;
        resultsStream.push( boost::iostreams::zlib_compressor() );
        resultsStream.push( resultsCompressedStream );
        SerializeResults<boost::archive::xml_oarchive>( resultsStream, objects, probes, walk, pdfAdjust );
    }
    else if ( endsWith( filename, ".xml.gz" ) ) {
        std::ofstream resultsCompressedStream( filename, std::ios_base::out | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::output> resultsStream;
        resultsStream.push( boost::iostreams::gzip_compressor() );
        resultsStream.push( resultsCompressedStream );
        SerializeResults<boost::archive::xml_oarchive>( resultsStream, objects, probes, walk, pdfAdjust );
    }
    else if ( endsWith( filename, ".xml.bz2" ) ) {
        std::ofstream resultsCompressedStream( filename, std::ios_base::out | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::output> resultsStream;
        resultsStream.push( boost::iostreams::bzip2_compressor() );
        resultsStream.push( resultsCompressedStream );
        SerializeResults<boost::archive::xml_oarchive>( resultsStream, objects, probes, walk, pdfAdjust );
    }
    else if ( endsWith( filename, ".bar" ) ) {
        std::ofstream resultsStream( filename, std::ios_base::out | std::ios_base::binary );
        SerializeResults<boost::archive::binary_oarchive>( resultsStream, objects, probes, walk, pdfAdjust );
    }
    else if ( endsWith( filename, ".bar.z" ) ) {
        std::ofstream resultsCompressedStream( filename, std::ios_base::out | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::output> resultsStream;
        resultsStream.push( boost::iostreams::zlib_compressor() );
        resultsStream.push( resultsCompressedStream );
        SerializeResults<boost::archive::binary_oarchive>( resultsStream, objects, probes, walk, pdfAdjust );
    }
    else if ( endsWith( filename, ".bar.gz" ) ) {
        std::ofstream resultsCompressedStream( filename, std::ios_base::out | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::output> resultsStream;
        resultsStream.push( boost::iostreams::gzip_compressor() );
        resultsStream.push( resultsCompressedStream );
        SerializeResults<boost::archive::binary_oarchive>( resultsStream, objects, probes, walk, pdfAdjust );
    }
    else if ( endsWith( filename, ".bar.bz2" ) ) {
        std::ofstream resultsCompressedStream( filename, std::ios_base::out | std::ios_base::binary );
        boost::iostreams::filtering_stream<boost::iostreams::output> resultsStream;
        resultsStream.push( boost::iostreams::bzip2_compressor() );
        resultsStream.push( resultsCompressedStream );
        SerializeResults<boost::archive::binary_oarchive>( resultsStream, objects, probes, walk, pdfAdjust );
    }
    else {
        THROW_RUNTIME_ERROR( "Unsupported extension of results file" );
    }
}
