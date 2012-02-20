#include "cemm/bimap/OPADataImportCSV.h"

#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>

namespace cemm { namespace bimap {

typedef boost::escaped_list_separator<char> csv_listsep_t;
typedef boost::tokenizer< boost::escaped_list_separator<char> > csvrow_tokenizer_t;

BIMAPIOParams::BIMAPIOParams()
    : minCrossClusRefCount( 1 )
    , minObjectsPtnRefCount( 1 )
    , minProbesPtnRefCount( 1 )
    , csvColumnSeparator( '\t' )
    , mapBaitsToObjects( true )
    , objectsUniverseSize( 25000 )
{
}

void BIMAPIOParams::checkFilenames()
{
    boost::filesystem::path proteins_file_path( proteinsFilename );
    if ( !boost::filesystem::exists( proteins_file_path ) ) {
        THROW_RUNTIME_ERROR( "Proteins file not found: " << proteins_file_path );
    }
    proteinsFilename = proteins_file_path.string();
    boost::filesystem::path exp_design_file_path( expDesignFilename );
    if ( !boost::filesystem::exists( exp_design_file_path ) ) {
        THROW_RUNTIME_ERROR( "Experimental design file not found: " << exp_design_file_path );
    }
    expDesignFilename = exp_design_file_path.string();
    boost::filesystem::path measurements_file_path( measurementsFilename );
    if ( !boost::filesystem::exists( measurements_file_path ) ) {
        THROW_RUNTIME_ERROR( "Measurements file not found: " << measurements_file_path );
    }
    measurementsFilename = measurements_file_path.string();
}

csvrow_tokenizer_t RowTokenizer( const std::string& row, char sep = '\t' )
{
    return ( csvrow_tokenizer_t( row, csv_listsep_t( '\\', sep ) ) ); 
}

/**
 *  Imports proteins information from CSV file.
 *  Columns: protein AC, sequence length.
 */
void ImportProteins(
    OPAData&    data,
    const char* proteinsFilename,
    char sep = '\t'
){
    LOG_INFO( "Loading proteins from CSV file " << proteinsFilename << "..." );
    std::ifstream proteinsStream( proteinsFilename, std::ios_base::in );
    std::string row;

    // skip header
    std::getline( proteinsStream, row );
    while ( std::getline( proteinsStream, row ) ) {
        csvrow_tokenizer_t tokenizer = RowTokenizer( row, sep );
        csvrow_tokenizer_t::const_iterator tokenIt = tokenizer.begin();
        std::string ac = *(tokenIt++);
        boost::trim( ac );
        size_t      seqlen = boost::lexical_cast<size_t>( boost::trim_copy( *tokenIt ) );
        data.addObject( ac, seqlen );
    }
}

/**
 *  Imports experimental design information from CSV file.
 *  Columns: bait Accession Code, sample ID, MS run ID, MS run multiplier
 */
void ImportExpDesign(
    OPAData&    data,
    const char* expDesignFilename,
    char sep = '\t'
){
    LOG_INFO( "Loading experimental design from CSV file " << expDesignFilename << "..." );
    std::ifstream expDesignStream( expDesignFilename, std::ios_base::in );
    std::string row;

    // skip header
    std::getline( expDesignStream, row );
    while ( std::getline( expDesignStream, row ) ) {
        csvrow_tokenizer_t tokenizer = RowTokenizer( row, sep );
        csvrow_tokenizer_t::const_iterator tokenIt = tokenizer.begin();
        std::string baitAc = *(tokenIt++);
        boost::trim( baitAc );
        std::string sample = *(tokenIt++);
        boost::trim( sample );
        std::string msrun = *(tokenIt++);
        boost::trim( msrun );
        double multiplier = boost::lexical_cast<double>( boost::trim_copy( *tokenIt ) );
        OPAData::const_probe_ptr_t probe = ((const OPAData&)data).probe( sample );
        if ( !probe ) {
            // create probe if not exist yet
            // bait protein should be already defined
            probe = data.addProbe( sample, baitAc );
        } else {
            /// @todo check bait protein in CSV file matches one already defined in OPAData for given probe
        }
        data.addAssay( msrun, sample, multiplier );
    }
}

/**
 *  Imports measurements information from CSV file.
 *  Columns: MS run ID, protein AC, spectral counts, peptide counts (optional).
 */
void ImportMeasurements(
    OPAData&    data,
    const char* measurementsFilename,
    char sep = '\t'
){
    LOG_INFO( "Loading measurements from CSV file " << measurementsFilename << "..." );
    std::ifstream measurementsStream( measurementsFilename, std::ios_base::in );
    std::string row;

    // skip header
    std::getline( measurementsStream, row );
    while ( std::getline( measurementsStream, row ) ) {
        csvrow_tokenizer_t tokenizer = RowTokenizer( row, sep );
        csvrow_tokenizer_t::const_iterator tokenIt = tokenizer.begin();
        std::string msrun = *(tokenIt++);
        boost::trim( msrun );
        std::string proteinAc = *(tokenIt++);
        boost::trim( proteinAc );
        size_t sc = boost::lexical_cast<size_t>( boost::trim_copy( *tokenIt ) );
        size_t pc = tokenIt != tokenizer.end() ? boost::lexical_cast<size_t>( boost::trim_copy( *tokenIt ) ) : 0;
        OPAData::celldata_t mes( sc, pc );
        data.addMeasurement( proteinAc, msrun, mes );
    }
}

/**
 *  Imports CSV files into OPAData.
 *  @see ImportProteins()
 *  @see ImportExpDesign()
 *  @see ImportMeasurements()
 */
OPAData OPADataImportCSV(
    BIMAPIOParams&  ioParams
){
    ioParams.checkFilenames();
    OPAData data( ioParams.mapBaitsToObjects, ioParams.objectsUniverseSize );
    ImportProteins( data, ioParams.proteinsFilename.c_str(), ioParams.csvColumnSeparator );
    ImportExpDesign( data, ioParams.expDesignFilename.c_str(), ioParams.csvColumnSeparator );
    ImportMeasurements( data, ioParams.measurementsFilename.c_str(), ioParams.csvColumnSeparator );
    return ( data );
}

} }