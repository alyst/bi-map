# BIMAP-related math functions.
# 
# Author: Alexey Stukalov (astukalov@cemm.oeaw.ac.at)
###############################################################################

require( 'Rcpp' )

# load RBIMAP library
# Note: RBIMAP.libpath should be set before executing this script, e.g.:
#       RBIMAP.libpath <- file.path( bimap_scripts_path, "build/release/src/R" )
try( dyn.unload( file.path( RBIMAP.libpath, paste("libRBIMAP-math", .Platform$dynlib.ext, sep="")) ) )
dyn.load( file.path( RBIMAP.libpath, paste("libRBIMAP-math", .Platform$dynlib.ext, sep="")), 
          type = "Call" ) 

dlagpois <- function( x, lnRate, shape, log = FALSE )
{
    res <- .Call( "LagrangianPoissonLnPdf", x, lnRate, shape )
    return ( ifelse( log, res, exp( res ) ) )
}

plagpois <- function( x, lnRate, shape, upper = FALSE, log = FALSE )
{
    fname <- ifelse( upper, "LagrangianPoissonLnCdfQ", "LagrangianPoissonLnCdfP" )
    res <- .Call( fname, x, lnRate, shape )
    return ( ifelse( log, res, exp( res ) ) )
}

rlagpois <- function( size, lnRate, shape )
{
    return ( .Call( "LagrangianPoissonRandom", size, lnRate, shape ) )
}

rPitmanYor <- function( nElements, concentration = 0.5, discount = 0.0 )
{
    return ( .Call( "PitmanYorRandomSample", nElements, concentration, discount ) )
}

countPartitionMismatches <- function(
    templatePartitionCollection,
    partitionCollection,
    templatePartitionId_col = 'template.partition',
    partitionId_col = 'partition',
    clusterId_col = 'cluster',
    elementId_col = 'element'
){
    return ( within( .Call( "CountPartitionMismatches",
                    templatePartitionCollection, partitionCollection, 
                    templatePartitionId_col, partitionId_col,
                    clusterId_col, elementId_col ), {
        exp.co.tco <- co.co * tco.tco / totalPairs
        rand.index <- ( co.tco + mismatch.tmismatch ) / totalPairs
        adj.rand.index <- ( co.tco - exp.co.tco ) / ( tco.tco - exp.co.tco )
    } ) )
}

countPartsMismatches <- function(
    templatePartitionCollection,
    partitionCollection,
    templatePartitionId_col = 'template.partition',
    partitionId_col = 'partition',
    clusterId_col = 'cluster',
    elementId_col = 'element'
){
    return ( .Call( "CountPartsMismatches",
                    templatePartitionCollection, partitionCollection, 
                    templatePartitionId_col, partitionId_col,
                    clusterId_col, elementId_col ) )
}

mutualInformation <- function(
        templatePartition,
        partition,
        clusterId_col = 'cluster',
        elementId_col = 'element',
        normalize = TRUE
){
    templatePtnSizes <- table( templatePartition[,clusterId_col] )
    ptnSizes <- table( partition[,clusterId_col] )
    ptnMerged <- merge( templatePartition, partition, by = elementId_col,
                        suffixes = c(".t", "" ) )
    ptnMerged$clu_isect <- paste( ptnMerged[,paste(clusterId_col,'t',sep='.')],
                                  ptnMerged[,clusterId_col] )
    ptnMerged <- ddply( ptnMerged, .( clu_isect ), function( cluIsect ) {
        data.frame( size.isect = nrow( cluIsect ),
                    size.t =  templatePtnSizes[ as.character(cluIsect[1, paste(clusterId_col,'t',sep='.')]) ],
                    size =  ptnSizes[ as.character(cluIsect[1, clusterId_col]) ],
                    stringsAsFactors = FALSE )
        } )
    res <- sum( ptnMerged$size.isect * ( log( ptnMerged$size.isect ) - log(ptnMerged$size.t) - log( ptnMerged$size ) ) )  / sum( ptnMerged$size.isect ) +
           log( sum( ptnMerged$size.isect) )
    if ( normalize ) {
        ent <- -sum( log( ptnSizes ) * ptnSizes ) / sum( ptnSizes ) + log( sum( ptnSizes ) ) 
        ent.t <- -sum( log( templatePtnSizes ) * templatePtnSizes ) / sum( templatePtnSizes ) + log( sum( templatePtnSizes ) )
        res <- res / ( ent + ent.t )
        print( min( ent, ent.t ) / ( ent + ent.t ) )
    }
    return ( res )
}
