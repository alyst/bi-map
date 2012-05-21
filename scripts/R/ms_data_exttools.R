# Support for external MS clustering and filtering tools.
# 
# Author: Alexey Stukalov (astukalov@cemm.oeaw.ac.at)
#
###############################################################################

# Read NestedCluster results
import.nestedcluster <- function(
    path_prefix,
    data_filename,
    max.clusterings = NA,
    min.intensity = 0.1,
    prior.proteins.clusterings = FALSE
){
    # original bait x prey abundance matrix
    message('Reading Input data...')
    dheader = read.delim( file.path( path_prefix, data_filename ),
                          header=FALSE, sep="\t", as.is=TRUE, skip=0, nrows = 2)
    d = read.delim( file.path( path_prefix, data_filename ),
                    header=FALSE, sep="\t", as.is=TRUE, skip=2)
    colnames(d) <- dheader[2,]
    # first column is protein names (Experiment is description for corresponding column)
    rownames( d ) <- d$Experiment
    d$Experiment <- NULL

    message('Reading Bait Clusters...')
    # bait clusterings
    # colnames: logLik, bait IDs
    # rows, sorted by likelihood: likelihood, indices of bait clusters
    # nrows = number of records
    clus.header = read.delim( file.path( path_prefix, "Clusters" ),
                       header=FALSE, sep="\t", as.is=TRUE, skip=0, nrows=1 )
    clus = read.delim( file.path( path_prefix, "Clusters" ),
        header=FALSE, sep="\t", as.is=TRUE, skip=1 )
    colnames( clus ) <- clus.header
    nclusterings <- nrow( clus )
    message( "Total ", nclusterings, " clusterings" )
    clusterings = data.frame(
        clustering = 1:nrow( clus ),
        llh = as.numeric( clus[,1] ),
        stringsAsFactors = FALSE )
    rownames( clusterings ) <- clusterings$clustering
    baits.clusters = data.frame(
        clustering = rep( 1:nclusterings, each = ncol(clus)-1),
        bait_ac = rep( colnames(clus)[-1], times = nclusterings ),
        baits.cluster = as.vector( t( clus[,-1] ) ),
        stringsAsFactors = FALSE
    )

    # prey clusterings
    # integer matrix
    # columns: max.possible bait clusters (nb parameter)
    #          column position = baits cluster index
    # rows = (number of records, ordered by likelihood) x prey
    # cells: index of prey(row) cluster of corresponding bait(col) from given record(row - prey)
    message('Reading Nested clusters...')
    nestedclusters = read.delim( file.path( path_prefix, "NestedClusters" ),
        header=FALSE, sep="\t", as.is=TRUE )
    max.baits.clusters <- ncol( nestedclusters )
    message( nrow(nestedclusters), " rows in NestedClusters table" )
    message( ncol(nestedclusters), " max possible bait clusters" )
    proteins.clusters <- data.frame(
            clustering = rep( 1:nclusterings, each = nrow(d) * max.baits.clusters ),
            baits.cluster = rep( 1:max.baits.clusters, times = nrow(d) * nclusterings ),
            proteins.cluster = as.vector( t( nestedclusters ) ), # row-wise conversion to vector
            protein_ac = rep( rownames(d), each = max.baits.clusters, times = nclusterings ),
            stringsAsFactors = FALSE
     )

    # prey clusters intensities
    # numeric matrix
    # columns: max.possible bait clusters (nb parameter)
    # rows = (number of records, ordered by likelihood) x prey
    # cells: abundance of prey(row) in cluster of corresponding bait(col) from given record(row)
    message('Reading Phi (bicluster intensity)...')
    preys.phi = as.matrix( read.delim( file.path( path_prefix, "NestedMu" ),
        header=FALSE, sep="\t", as.is=TRUE ) )
    message( nrow(preys.phi), ' Phi values' )
    proteins.clusters$intensity = as.vector( t( preys.phi ) )
    # all clusters below intensity threshold are considered empty
    if ( is.numeric( min.intensity ) ) {
        proteins.clusters$proteins.cluster[ proteins.clusters$intensity < min.intensity ] <- -1
    }

    # prey clusters variation
    # numeric matrix
    # columns: max.possible bait clusters (nb parameter)
    # rows = (number of records, ordered by likelihood) x prey
    # cells: variation of abundance of prey(row) in cluster of corresponding bait(col) from given record(row)
    message('Reading Sigma^2 (bicluster intensity variation)...')
    preys.sigma = as.matrix( read.delim( file.path( path_prefix, "NestedSigma2" ),
        header=FALSE, sep="\t", as.is=TRUE ) )
    message( nrow(preys.sigma), ' Sigma^2 values' )
    proteins.clusters$sigma= as.vector( t( preys.sigma ) )

    if ( !prior.proteins.clusterings ) {
        # remove preys clusterings not attached to baits cluster
        nclus.old <- nrow( proteins.clusters )
        baits.clusters.ids <- unique( baits.clusters[ ,c( 'baits.cluster', 'clustering') ] )
        proteins.clusters <- merge( baits.clusters.ids, proteins.clusters,
            by = c( 'baits.cluster', 'clustering' ),
            all.x = FALSE, all.y = FALSE )
        message( "Removed prior-only prey cluster records: before=", nclus.old, " after=", nrow(proteins.clusters), sep = '' )
    }
    if ( !is.na( max.clusterings ) ) {
        clusterings <- subset( clusterings, clustering <= max.clusterings )
        baits.clusters <- subset( baits.clusters, clustering <= max.clusterings )
        proteins.clusters <- subset( proteins.clusters, clustering <= max.clusterings )
    }

    # MCMC parameters
    # iteration, likelihood, number of bait clusters, number of prey clusters for each bait
    message('Reading MCMC walk...')
    mcmc.params.file <- file.path( path_prefix, "MCMCparameters" )
    if ( file.exists( mcmc.params.file ) ) {
        mcmc.params <- read.delim( mcmc.params.file,
                                   header=FALSE, sep="\t", as.is=TRUE )
        colnames( mcmc.params ) <- c( 'step', 'llh', 'nbait.clusters', 'move.type', colnames(clus)[-1] )
    } else {
        mcmc.params <- NA
    }

    # optimal clusters
    # number of files = 2x number of records
    # optcluN: record blocks delimeterd by UNIQID,
    #          baits list
    # prey per line: prey's cluster signal
    # optclu-mu: same as above, but data precedes each mean

    return ( list(
            data = d,
            clusterings = clusterings,
            baits.clusters = baits.clusters,
            proteins.clusters = proteins.clusters,
            walk = mcmc.params
        )
    )
}

#' Multiple clusters intersections.
#' Intersects all clusters of all partitions
#' and return the resulting intersections as a new partition.
#' @param partition dataset of partitions
#' @param ptn_id_col name of partition ID column
#' @param clu_id_col name of cluster ID column
#' @param elem_id_col name of element ID column
#' @returnType 
#' @return 
#' @author astukalov
#' @export
partition.multisect <- function( partition, ptn_id_col, clu_id_col, elem_id_col )
{
    ptn = partition[, c(ptn_id_col,clu_id_col,elem_id_col) ]
    colnames(ptn) <- c( 'ptn_id', 'clu_id', 'elem_id')
    # clustering signa
    ptn.isect <- ddply( ptn, .(elem_id), function( elm.clus ) {
        data.frame( clu_id.list = paste(
                            elm.clus[order(elm.clus$ptn_id),'ptn_id'],
                            elm.clus[order(elm.clus$ptn_id),'clu_id'],
                            sep='::', collapse=' ' ),
                    stringsAsFactors = FALSE )
    } )
    # convert strings to numbers
    ptn.isect$isect.clu_id <- as.character( as.integer( as.factor( ptn.isect$clu_id.list ) ) )
    colnames( ptn.isect ) <- c( elem_id_col,
                                paste( clu_id_col, 'list', sep = '.' ),
                                paste( clu_id_col, 'isect', sep = '.' ) )
    return ( ptn.isect )
}

nestedcluster.extract <- function( clusterings, clustering_id )
{
    return ( list( data = clusterings$data,
                   clusterings = subset( clusterings$clusterings, clustering == clustering_id ),
                   baits.clusters = subset( clusterings$baits.clusters, clustering == clustering_id ),
                   proteins.clusters = subset( clusterings$proteins.clusters, clustering == clustering_id )
           ) )
}

#' Convert NestedCluster biclustering to Chessboard Biclustering
#' @param nested.clustering 
#' @param proteins 
#' @param samples 
#' @param min.intensity signal/noise cutoff threshold
#' @returnType list of data.frames compatible with chessboard biclustering representation
#' @return corresponding chessboard biclustering
#' @author astukalov
#' @see import.nestedcluster()
#' @export
nestedcluster.to.bimap <- function(
    nested.clustering,
    proteins,
    samples,
    min.intensity = 0.01
){
    # compose proteins.clusters -- intersect all nested clusters
    preys.isect.clusters <- partition.multisect( nested.clustering$proteins.clusters,
                                                 'baits.cluster', 'proteins.cluster', 'protein_ac' )
    proteins.clusters <- preys.isect.clusters[,c('proteins.cluster.isect', 'protein_ac')]
    colnames(proteins.clusters) <- c('proteins.cluster','protein_ac')
    proteins.clusters$proteins.cluster <- as.character( proteins.clusters$proteins.cluster )

    # compose samples clusters
    samples$bait_ac <- as.character( samples$bait_ac )
    samples.clusters <- merge( samples, nested.clustering$baits.clusters, by = 'bait_ac' )[,c('sample','baits.cluster')]
    colnames( samples.clusters ) <- c('sample','samples.cluster' )
    samples.clusters$samples.cluster <- as.character( samples.clusters$samples.cluster )

    # merge nested clusters with intersected prey clusters to get chessboard biclustering blocks
    # FIXME: non-unique intensity seems to be internal nestedcluster error
    preys2isect <- ddply( unique( merge( nested.clustering$proteins.clusters,
                                       preys.isect.clusters, by = 'protein_ac' )[
        , c('baits.cluster','proteins.cluster','proteins.cluster.isect','intensity')] ),
        c('baits.cluster','proteins.cluster','proteins.cluster.isect'), function( rows ) {
            if ( nrow( rows ) > 1 ) warning( 'Non unique row ', length( rows ) )
            if ( length( unique( rows$intensity) ) > 1 ) warning( 'Non unique intensitiy ', sort(unique( rows$intensity )) )
            data.frame( intensity = mean( rows$intensity ), stringsAsFactors = FALSE )
        } )
    # remove off-blocks by intensity threshold
    preys2isect <- subset( preys2isect, intensity >= min.intensity )
    preys2isect$baits.cluster <- as.character( preys2isect$baits.cluster )
    blocks <- preys2isect[,c('baits.cluster','proteins.cluster.isect','intensity','proteins.cluster')]
    colnames(blocks) <- c('samples.cluster','proteins.cluster','signal','nested.preys.cluster')
    return ( list(
            samples = samples,
            proteins = proteins,
            samples.clusters = samples.clusters,
            proteins.clusters = proteins.clusters,
            blocks = blocks,
            signals.mean = daply( blocks, .( proteins.cluster, samples.cluster ), function( data ) {
                        return ( mean(data$signal) )
                    } ),
            signals.sd= daply( blocks, .( proteins.cluster, samples.cluster ), function( data ) {
                        return ( mean( data$sigma ) )
                    } )
        )
    )
}

# Export data to NestedCluster format
ms_data.export_nestedcluster <- function(
    ms_data, protein_info, export_file, nsaf_avg = 10,
    bait_ac_col = 'bait_ac', sample_col = 'sample', prey_ac_col = 'prey_ac', sc_col = 'sc'
){
    export.df <- data.frame(
        bait_ac = ms_data[,bait_ac_col],
        sample = ms_data[,sample_col],
        prey_ac = ms_data[,prey_ac_col],
        sc_factor = ms_data[,sc_col] / protein_info[ms_data[,prey_ac_col], 'seqlength'],
        stringsAsFactors = FALSE )
    sc_factor.df <- ddply( export.df, c( 'bait_ac', 'sample' ), function ( sample_sc_factors ) {
            data.frame( mean_sc_factor = mean( sample_sc_factors$sc_factor ),
                stringsAsFactors = FALSE )
        } )
    rownames( sc_factor.df ) <- sc_factor.df$sample

    export.df$nsaf <- nsaf_avg * export.df$sc_factor /
        sc_factor.df[ export.df$sample, 'mean_sc_factor' ]
    
    nsaf_table <- reshape( export.df,
        v.names = c('nsaf'),
        timevar = c( 'sample' ),
        idvar = c( 'prey_ac'),
        direction = 'wide',
        drop = c( 'bait_ac', 'sc_factor' ) )
    rownames(nsaf_table) <- nsaf_table$prey_ac
    nsaf_table$prey_ac <- NULL
    # substitute NA with 0
    for ( sample_col in colnames( nsaf_table ) ) {
        nsaf_table[ is.na( nsaf_table[ , sample_col ] ), sample_col ] <- 0 
    }
    # prepare header
    sample_order <- sapply( colnames(nsaf_table),
        function( colname ) strsplit( colname, '.', fixed = TRUE )[[1]][[2]] )
    baits_order <- sc_factor.df[ sample_order, 'bait_ac' ]
    header <- matrix( c( baits_order, sample_order ),
        nrow = 2, ncol = length( sample_order ), byrow = TRUE )
    rownames(header) <- c('Bait','Experiment')
    #print( header )

    # write header and data
    write.table( header, file = export_file, sep = '\t',
        quote = FALSE, row.names = TRUE, col.names = FALSE )
    write.table( nsaf_table, file = export_file, sep = '\t',
        quote = FALSE, row.names = TRUE, col.names = FALSE, append = TRUE )
    return ( nsaf_table )
}

# export MS data to SAINT v.2+
ms_data.export_saint.v2 <- function(
    ms_data, protein_info, export.path,
    bait_col = 'bait_ac', prey_col = 'prey_ac',
    neg.ctrl.baits = c( 'nobait', 'gfp' )
){
    preys <- unique(ms_data[,prey_col])
    saint.data <- list( 
        preys = data.frame(
            protein_ac = preys,
            seqlength = protein_info[preys, 'seqlength' ],
            stringsAsFactors = FALSE ),
        baits = unique( subset( ms_data, select = c( 'msrun', bait_col ), msrun != '#glob#' ) ),
        ms_data = subset( ms_data, select = c( 'msrun', bait_col, prey_col, 'sc' ), msrun != '#glob#' )
    )
    saint.data$baits$type <- ifelse( saint.data$baits$bait_ac %in% neg.ctrl.baits, 'C', 'T' )
    write.table( saint.data$preys, file.path( export.path, 'preys.csv' ),
        sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE )
    write.table( saint.data$baits, file.path( export.path, 'baits.csv' ),
        sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE )
    write.table( saint.data$ms_data, file.path( export.path, 'ms_data.csv' ),
        sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE )

    return ( saint.data )
}

# export MS data to SAINT v.1+
ms_data.export_saint.v1 <- function(
    ms_data, protein_info, export.filename,
    bait_col = 'bait_ac', prey_col = 'prey_ac',
    data_col = 'sc',
    neg.ctrl.baits = c( 'nobait', 'gfp' ),
    contaminants.extra = NULL
){
    # prewide table with protein information
    ms_data$ext_id <- paste( 
        ms_data[,bait_col],
        ms_data$sample_origin,
        ms_data$msrun, sep = '.' )
    ms_data.prewide <- subset( ms_data, msrun != '#glob#')
    ms_data.prewide <- subset( ms_data.prewide, select = c(prey_col, 'ext_id', data_col ) )
    ms_data.prewide[,c('protein_name' , 'seqlength' , 'description')] <- protein_info[ 
            ms_data.prewide$prey_ac,
            c( 'protein_id' , 'seqlength' , 'description' ) ]

    msruns.df <- unique( subset( ms_data, msrun != '#glob#', select = c( 'ext_id', 'msrun', bait_col, 'sample_origin', 'msrun_origin' ) ) )
    rownames( msruns.df ) <- msruns.df$ext_id
    msruns.df <- msruns.df[ order(msruns.df$ext_id), ]
    msruns.df$bait.seqcov <- sapply( msruns.df$ext_id, function( msrun_ext_id ) {
        curbait_ac <- uniprotac.removeiso( msruns.df[ msrun_ext_id, bait_col ] )
        curbait_ac <- ifelse( curbait_ac == 'gfp', 'mgtagzh', curbait_ac ) # fix due to buggy mgtagzh sequence
        bait.ms_data <- subset( ms_data, uniprotac.removeiso( prey_ac ) == curbait_ac & ext_id == msrun_ext_id )
        if ( nrow( bait.ms_data ) > 0 ) {
            return ( max( bait.ms_data$seqcov ) )
        } else if ( curbait_ac == 'nobait' ){
            return ( 0.5 )
        } else {
            warning( 'No bait ', curbait_ac, ' found in pulldown ', msrun_ext_id )
            return ( 0 )
        }
    } )
    
    ms_data.wide <- reshape( ms_data.prewide, timevar = c( 'ext_id' ), 
        idvar = c( prey_col ), 
        ids = sort( unique( as.character( ms_data.prewide[,c(prey_col)] ) ) ),
        times = rownames( msruns.df ),
        v.names = c( data_col ),
        direction = 'wide', sep = '.' )
    ms_data.wide <- ms_data.wide[
            order( ms_data.wide$protein_name, ms_data.wide$prey_ac ), ]

    out <- file( export.filename, 'w' )

    baits <- unique( uniprotac.removeiso( ms_data$bait_ac ) )
    contaminants <- unique( uniprotac.removeiso( subset( ms_data, bait_ac %in% neg.ctrl.baits )$bait_ac ) )
    if ( !is.null(contaminants.extra) ) contaminants <- union( contaminants, unique( uniprotac.removeiso( contaminants.extra ) ) )
    # fill addition columns with neutral information
    ms_data.wide$prey_abundance <- 1 # should get it from PepAtlas
    ms_data.wide$prey_type <- ifelse( uniprotac.removeiso( ms_data.wide$prey_ac ) %in% contaminants, 'C',
                                      ifelse( uniprotac.removeiso( ms_data.wide$prey_ac ) %in% baits, 'R', 'N' ) )
    prey_info_cols <- c( 'prey_ac', 'prey_abundance', 'seqlength', 'prey_type' )
    ms_data.wide.saint <- ms_data.wide[ , c( prey_info_cols, setdiff( colnames( ms_data.wide ),
                                             c( prey_info_cols, 'protein_name', 'description') ) ) ]

    # write header
    # experiment IDs
    cat( '?\t?\t?\t?\t', file=out )
    cat( paste( msruns.df$msrun, collapse = '\t' ), file=out )
    cat( '\n', file=out)

    # bait ACs
    cat( '?\t?\t?\t?\t', file=out )
    cat( paste( msruns.df$bait_ac, collapse = '\t' ), file=out )
    cat( '\n', file=out)

    # bait coverage
    cat( '?\t?\t?\t?\t', file=out )
    cat( paste( msruns.df$bait.seqcov, collapse = '\t' ), file=out )
    cat( '\n', file=out )

    write.table( ms_data.wide.saint, quote = FALSE, file = out,
                 sep='\t', na = '0', row.names = FALSE, col.names = FALSE )
    close(out)

    return ( ms_data.wide )
}
