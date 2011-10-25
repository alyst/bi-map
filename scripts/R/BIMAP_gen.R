# Generation of signals/MS data matrices for synthetic tests.
# 
# Author: Alexey Stukalov
###############################################################################

source( paste( ppi_scripts_path, "BIMAP.R", sep = .Platform$file.sep ) )
source( paste( ppi_scripts_path, "BIMAPMath.R", sep = .Platform$file.sep ) )

#' Loads cross-clusters matrix from file
#' @param matrix_file 
#' @param replicate_msruns 
#' @param multiplier.sd 
#' @returnType list of dataframes
#' @return proteins, samples, cross-clusters, signals dataframe
#' @author Alexey Stukalov
#' @export
BIMAP.read_cc_matrix <- function(
    matrix.file
){
    bimap.matrix <- read.table( matrix.file, sep = '\t', header = FALSE, stringsAsFactors = FALSE )
    print(bimap.matrix)

    make.partition <- function( elements, start_cluster_indicator, 
        element_col, cluster_col
    ){
        last_cluster <- 0
        res <- do.call( 'rbind', lapply( 1:length(elements), function (ix) {
                    if ( start_cluster_indicator[ix] ) {
                        last_cluster <<- last_cluster+1
                    }
                    data.frame(
                        cluster = last_cluster, 
                        element = elements[ix],
                        stringsAsFactors = FALSE
                    )
                } ) )
        rownames( res ) <- res$element
        colnames( res ) <- c( cluster_col, element_col )
        res
    }

    # read samples (colnames)
    samples <- as.character(bimap.matrix[1,2:ncol(bimap.matrix)])
    first.samples <- grepl('*', samples, fixed = TRUE )
    # starts of clusters are marked by *
    samples <- gsub( '*', '', samples, fixed = TRUE )
    
    # make samples data.frame
    samples.df <- do.call( 'rbind', lapply( samples, function( sample_code ) {
            sample_props <- gsub( ' ', '', unlist(strsplit( sample_code, ',' )), fixed = TRUE )
            data.frame(
                sample = sample_props[1],
                bait_ac = sample_props[2],
                stringsAsFactors = FALSE
            )
        } ) )
    rownames( samples.df ) <- samples.df$sample
    samples.clusters.df <- make.partition( samples.df$sample, first.samples,
                                           'sample', 'samples.cluster' )

    # read proteins (rownames)
    proteins <- as.character(bimap.matrix[2:nrow(bimap.matrix),1])
    # starts of clusters are marked by *
    first.proteins <- grepl('*', proteins, fixed = TRUE )
    proteins <- gsub( '*', '', proteins, fixed = TRUE )
    proteins.df <- do.call( 'rbind', lapply( proteins, function( protein_code ) {
            protein_props <- gsub( ' ', '', unlist(strsplit( protein_code, ',' )), fixed=TRUE )
            data.frame(
                protein_ac = protein_props[1],
                seqlength = as.integer( protein_props[2] ),
                multiple = ifelse( length( protein_props ) >= 3, protein_props[3], 1 ),
                stringsAsFactors = FALSE
            )
        } ) )
    rownames( proteins.df ) <- proteins.df$protein_ac
    proteins.clusters.df <- make.partition( proteins.df$protein_ac, first.proteins,
        'protein_ac', 'proteins.cluster' )

    # read cross clusters
    cross.clusters.df <- do.call( 'rbind', lapply( 2:nrow(bimap.matrix), function ( row.ix ) {
        if ( first.proteins[ row.ix - 1 ] ) {
            protein.row <- proteins.clusters.df[ row.ix - 1, ]
            return ( do.call( 'rbind', lapply( 2:ncol(bimap.matrix), function ( col.ix ) {
                    cell.val <- bimap.matrix[row.ix,col.ix]
                    if ( first.samples[ col.ix - 1 ] && !is.na( cell.val ) && nchar( cell.val ) > 0 ) {
                        sample.row <- samples.clusters.df[ col.ix - 1, ]
                        print(sample.row)
                        print(protein.row)
                        return( data.frame( proteins.cluster = protein.row$proteins.cluster,
                                            samples.cluster = sample.row$samples.cluster,
                                            signal = as.numeric( cell.val ),
                                            stringsAsFactors = FALSE
                                ) )
                    } else {
                        return ( NULL )
                    }
                } ) ) )
        } else {
            return ( NULL )
        }
    } ) )

    return ( list(
            proteins = proteins.df,
            proteins.clusters = proteins.clusters.df,
            samples = samples.df,
            samples.clusters = samples.clusters.df,
            cross.clusters = cross.clusters.df
    ) )
}


#' Generate partition as a sample from Pitman-Yor process
#' @returnType data.frame
BIMAP.generate_partition <- function(
    nElements, concentration = 0.5, discount = 0.0,
    elementPrefix = 'E'
){
    clusters <- rPitmanYor( nElements, concentration, discount )
    res <- data.frame( 
        element = paste( elementPrefix, 1:nElements, sep = '' ),
        cluster = rep( as.character(1:length( clusters )), clusters ),
        stringsAsFactors = FALSE
    )
    rownames( res ) <- res$element
    return ( res )
}

#' Generate partition with Poisson-distributed cluster sizes
#' @returnType data.frame
BIMAP.Poisson_distributed_partition <- function(
    nElements, cluster.size.rate = 4.0, cluster.size.shape = NA,
    elementPrefix = 'E'
){
    elements_not_used <- paste( elementPrefix, 1:nElements, sep = '' )
    clu_ix = 1
    res = NULL
    while ( length(elements_not_used) > 0 ) {
        clu_size <- min( 1 + ifelse( is.na( cluster.size.shape ),
                                     rpois( 1, cluster.size.rate - 1 ),
                                     rlagpois( 1, log( cluster.size.rate - 1 ), cluster.size.shape ) ),
                         length( elements_not_used ) )
        print( paste( 'Cluster #',clu_ix,', size=',clu_size, sep='' ) )
        cluster.df <- data.frame(
            element = sample( elements_not_used, size = clu_size ),
            cluster = as.character(clu_ix),
            stringsAsFactors = FALSE
        )
        clu_ix <- clu_ix + 1
        if ( length(res)>0 ) {
            res <- rbind( res, cluster.df )
        } else {
            res <- cluster.df
        }
        elements_not_used <- setdiff( elements_not_used, cluster.df$element )
    }
    rownames( res ) <- res$element
    return ( res )
}

BIMAP.generate_proteins <- function(
    nProteins,
    cluster.size.rate = NULL, cluster.size.shape = NA,
    protein.concentration = 0.5, protein.discount = 0.5,
    nProteinLengthRate = 500
){
    if ( is.null( cluster.size.rate ) ) {
        proteins.partition <- BIMAP.generate_partition( nProteins, protein.concentration, protein.discount, 'P' )
    } else {
        proteins.partition <- BIMAP.Poisson_distributed_partition( nProteins, cluster.size.rate, cluster.size.shape, 'P' )
    }
    colnames( proteins.partition ) <- c( 'protein_ac', 'proteins.cluster' )
    nProtClusters <- length( unique( proteins.partition$proteins.cluster ) )
    
    proteins.df <- data.frame(
            protein_ac = proteins.partition$protein_ac,
            seqlength = rpois( nProteins, nProteinLengthRate ),
            multiple = 1,
            stringsAsFactors = FALSE
        )
    rownames( proteins.df ) <- proteins.df$protein_ac

    return ( list( proteins = proteins.df,
                   proteins.clusters = proteins.partition
    ) )
}

# convert cross-cluster signals to a matrix
BIMAP.cross_clusters.to.signals_matrix <- function( cc_df, proteinsClusters, samplesClusters )
{
    signals.mean <- matrix( rep( NA, length( proteinsClusters ) * length( samplesClusters ) ), 
                                 nrow = length( proteinsClusters ) )
    rownames( signals.mean ) <- proteinsClusters
    colnames( signals.mean ) <- samplesClusters
    apply( cc_df, 1, function ( cc ) {
        signals.mean[ cc[['proteins.cluster']], cc[['samples.cluster']] ] <<- as.numeric( cc[['signal']] )
    } )
    return ( signals.mean )
}

BIMAP.generate_clustering <- function(
    nProteins, nSamples,
    protein.concentration = 0.5, protein.discount = 0.5,
    sample.concentration = 0.5, sample.discount = 0.3,
    nProteinLengthRate = 500,
    bimap.prob = 0.1, signal.mean = 0.0, signal.sd = 1.0
){
    # generate proteins and samples partitions
    proteins.info <- BIMAP.generate_proteins( nProteins,
                protein.concentraion = protein.concentration,
                protein.discount = protein.discount,
                nProteinLengthRate = nProteinLengthRate )
    proteins.partition <- proteins.info$proteins.clusters
    proteins.df <- proteins.info$proteins
    samples.partition <- BIMAP.generate_partition( nSamples, sample.concentration, sample.discount, 'S' )
    colnames( samples.partition ) <- c( 'sample', 'samples.cluster' )
    nSampleClusters <- length( unique( samples.partition$samples.cluster ) )

    # generate cross-clusters
    # at least one cc per sample cluster
    bimap.df <- data.frame(
        proteins.cluster = sample( nProtClusters, nSampleClusters, replace = TRUE ),
        samples.cluster = 1:nSampleClusters )
    # at least one cc per protein cluster
    bimap.df <- rbind( bimap.df, data.frame(
        proteins.cluster = 1:nProtClusters,
        samples.cluster = sample( nSampleClusters, nProtClusters, replace = TRUE ) ) )
    # turn on every cc with bimap.prob
    bimap.all.df <- merge( data.frame( proteins.cluster = 1:nProtClusters ),
                     data.frame( samples.cluster = 1:nSampleClusters ), all = TRUE )
    bimap.on.df <- bimap.all.df[ sample( nrow(bimap.all.df), bimap.prob * nrow(bimap.all.df) ), ]
    bimap.df <- unique( rbind( bimap.df, bimap.on.df ) )
    # generate cc signals
    bimap.df$signal <- rnorm( nrow(bimap.df), mean = signal.mean, sd = signal.sd )

    # map bait of clusters, protein from cc with strongest signal becomes bait
    samples.df <- data.frame( sample = samples.partition$sample,
                              bait_ac = NA,
                              stringsAsFactors = FALSE )
    rownames( samples.df ) <- samples.df$sample
    for ( sample in samples.df$sample ) {
        baits <- unique( subset( samples.df, !is.na(bait_ac) )$bait_ac )
        s.clu <- samples.partition[ sample, 'samples.cluster' ]
        sample.cc <- subset( bimap.df, samples.cluster == s.clu )
        cosamples <- subset( samples.partition, samples.cluster == s.clu )$sample
        cobaits <- subset( samples.df, sample %in% cosamples )$bait_ac
        # count number of baits per cluster
        cobait.clusters <- table( proteins.partition[ cobaits, 'proteins.cluster' ] )
        # build table of protein cluster scores (signal + number of baits)
        clu.scores <- sample.cc$signal
        names( clu.scores ) <- sample.cc$proteins.cluster
        # proteins not used as baits yet per cluster
        navail.proteins <- sapply( names( clu.scores ), function( clu.id ) {
                nrow( subset( proteins.partition, !(protein_ac %in% baits) & proteins.cluster == clu.id ) )
            } )
        names( navail.proteins ) <- names( clu.scores )
        print( navail.proteins )
        clu.scores[ names( cobait.clusters ) ] <- clu.scores[ names( cobait.clusters ) ] +
                                                  cobait.clusters
        clu.scores <- ( clu.scores - min(clu.scores) + 0.1 ) * ( navail.proteins + 1E-3)^1.5 
        print(clu.scores)
        # randomly select protein cluster with probability proportional to score
        p.clu <- sample( names(clu.scores), 1, prob = clu.scores )
        # randomly select protein from the cluster, preferably one which was not used as bait before
        avail.proteins <- subset( proteins.partition, proteins.cluster == p.clu )$protein_ac
        p.scores <- table(avail.proteins) - 1E-5 # -1E-5 is a hack to make it numeric
        p.scores[ intersect( cobaits, names( p.scores ) ) ] <- 0.1
        print(p.scores)
        bait <- sample( names(p.scores), 1, prob = p.scores )
        samples.df[ sample, 'bait_ac' ] <- bait
    }

    return ( list( proteins = proteins.df,
                   proteins.clusters = proteins.partition,
                   samples.clusters = samples.partition,
                   samples = samples.df,
                   cross.clusters = bimap.df,
                   signals.mean = BIMAP.cross_clusters.to.signals_matrix( bimap.df, 1:nProtClusters, 1:nSampleClusters )
    ) )
}

BIMAP.generate.msruns <- function(
    samples, nreplicates = 1, lnmultiplier.sd = 0.5
){
    res <- data.frame(
        msrun = paste( rep( samples$sample, each = nreplicates ),
                    rep( 1:nreplicates, nrow(samples) ), sep = '-' ),
        sample = rep( samples$sample, each = nreplicates ),
        multiplier = exp( rnorm( nrow(samples)*nreplicates,
                            mean = 0, sd = lnmultiplier.sd ) ),
                    stringsAsFactors = FALSE )
    rownames( res ) <- res$msrun
    return ( res )
}

BIMAP.generate.ms_data <- function(
    proteins.clusters, samples.clusters, cross.clusters,
    proteins, msruns, samples = NULL,
    seq.length.factor = 0.5,
    signal.shape = 0.1, noise.rate = 1E-3,
    rate.factor = 3
){
    res <- ddply( proteins, c('protein_ac'), function( protein ) {
        protsClu <- subset( proteins.clusters, protein_ac == protein$protein_ac )$proteins.cluster
        # signal factor of the protein
        protein.factor.log <- log( protein$seqlength ) * seq.length.factor + log( protein$multiple )

        do.call('rbind', lapply( unique( samples.clusters$samples.cluster ), function( samplesClu ) {
            cluSamples <- subset( samples.clusters, samples.cluster == samplesClu )$sample # samples of cluster
            cluMsRuns <- subset( msruns, sample %in% cluSamples ) # msruns of cluster
            crossClu <- subset( cross.clusters, proteins.cluster == protsClu & samples.cluster == samplesClu )
            if ( nrow( crossClu ) == 0 ) {
                # no signal, generate noise
                counts <- rgeom( nrow(cluMsRuns), 1-noise.rate )
            }
            else {
                # signal
                counts <- sapply( cluMsRuns$multiplier, function( mult ) {
                            lnRate = crossClu$signal + log( mult )
                            ifelse( is.na(rate.factor) || rbinom( 1, 1, max(0.01, 1 - exp( -rate.factor* exp(lnRate) )) ),
                                    rlagpois( 1, lnRate + protein.factor.log, signal.shape ),
                                    0 ) } )
            }
            data.frame(
                msrun = cluMsRuns$msrun,
                sc = counts,
                pc = rep( 1, nrow(cluMsRuns) ),
                stringsAsFactors = FALSE )
        } ) )
    } )
    # fix column names and order
    res <- subset( res, sc > 0 )
    colnames( res ) <- c( 'prey_ac', 'msrun', 'sc', 'pc' )
    res$sample <- msruns[ res$msrun, 'sample' ]
    if ( !is.null( samples ) ) {
        res$bait_ac <- samples[ res$sample, 'bait_ac' ]
    }
    return ( subset( res, sc > 0 ) )
}
