# Exporting PSI bbb results to graphml
# 
# Author: astukalov
###############################################################################

require( plyr )

source( file.path( bimap_scripts_path, "graphml_support.R" ) )

BIMAP.graphML.dataframe <- function( bimap.props,
    bimap.data,
    protein_export_cols = c( name = 'protein_ac' ),
    sample_export_cols = c( name = 'sample' ),
    proteins_cluster_label_col = 'proteins.cluster',
    proteins_cluster_export_cols = c( 'avg.pairs.freq' ),
    min.intensity = 5, hex.intensity.power = 5
){
    if ( is.null( names( protein_export_cols ) ) ) {
        names( protein_export_cols ) <- protein_export_cols
    }
    protein_label_col <- ifelse( any( 'name' %in% names( protein_export_cols ) ),
                                 protein_export_cols[[ 'name' ]], 'protein_ac' )
    if ( is.null( names( sample_export_cols ) ) ) {
        names( sample_export_cols ) <- sample_export_cols
    }
    message( 'Labeling protein nodes using: ', protein_label_col )
    sample_label_col <- ifelse( any( 'name' %in% names( sample_export_cols ) ),
                                 sample_export_cols[[ 'name' ]], 'sample' )
    message( 'Labeling sample nodes using: ', sample_label_col )

    # add to sample clusters dataframe
    bimap.props$samples.clusters$bait_ac <- bimap.data$samples[ bimap.props$samples.clusters$sample, 'bait_ac' ]
    bimap.props$signal.max <- max( bimap.props$signals.mean, na.rm = TRUE )
    bimap.props$signal.min <- min( bimap.props$signals.mean, na.rm = TRUE )

    # create fake protein clusters and protein information for baits that are not proteins
    alien.baits <- setdiff( bimap.props$samples.clusters$bait_ac,
                            bimap.props$proteins.clusters$protein_ac )
    if ( length( alien.baits ) > 0 ) {
        # protein clusters that match sample clusters for "alien" baits
        fake_protein_clusters <- subset( bimap.props$samples.clusters, bait_ac %in% alien.baits,
                                         select = c( 'samples.cluster', 'bait_ac' ) )
        colnames( fake_protein_clusters ) <- c( 'proteins.cluster', 'protein_ac' )
        fake_protein_clusters$proteins.cluster <- paste( 'bait_clu_', fake_protein_clusters$proteins.cluster, sep = '' )
        bimap.props$proteins.clusters <- rbind( bimap.props$proteins.clusters,
                                                fake_protein_clusters )
		# fake protein information for baits that are not proteins
        bait.proteins.template <- bimap.data$proteins[ 1, , drop=FALSE ]
        bait.proteins<- do.call( 'rbind', lapply( alien.baits, function( bait_ac ) {
            res <- bait.proteins.template
            res[,intersect(colnames(bimap.data$proteins),proteins_cluster_export_cols)] <- NA
            res[,protein_label_col] <- bait_ac
            return ( res )
        } ) )
        rownames( bait.proteins ) <- alien.baits
        bimap.data$proteins <- rbind( bimap.data$proteins, bait.proteins )
    }

    # merge samples and proteins clusters via bait-prey
    bait_prey.clusters <- merge( bimap.props$proteins.clusters, bimap.props$samples.clusters,
                                 by.x = 'protein_ac', by.y = 'bait_ac',
                                 all.x = TRUE, all.y = FALSE )
    # merge sample cluster ids for those sample clusters that share bait
    group_union <- '__'
    message( 'Generating agg bait prey dataframe...' )
    agg_samples_bait_prey.clusters <- ddply( bait_prey.clusters, c( 'proteins.cluster', 'protein_ac' ), function( baits_cluster_data ) {
            samples.clusters.serials <- sort( unique( baits_cluster_data$samples.cluster ) )
            samples <- sort( unique( baits_cluster_data$sample ) )
            data.frame( samples.clusters.serials = ifelse( length( samples.clusters.serials ) > 0,
                                                           paste( samples.clusters.serials, collapse = group_union ),
                                                           NA ),
                        samples = paste( samples, collapse = group_union ),
                        stringsAsFactors = FALSE )
        })
    gen_proteins_group_id <- function( proteins_cluster_serial ) paste( 'prot_grp', proteins_cluster_serial, sep = '_' )
    gen_samples_group_id <- function( samples_clusters_serials ) paste( 'sample_grp', paste( samples_clusters_serials, collapse = group_union ), sep = '_' )
    #print(agg_samples_bait_prey.clusters)
    message( 'Generating nodes...' )
    nodes.df <- ddply( agg_samples_bait_prey.clusters, c( 'proteins.cluster' ), function( proteins_cluster_data ) {
            prot_clu_serial <- unique( proteins_cluster_data$proteins.cluster )
            prot_group_id <- gen_proteins_group_id( prot_clu_serial )
            proteins_inside <- unique( proteins_cluster_data$protein_ac )
            nproteins <- length( proteins_inside )
            baits_inside <- unique( subset( proteins_cluster_data, !is.na(samples.clusters.serials) )$protein_ac )
            preys_inside <- setdiff( proteins_inside, baits_inside )
            npreys <- length( preys_inside )
            # define bait clusters inside protein cluster
            # this bait clusters are, generally, non-overlapping partitions of sample clusters:
            # if 2 sample clusters share baits in common --there would be 3 baits clusters -- for 1st only, 2nd only and overlapping baits
            res <- ddply( subset( proteins_cluster_data, !is.na( samples.clusters.serials ) ), c( 'samples.clusters.serials' ), function( samples_clusters_data ) {
                samples_clusters_serials <- unique(unlist(strsplit(samples_clusters_data$samples.clusters.serials,group_union)))
                samples <- unique(unlist(strsplit(samples_clusters_data$samples,group_union)))
                scs_regex <- paste( '(\\b', samples_clusters_serials, '\\b)', sep='', collapse='|' )
                #print(scs_regex)
                samples_group_id <- gen_samples_group_id( samples_clusters_serials )
                # check if proteins of proteins cluster are baits of samples cluster
                # find all baits of given samples clusters part
                all_baits <- unique( subset( agg_samples_bait_prey.clusters, grepl( scs_regex, samples.clusters.serials ) )$protein_ac )
                nall_baits = length( all_baits )
                # baits inside this proteins cluster
                baits <- unique( samples_clusters_data$protein_ac )
                nbaits <- length( baits )
                if ( nbaits == 0 ) stop( 'No baits inside proteins cluster, linked to this bait' )
                #print( all_baits )
                #print( baits )

                # calculate type of graph node
                is_group <- nbaits > 1
                is_embedded <- nproteins > nbaits # proteins cluster larger than baits cluster [part]
                is_partial <- nbaits < nall_baits  # proteins cluster does not fully contain baits cluster

                parent_cluster_id <- ifelse( is_embedded, prot_group_id, NA )
                if ( is_group ) {
                    cluster_node_id = paste( prot_group_id, samples_group_id, sep = ':' ) #subgroup
                    # baits cluster -- create nodes for bait and group node for cluster
                    return( rbind(
                        # create nodes for baits
                        data.frame( node_type = 'bait', node_id = baits,
                            parent_node_id = cluster_node_id,
                            is_group = FALSE, is_embedded = TRUE, is_partial = FALSE,
                            proteins_cluster = NA,
                            samples_clusters = NA, samples = NA,
                            stringsAsFactors = FALSE ),
                        # group node for cluster
                        data.frame( node_type = 'baits_cluster', node_id = cluster_node_id,
                            parent_node_id = ifelse( is_embedded, parent_cluster_id, NA ),
                            is_group = TRUE, is_embedded = is_embedded, is_partial = is_partial,
                            proteins_cluster = ifelse( is_embedded, NA, prot_clu_serial ),
                            samples_clusters = paste( samples_clusters_serials, sep='', collapse=group_union ),
                            samples = paste( samples, sep='', collapse=group_union ),
                            stringsAsFactors = FALSE )
                    ) )
                }
                else {
                    # single bait -- do not do grouping
                return ( data.frame( node_type = 'bait',
                        node_id = baits, parent_node_id = parent_cluster_id,
                        is_group = is_group, is_embedded = is_embedded, is_partial = is_partial,
                        proteins_cluster = ifelse( is_embedded, NA, prot_clu_serial ),
                        samples_clusters = paste( samples_clusters_serials, sep='', collapse=group_union ),
                        samples = paste( samples, sep='', collapse=group_union ),
                        stringsAsFactors = FALSE ) )
                }
            } )
            res$samples.clusters.serials <- NULL
            # add pure preys
            if ( npreys > 0 ) {
                # add preys to the cluster
                res <- rbind( res,
                    # prey nodes
                    data.frame( node_type = 'prey',
                        node_id = preys_inside, parent_node_id = ifelse( nproteins > 1, prot_group_id, NA ),
                        is_group = FALSE, is_embedded = nproteins > 1, is_partial = FALSE,
                        proteins_cluster = ifelse( nproteins > 1, NA, prot_clu_serial ),
                        samples_clusters = NA, samples = NA,
                        stringsAsFactors = FALSE )
                )
            }
            # add group node, if it's referenced by created nodes and is not yet added
            if ( any( res$parent_node_id == prot_group_id, na.rm = TRUE ) && sum( res$node_id == prot_group_id, na.rm = TRUE ) == 0 ) {
                res <- rbind( res,
                    data.frame( node_type = 'proteins_cluster', node_id = prot_group_id,
                        parent_node_id = NA,
                        is_group = TRUE, is_embedded = FALSE, is_partial = FALSE,
                        proteins_cluster = prot_clu_serial,
                        samples_clusters = NA, samples = NA,
                        stringsAsFactors = FALSE )
                )
            }
            return ( res )
        } )
    nodes.df$proteins.cluster <- NULL
    #print( colnames(nodes.df))
    # add nodes for samples clusters, that are partitioned across protein clusters
    # or overlap with other samples clusters
    split_samples_clusters_serials <- unique( unlist( strsplit( paste( subset( nodes.df, is_partial | grepl( group_union, samples_clusters, fixed = TRUE ) )$samples_clusters,
                                                       collapse = group_union ), group_union, fixed = TRUE ) ) )
    split_samples_clusters_serials <- split_samples_clusters_serials[ split_samples_clusters_serials != 'NA' ]
    #print( split_samples_clusters_serials )
    samples_clusters_nodes.df <- do.call( 'rbind', lapply( split_samples_clusters_serials,
        function ( samples_cluster_serial ) {
            samples_group_id <- gen_samples_group_id( samples_cluster_serial )
            data.frame( node_type = 'samples_cluster', node_id = samples_group_id,
                        parent_node_id = NA,
                        is_group = FALSE, is_embedded = FALSE, is_partial = FALSE,
                        proteins_cluster = NA,
                        samples_clusters = samples_cluster_serial,
                        samples = paste( subset( bimap.props$samples.clusters, samples.cluster == samples_cluster_serial )$sample,
                                         collapse = group_union ),
                        stringsAsFactors = FALSE )
        } ) )
    #print( samples_clusters_nodes.df )
    if ( !is.null( samples_clusters_nodes.df ) && nrow( samples_clusters_nodes.df ) > 0 ) {
        nodes.df <- rbind( nodes.df, samples_clusters_nodes.df )
    }

    # add information
    # add protein information
    message( "Adding protein information...")
    protein_node_mask <- nodes.df$node_type %in% c( 'bait', 'prey' ) &
                         nodes.df$node_id %in% rownames( bimap.data$proteins )
    nodes.df$short_id <- NA
    if ( !is.null( protein_label_col ) ) {
        nodes.df[ protein_node_mask, 'short_id' ] <-
            as.character( bimap.data$proteins[ nodes.df[ protein_node_mask, 'node_id' ], protein_label_col ] )
    }
    # export protein columns
    exported_cols <- names( protein_export_cols )
    if ( is.null( exported_cols ) ) exported_cols <- protein_export_cols
    nodes.df[,exported_cols] <- NA
    for ( col_ix in 1:length( protein_export_cols ) ) {
        external_col <- exported_cols[[ col_ix ]]
        internal_col <- protein_export_cols[[ col_ix ]]
        nodes.df[ protein_node_mask, external_col ] <-
            bimap.data$proteins[ nodes.df[ protein_node_mask, 'node_id' ], internal_col ]
    }
    # add proteins cluster information
    proteins_clu_node_mask <- nodes.df$node_type == 'proteins_cluster'
    if ( !is.null( proteins_cluster_label_col ) ) {
        nodes.df[ proteins_clu_node_mask, 'short_id' ] <-
                as.character( bimap.props$proteins.clusters.info[ nodes.df[ proteins_clu_node_mask, 'proteins_cluster' ], proteins_cluster_label_col ] )
        # print( bimap.props$proteins.clusters.info[ nodes.df[ proteins_clu_node_mask, 'proteins_cluster' ], ] )
    } else {
        nodes.df[ proteins_clu_node_mask, 'short_id' ] <- nodes.df[ proteins_clu_node_mask, 'proteins_cluster' ]
    }
    for ( clu_col in colnames( bimap.props$proteins.clusters.info ) ) {
        nodes.df[ proteins_clu_node_mask, clu_col ] <-
                bimap.props$proteins.clusters.info[ nodes.df[ proteins_clu_node_mask, 'proteins_cluster' ], clu_col ]
    }
    nodes.df$experiment_description <- NA
    nodes.df[ proteins_clu_node_mask, 'experiment_description' ] <- paste( 'Proteins cluster #',
              nodes.df[ proteins_clu_node_mask, 'proteins_cluster' ], sep = '' )
    # add samples cluster information
    if ( is.character( sample_label_col ) ) {
        samples_clu_node_mask <- nodes.df$node_type %in% c( 'baits_cluster', 'samples_cluster' ) & !is.na( nodes.df$samples_clusters )
        nodes.df[ samples_clu_node_mask, 'experiment_description' ] <- sapply( nodes.df[ samples_clu_node_mask, 'samples_clusters' ],
            function( samples_clusters_concat )
            {
                samples_clusters <- unlist( strsplit( samples_clusters_concat, group_union, fixed = TRUE ) )
                samples_stats <- table( subset( bimap.props$samples.clusters, samples.cluster %in% samples_clusters )$sample )
                # samples in the intersection of all current clusters
                intersect_samples <- names( samples_stats )[samples_stats == length( samples_clusters )]
                intersect_samples <- sort( unique( bimap.data$samples[ intersect_samples, sample_label_col] ) )
                return ( paste( intersect_samples, collapse = '\n' ) )
            } )
        nodes.df[ samples_clu_node_mask, 'short_id' ] <- nodes.df[ samples_clu_node_mask, 'experiment_description' ]
        samples_bait_node_mask <- nodes.df$node_type == 'bait'
        nodes.df[ samples_bait_node_mask, 'experiment_description' ] <- sapply( nodes.df[ samples_bait_node_mask, 'node_id' ], function( cur_bait_ac ) {
                samples <- subset( bimap.data$samples, bait_ac == cur_bait_ac )
                sample_labels <- sort( unique( samples[,sample_label_col] ) )
                return ( paste( sample_labels, collapse = '\n' ) )
            } )
        #print( nodes.df[ samples_clu_node_mask,] )
    }

    # add edges for each BIMAP block
    message( 'Generating edges...' )
    signal.max <- max( bimap.props$signals.mean, na.rm = TRUE )
    signal.min <- min( bimap.props$signals.mean, na.rm = TRUE )

    edges.df <- ddply( bimap.props$blocks, c( 'proteins.cluster', 'samples.cluster' ), function( blk ) {
            # TODO: if directed, there could be a link from inner to outer node
            source_node <- subset( nodes.df, blk$samples.cluster == samples_clusters
                                            #& ( proteins_cluster != blk$proteins.cluster | is.na( proteins_cluster ) )
                                            & !is_partial )
            target_node <- subset( nodes.df, proteins_cluster == blk$proteins.cluster & !is_embedded )
            if ( nrow( source_node ) == 1 && nrow( target_node ) == 1 ) {
                ocId <- as.character( blk$proteins.cluster )
                scId <- as.character( blk$samples.cluster )
                #print( paste( '[', ocId, ',', blk_samples, ']' ) )
                #print( bimap.props$signals.mean[ ocId, blk_samples ] )
                #print( bimap.props$signals.sd[ ocId, blk_samples ] )
                if ( ( !any( source_node$proteins_cluster == target_node$proteins_cluster, na.rm = TRUE )
                     || any( source_node$is_partial & !target_node$is_partial, na.rm = TRUE ) )
                  &&  !any( grepl( paste( group_union, blk$samples.cluster, group_union, sep=''),
                                   paste( group_union, target_node$samples_clusters, group_union, sep='' ),
                                          fixed = TRUE ) )
                ){
                    res <- data.frame( source_id = source_node$node_id, target_id = target_node$node_id,
                                type = 'bimap_block',
                                abundance = bimap.props$signals.mean[ ocId, scId ],
                                abundance_sd = NA,
                                intensity.scaled = (bimap.props$signals.mean[ ocId, scId ] - signal.min) / (signal.max - signal.min),
                                stringsAsFactors = FALSE )
                    res$hex.intensity = format( as.hexmode( min.intensity + as.integer( (255-min.intensity) * res$intensity.scaled^hex.intensity.power ) ),
                                                width=2, upper.case=TRUE)
                    if ( !is.null( bimap.props$signals.sd ) ) {
                        res$abundance_sd = bimap.props$signals.sd[ ocId, scId ]
                    }
                    return( res )
                } else {
                    # don't create edge from the node to itself
                    data.frame( source_id = list(), target_id = list(), type = list(),
                                stringsAsFactors = FALSE )
                }
            } else {
                if ( nrow( source_node ) > 1 ) {
                    warning( "Samples cluster #", blk$samples.cluster, ": ", nrow( source_node ), " source nodes" )
                } else if ( nrow( source_node ) == 0 ) {
                    warning( "Samples cluster #", blk$samples.cluster, ": no source nodes. Is bait protein present in proteins clusters?" )
                }
                if ( nrow( target_node ) != 1 ) {
                    warning( "Proteins cluster #", blk$proteins.cluster, ": ", nrow( target_node ), " target nodes" )
                }
                # no edges for embedded nodes
                data.frame( source_id = list(), target_id = list(), type = list(),
                            stringsAsFactors = FALSE )
            }
        } )
    # remove tentative columns
    edges.df$proteins.cluster <- NULL
    edges.df$samples.cluster <- NULL
    #print(edges.df)

    # add edges for split-sample clusters
    if ( !is.null( samples_clusters_nodes.df) && nrow( samples_clusters_nodes.df ) > 0 ) {
        message( 'Generating split-sample edges...' )
        split_sample_edges.df <- ddply( samples_clusters_nodes.df, c( 'node_id '), function( node_data ) {
            samples_cluster_serial <- node_data$samples_clusters
            split_nodes <- subset( nodes.df, grepl( paste(group_union, samples_cluster_serial, group_union, sep = '' ),
                                                    paste(group_union,samples_clusters,group_union,sep='') ) & node_id != node_data$node_id )
            data.frame( source_id = gen_samples_group_id( samples_cluster_serial ),
                        target_id = split_nodes$node_id, type = 'samples_cluster_composition',
                        abundance = NA, abundance_sd = NA,
                        intensity.scaled = NA,
                        hex.intensity = NA,
                        stringsAsFactors = FALSE )
            } )
        split_sample_edges.df$node_id <- NULL
        edges.df <- rbind( edges.df, split_sample_edges.df )
    }
    message( 'BIMAP GraphML data.frame generation done' )
    return ( list( nodes = nodes.df, edges = edges.df ) )
}

BIMAP.extract_graphML.dataframe <- function( bimap.walk, bimapId,
    bimap.data,
    onblock.threshold = 0.6, ...
){
    bimap.props <- BIMAP.extract_clustering( bimap.walk, bimapId, TRUE, 
        onblock.threshold = onblock.threshold )
    BIMAP.graphML.dataframe( bimap.props, bimap.data, ... )
}

BIMAP.graphML <- function( bimap.props, bimap.data, ...
){
    dfs <- BIMAP.graphML.dataframe( bimap.props, bimap.data, ... )
    return ( generateGraphML( dfs$nodes, dfs$edges,
            node_col = 'node_id', parent_col = 'parent_node_id',
            source_col = 'source_id', target_col = 'target_id',
            directed = FALSE ) )
}

BIMAP.extract_graphML <- function( bimap.walk, bimapId,
    bimap.data,
    onblock.threshold = 0.5,
    ...
){
    dfs <- BIMAP.extract_graphML.dataframe( bimap.walk, bimapId, bimap.data, ... )
    return ( generateGraphML( dfs$nodes, dfs$edges,
                    node_col = 'node_id', parent_col = 'parent_node_id',
                    source_col = 'source_id', target_col = 'target_id',
                    directed = FALSE ) )
}
