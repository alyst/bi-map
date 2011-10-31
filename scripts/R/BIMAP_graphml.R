# Exporting PSI bbb results to graphml
# 
# Author: astukalov
###############################################################################

require( plyr )

source( file.path( bimap_scripts_path, "graphml_support.R" ) )

bimap.graph.dataframe <- function( bimap.props,
    protein.info, sample.info, msrun.info,
    protein_label_col = 'short_label',
    protein_export_cols = c(),
    proteins_cluster_label_col = 'proteins.cluster',
    proteins_cluster_export_cols = c( 'avg.pairs.freq' ),
    min.intensity = 5, hex.intensity.power = 5
){
    # add to sample clusters dataframe
    rownames( sample.info ) <- sample.info$sample
    bimap.props$samples.clusters$bait_ac <- sample.info[ bimap.props$samples.clusters$sample, 'bait_ac' ]
    bimap.props$signal.max <- max( bimap.props$signals.mean, na.rm = TRUE )
    bimap.props$signal.min <- min( bimap.props$signals.mean, na.rm = TRUE )

    # merge samples and proteins clusters via bait-prey
    bait_prey.clusters <- merge( bimap.props$proteins.clusters, bimap.props$samples.clusters,
                                 by.x = 'protein_ac', by.y = 'bait_ac',
                                 all.x = TRUE, all.y = FALSE )
    # merge sample cluster ids for those sample clusters that share bait
    group_union <- '_'
    print( 'Generating agg bait prey dataframe' )
    agg_samples_bait_prey.clusters <- ddply( bait_prey.clusters, c( 'proteins.cluster', 'protein_ac' ), function( baits_cluster_data ) {
            samples.clusters.serials <- sort( unique( baits_cluster_data$samples.cluster ) )
            samples <- sort( unique( baits_cluster_data$sample ) )
            data.frame( samples.clusters.serials = ifelse( length( samples.clusters.serials ) > 0,
                                                           paste( samples.clusters.serials, collapse = group_union ),
                                                           NA ),
                        samples = paste( samples, collapse = group_union ),
                        stringsAsFactors = FALSE )
        })
    gen_proteins_group_id <- function( proteins_cluster_serial ) paste( 'proteinsgrp', proteins_cluster_serial, sep = '' )
    gen_samples_group_id <- function( samples_clusters_serials ) paste( 'samplesgrp', paste( samples_clusters_serials, collapse = group_union ), sep = '' )
    #print(agg_samples_bait_prey.clusters)
    print( 'Generating nodes' )
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
    protein_node_mask <- nodes.df$node_type %in% c( 'bait', 'prey' ) &
                         nodes.df$node_id %in% rownames( protein.info )
    nodes.df$short_id <- NA
    if ( !is.null( protein_label_col ) ) {
        nodes.df[ protein_node_mask, 'short_id' ] <-
            as.character( protein.info[ nodes.df[ protein_node_mask, 'node_id' ], protein_label_col ] )
    }
    # export protein columns
    for ( exported_col in names( protein_export_cols ) ) {
        internal_col = protein_export_cols[[ exported_col ]]
        nodes.df[,exported_col] <- NA
        nodes.df[ protein_node_mask, exported_col ] <-
            protein.info[ nodes.df[ protein_node_mask, 'node_id' ], internal_col ]
        
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
    samples_clu_node_mask <- nodes.df$node_type %in% c( 'baits_cluster', 'samples_cluster' ) & !is.na( nodes.df$samples_clusters )
    nodes.df[ samples_clu_node_mask, 'experiment_description' ] <- sapply( nodes.df[ samples_clu_node_mask, 'samples_clusters' ], function( samples_clusters ) {
            samples_clusters_list <- unlist( strsplit( samples_clusters, group_union, fixed = TRUE ) )
            samples_stats <- table( unlist( lapply( samples_clusters_list, function( sample_cluster ) {
                    subset( bimap.props$samples.clusters, samples.cluster == sample_cluster )$sample
                } ) ) )
            intersect_samples <- names( samples_stats )[samples_stats == length( samples_clusters_list )]
            return ( paste( 'Samples', paste( intersect_samples, collapse = ', ' ) ) )
        } )
    nodes.df[ samples_clu_node_mask, 'short_id' ] <- nodes.df[ samples_clu_node_mask, 'experiment_description' ]
    samples_bait_node_mask <- !is.na( nodes.df$node_type == 'bait' )
    nodes.df[ samples_bait_node_mask, 'experiment_description' ] <- sapply( nodes.df[ samples_bait_node_mask, 'node_id' ], function( cur_bait_ac ) {
            samples <- subset( sample.info, bait_ac == cur_bait_ac )$sample
            return ( paste( 'Samples', paste( samples, collapse = ', ' ) ) )
        } )
    #print( nodes.df )

    # add edges for each biclustering block
    print( 'Generating edges' )
    signal.max <- max( bimap.props$signals.mean, na.rm = TRUE )
    signal.min <- min( bimap.props$signals.mean, na.rm = TRUE )

    edges.df <- ddply( bimap.props$blocks, c( 'proteins.cluster', 'samples.cluster' ), function( cc ) {
            # TODO: if directed, there could be a link from inner to outer node
            source_node <- subset( nodes.df, cc$samples.cluster == samples_clusters
                                            & ( proteins_cluster != cc$proteins.cluster | is.na( proteins_cluster ) )
                                            & !is_partial )
            target_node <- subset( nodes.df, !grepl( paste(group_union,cc$samples.cluster,group_union,sep=''),
                                                     paste( group_union, samples_clusters, group_union, sep='') )
                    & proteins_cluster == cc$proteins.cluster & !is_embedded )
            if ( nrow( source_node ) == 1 && nrow( target_node ) == 1 ) {
                ocId <- as.character( cc$proteins.cluster )
                scId <- as.character( cc$samples.cluster )
                #print( paste( '[', ocId, ',', cc_samples, ']' ) )
                #print( bimap.props$signals.mean[ ocId, cc_samples ] )
                #print( bimap.props$signals.sd[ ocId, cc_samples ] )
                res <- data.frame( source_id = source_node$node_id, target_id = target_node$node_id,
                            type = 'block',
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
                if ( nrow( source_node ) > 1 ) {
                    warning( "Samples cluster #", cc$samples.cluster, ": ", nrow( source_node ), " source nodes" )
                } else if ( nrow( source_node ) == 0 ) {
                    warning( "Samples cluster #", cc$samples.cluster, ": no source nodes. Is bait protein present in proteins clusters?" )
                }
                if ( nrow( target_node ) != 1 ) {
                    warning( "Proteins cluster #", cc$proteins.cluster, ": ", nrow( target_node ), " target nodes" )
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
        print( 'Generating split-sample edges' )
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
    return ( list( nodes = nodes.df, edges = edges.df ) )
}

BIMAP.graph.dataframe <- function( bimap.walk, bimapId,
    protein.info, sample.info, msrun.info,
    onblock.threshold = 0.6, ...
){
    bimap.props <- BIMAP.extract_clustering( bimap.walk, bimapId, TRUE, 
        onblock.threshold = onblock.threshold )
    bimap.graph.dataframe( bimap.props, protein.info, sample.info, msrun.info, ... )
}

bimap.graphML <- function( bimap.props,
    protein.info, sample.info, msrun.info, ...
){
    dfs <- bimap.graph.dataframe( bimap.props, protein.info, sample.info, msrun.info, ... )
    return ( generateGraphML( dfs$nodes, dfs$edges,
            node_col = 'node_id', parent_col = 'parent_node_id',
            source_col = 'source_id', target_col = 'target_id',
            directed = FALSE ) )
}

BIMAP.graphML <- function( bimap.walk, bimapId, 
    protein.info, sample.info, msrun.info,
    onblock.threshold = 0.5,
    ...
){
    dfs <- BIMAP.graph.dataframe( bimap.walk, bimapId, protein.info, sample.info, msrun.info, ... )
    return ( generateGraphML( dfs$nodes, dfs$edges,
                    node_col = 'node_id', parent_col = 'parent_node_id',
                    source_col = 'source_id', target_col = 'target_id',
                    directed = FALSE ) )
}
