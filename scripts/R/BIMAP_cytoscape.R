# TODO: Add comment
# 
# Author: astukalov
###############################################################################

require(plyr)

source( paste( bimap_scripts_path, "cytoscape_support.R", sep = .Platform$file.sep ) )
source( paste( bimap_scripts_path, "BIMAP.R", sep = .Platform$file.sep ) )

setClass( "BIMAPCytoscapeFilesWriter",
    contains = c( 'CytoscapeFilesWriter' ),
    representation = representation(
        walk = "BIMAPWalk",
        clusteringId = "integer",
        objClusters = "data.frame",
        stateToBait = "character",
        simpleNodes = "character",
        stateClusters = "data.frame",
        blocks = "data.frame"
    )
)

setGeneric( "writeNNF", function( .Object, ... ) standardGeneric( "writeNNF" ) )

# write Cytoscape's NNF (nested networks format) file for given chessboard biclustering
# object clusters become nested networks
setMethod( 'writeNNF', "BIMAPCytoscapeFilesWriter", function( .Object ) {
    nnf_filename <- paste( .Object@filename_prefix, ".nnf", sep='' )
    nnf_file <- file( nnf_filename , "w")
    
    # output objects nodes
    ddply( .Object@objClusters, .( objects.cluster.serial ), function( cluster.rows ) {
        if ( nrow( cluster.rows ) == 1 ) {
            # simple node
            cat( 'outer', as.character( cluster.rows$object[[1]] ), sep='\t', file = nnf_file )
        } else {
            # nest objects of cluster within single node -
            #cat( cluster.rows$objects.cluster.serial[[1]], file = nnf_file )
            cat( paste( as.character( cluster.rows$objects.cluster.serial ), 
                    as.character( cluster.rows$object ), 
                    sep='\t', collapse='\n' ), 
                file = nnf_file )
        }
        cat( '\n', file = nnf_file )
        return ( 0 )
    } ) 
    
    # output states nodes and interactions
    ddply( .Object@blocks, .( objects.cluster.serial, state ), function( rows ) {
            row <- rows[1,]
            bait <- .Object@stateToBait[ as.character( row$state ) ]
            if ( bait %in% .Object@simpleNodes ) {
                bait_node <- bait
            }
            else {
                bait_node <- subset( .Object@objClusters, as.character( object ) == as.character( bait ) )$objects.cluster.serial[[1]]
            }
            prey_node <- as.character( row$objects.cluster.serial )
            if ( prey_node %in% names( .Object@simpleNodes ) ) {
                prey_node <- .Object@simpleNodes[ prey_node ]
            }
            cat( paste( 'outer', bait_node, 'pp', prey_node, sep='\t' ), file = nnf_file )
            cat( '\n', file = nnf_file )
            return ( 0 )
        })
    close( nnf_file )
    return ( nnf_filename )
} )

setMethod( 'initialize', "BIMAPCytoscapeFilesWriter",
    function( .Object, bimap.walk, clusteringId, samples, filename_prefix, vizmap_filename = NULL ) {
        # new empty description
        .Object@walk = bimap.walk
        .Object@clusteringId = as.integer( clusteringId )
        opId <- subset( bimap.walk@clusterings, clustering.serial == clusteringId )$objects.partition.serial
        spId <- subset( bimap.walk@clusterings, clustering.serial == clusteringId )$states.partition.serial
        ocIds <- subset( bimap.walk@objects.partitions, objects.partition.serial == opId )$objects.cluster.serial
        scIds <- subset( bimap.walk@states.partitions, states.partition.serial == spId )$states.cluster.serial
        .Object@objClusters <- subset( bimap.walk@objects.clusters, objects.cluster.serial %in% ocIds, 
            select = c( 'objects.cluster.serial', 'object' ) )
        simple_nodes_frame <- subset( ddply( .Object@objClusters, .(objects.cluster.serial), function( cluster.rows ) {
                    return ( data.frame( objects.cluster.serial = cluster.rows$objects.cluster.serial[[1]], 
                            object = cluster.rows$object[[1]], 
                            n = nrow( cluster.rows ) ) )
                } ), select = c("objects.cluster.serial", 'object'), n == 1 )
        .Object@simpleNodes <- as.character( simple_nodes_frame$object )
        names( .Object@simpleNodes ) <- as.character( simple_nodes_frame$objects.cluster.serial )
        .Object@stateClusters <- subset( bimap.walk@states.clusters, states.cluster.serial %in% scIds, 
            select = c( 'states.cluster.serial', 'state' ) )
        .Object@stateToBait = as.character( samples$bait_label )
        names( .Object@stateToBait ) = as.character( samples$label )
        #print( .Object@stateToBait )
        # not cross-clusters, actually, since state clusters are exploded into states
        .Object@blocks <- merge( subset( bimap.walk@blocks, clustering.serial == clusteringId,
            select = c( 'objects.cluster.serial', 'states.cluster.serial' ) ),
            subset( bimap.walk@states.clusters, states.cluster.serial %in% scIds ),
            by.x = c('states.cluster.serial'), by.y = c('states.cluster.serial') )[,c('objects.cluster.serial','state')]
        .Object@filename_prefix = filename_prefix
        .Object@description = new( 'CytoscapeViewDescription', 
            networkFile = writeNNF( .Object ),
            nodeAttributeFiles = list(),
            edgeAttributeFiles = list(),
            vizMapFile = vizmap_filename
        )
        return ( .Object )
    }
)

BIMAPCytoscapeFilesWriter <- function( bimap.walk, clusteringId, samples, 
        filename_prefix = NULL, vizmap_filename = NULL )
{
    if ( length( filename_prefix ) == 0 ) {
        filename_prefix <- paste('BIMAP', clusteringId, sep='-' )
    }
    return ( new( "BIMAPCytoscapeFilesWriter", bimap.walk, clusteringId, samples,
            filename_prefix = filename_prefix, vizmap_filename = vizmap_filename ) )
}

setMethod( "cytoWrite", signature( .Object = "BIMAPCytoscapeFilesWriter", 
                                   writable = 'CytoscapeNodeAttributeDescriptor' ), 
    function( .Object, writable, cclusOverride = NULL ) {
        # start node attribute file
        fileAndName <- openCytoscapeAttributeFile( writable, .Object@filename_prefix, ext = 'noa' )
        
        writeNodeAttrib <- function( object, isCluster, objectsClusterInfo, states ) {
            nodeValue <- writable@nodeValueFunc( object, objectsClusterInfo, isCluster = isCluster, states )
            if ( !is.null( nodeValue ) ) {
                cat( object, "=", nodeValue, '\n', file = fileAndName$file, sep='' )
            }
        }
        ddply( .Object@objClusters, .(objects.cluster.serial), function( rows ) {
                objectsClusterInfo <- writable@contextFunc( rows )
                lapply( 1:nrow(rows), function( ix ) {
                        obj <- as.character( rows$object[[ix]] )
                        states <- names( .Object@stateToBait )[ .Object@stateToBait %in% as.character( obj ) ]
                        writeNodeAttrib( obj, FALSE, objectsClusterInfo, states ) 
                } )
                if ( nrow( rows ) > 1 ) {
                    # cluster node
                    states <- names( .Object@stateToBait )[ .Object@stateToBait %in% as.character( rows$object ) ]
                    writeNodeAttrib( as.character( rows$objects.cluster.serial[[1]] ), TRUE, objectsClusterInfo, states ) 
                }
                return( 0 )
            } )
        close( fileAndName$file )
        .Object@description@nodeAttributeFiles <- append( .Object@description@nodeAttributeFiles, fileAndName$filename )
        return ( .Object )
    }
)

setMethod( "cytoWrite", signature( .Object = "BIMAPCytoscapeFilesWriter", 
                                   writable = 'CytoscapeEdgeAttributeDescriptor' ), 
    function( .Object, writable ) {
        cclus <- .Object@cclustering
        
        # start node attribute files
        fileAndName <- openCytoscapeAttributeFile( writable, .Object@filename_prefix, ext = 'eda' )
        
        apply( .Object@blocks, 1, function( row ) {
            bait <- .Object@stateToBait[ row$state ]
            if ( !( bait %in% .Object@simpleNodes ) ) {
                baitNode <- .Object@objClusters[ .Object@objClusters$object == bait, 'objects.cluster.serial' ]
            }
            else {
                baitNode <- bait
            }
            edgeValue <- writable@edgeValueFunc( bait, row$state, row$objects.cluster.serial )
            if ( !is.null( edgeValue ) ) {
                cat( baitNode, " (pp) ", 
                     row$objects.cluster.serial, "=", edgeValue, 
                     file = fileAndName$file, '\n', sep='' )
            }
        } )
        close( fileAndName$file )
        
        .Object@description@edgeAttributeFiles <- append( .Object@description@edgeAttributeFiles,
            fileAndName$filename )
        return ( .Object )
    }
)

BIMAPCytoWriteProteinInfo <- function( cytoWriter, protein_info )
{
    cytoWriter <- cytoWrite( cytoWriter, CytoscapeNodeAttributeDescriptor( "label", "Label", valueClass = "String",
            nodeValueFunc = function( node, contextInfo = objectsClusterInfo, isCluster, states ) {
                if ( isCluster ) {
                    return( NULL )
                } else {
                    return( smart_prot_id( node, protein_info ) )
                }
            } ) )
    cytoWriter <- cytoWrite( cytoWriter, CytoscapeNodeAttributeDescriptor( "uniprot", "UniProt Accession", valueClass = "String",
            nodeValueFunc = function( node, contextInfo = objectsClusterInfo, isCluster, states ) {
                if ( isCluster ) {
                    return( NULL )
                } else {
                    return( node )
                }
            } ) )
    cytoWriter <- cytoWrite( cytoWriter, CytoscapeNodeAttributeDescriptor( "description", "Description", valueClass = "String",
            nodeValueFunc = function( node, contextInfo = objectsClusterInfo, isCluster, states ) {
                if ( isCluster ) {
                    return( NULL )
                } else {
                    return( protein_info[ node, 'description' ] )
                }
            } ) )
    return ( cytoWriter )
}

BIMAPCytoscapeNodeDataSetAttribute <- function( node_data_set, node_id_field, info_field, description )
{
    filter_expr <- parse(text=paste('subset(node_data_set,',node_id_field,'==node@id)$', info_field, sep='' ))
    CytoscapeNodeAttributeDescriptor( info_field, description, 
        nodeValue = function( node, isHead, hyperedgeInfo ) {
            return( eval(filter_expr) ) } )
}

BIMAPCytoscapeEdgeDataSetAttribute <- function( edge_data_set, out_node_id_field, in_node_id_field, info_field, description )
{
    filter_expr <- parse(text=paste('as.character(subset(edge_data_set,',out_node_id_field,'==outNode & ',in_node_id_field,'==inNode)$', info_field,')', sep='' ))
    CytoscapeEdgeAttributeDescriptor( info_field, description, 
        edgeValue = function( outNode, inNode ) {
            return( eval(filter_expr) ) } )
}

BIMAPCytoWriteBasicAttributes <- function( cytoWriter )
{
    #cytoWriter <- cytoWrite( cytoWriter, CytoscapeNodeAttributeDescriptor( "detcnt", "Detections Count", valueClass="Integer",
    #        hyperedgeFunc = function( hyperedge ) return( length( hyperedge@intersectingEdges ) ), 
    #        nodeValue = function( node, isHead, hyperedgeInfo ) return( hyperedgeInfo ) ), hypergraphOverride = hgraphPtn )
    cytoWriter <- cytoWrite( cytoWriter, CytoscapeNodeAttributeDescriptor( "type", "Node Type", valueClass = "String",
            nodeValue = function( node, contextInfo = objectClusterInfo, isCluster, states ) {
                if ( isCluster ) {
                    return ( ifelse( length( states ) > 0, 'has_bait', 'preys' ) )
                }
                else {
                    return( ifelse( length( states ) > 0, "bait", "prey" ) )
                }
    } ) )
    cytoWriter <- cytoWrite( cytoWriter, CytoscapeNodeAttributeDescriptor( "samples", "Samples", valueClass = "String",
        nodeValue = function( node, contextInfo = objectClusterInfo, isCluster, states ) {
            return ( paste( states, collapse = ', ' ) )
        } ) )
    return ( cytoWriter )
}
