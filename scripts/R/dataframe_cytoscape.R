source( paste( ppi_scripts_path, "cytoscape_support.R", sep = .Platform$file.sep ) )

# write Cytoscape's SIF file for given data.frame
writeDataframeSIF <- function( dataframe, col_bait_ac, col_prey_ac, col_edge_type, filename_prefix = 'ppi', directed_edges = TRUE )
{
    sif_filename <- paste( filename_prefix, ".sif", sep='' )
    sif_file <- file( sif_filename , "w")
    edges.df <- unique( data.frame( bait_ac = as.character( dataframe[ , col_bait_ac ] ),
            prey_ac = as.character( dataframe[ , col_prey_ac ] ),
            edge_type = dataframe[ , col_edge_type ] ) )
    if ( !directed_edges ) {
        # remove self-edges and edges to alphabetically preceding baits
        baits_ac <- unique( as.character( edges.df$bait_ac ) )
        edges.df <- subset( edges.df, !(prey_ac %in% baits_ac) | ( as.character(prey_ac) > as.character(bait_ac) ) )
    }
    apply( edges.df, 1, function( edge ) {
            bait_code <- edge[[ 'bait_ac' ]]
            prey_code <- edge[[ 'prey_ac' ]]
            edge_type <- edge[[ 'edge_type' ]]
            cat( bait_code, edge_type, prey_code, "\n", file = sif_file, sep = ' ' )
            return ()
        } )
    close( sif_file )
    return ( sif_filename )
}

setClass( "PPIDataFrameCytoscapeFilesWriter",
    contains = 'CytoscapeFilesWriter',
    representation = representation(
        dataframe = "data.frame",
        col_bait_ac = "character",
        col_prey_ac = "character",
        col_edge_type = "character",
        directed_edges = "logical"
    )
)

setMethod( 'initialize', "PPIDataFrameCytoscapeFilesWriter",
    function( .Object, table,
        col_bait_ac = 'bait_ac', col_prey_ac = 'prey_ac', col_edge_type = NULL,
        filename_prefix, vizmap_filename = NULL, directed_edges = TRUE
    ){
        # new empty description
        .Object@dataframe <- table
        .Object@dataframe[,col_bait_ac] <- as.character(.Object@dataframe[,col_bait_ac])
        .Object@dataframe[,col_prey_ac] <- as.character(.Object@dataframe[,col_prey_ac])
        if ( !is.null(col_edge_type) ) {
            .Object@dataframe[,col_edge_type] <- as.character(.Object@dataframe[,col_edge_type])
            .Object@col_edge_type <-col_edge_type
        }
        else {
            .Object@col_edge_type <- '__int_edge_type__'
            .Object@dataframe[,.Object@col_edge_type] <- rep( 'pp', nrow( .Object@dataframe ) )
        }
        .Object@filename_prefix = filename_prefix
        .Object@col_bait_ac = col_bait_ac
        .Object@col_prey_ac = col_prey_ac
        .Object@description = new( 'CytoscapeViewDescription', 
            networkFile = writeDataframeSIF( .Object@dataframe, .Object@col_bait_ac, .Object@col_prey_ac, .Object@col_edge_type,
                filename_prefix = .Object@filename_prefix, directed_edges = directed_edges ),
            nodeAttributeFiles = list(),
            edgeAttributeFiles = list(),
            vizMapFile = vizmap_filename
        )
        .Object@directed_edges = directed_edges
        return ( .Object )
    }
)

PPIDataFrameCytoscapeFilesWriter <- function( table,
    col_bait_ac = 'bait_ac', col_prey_ac = 'prey_ac', col_edge_type = NULL, 
	filename_prefix = NULL, vizmap_filename = NULL, 
	directed_edges = TRUE )
{
    if ( is.null( filename_prefix ) || nchar( filename_prefix ) == 0 ) {
        filename_prefix <- 'PPI'
    }
    return ( new( "PPIDataFrameCytoscapeFilesWriter", table,
                  col_bait_ac = col_bait_ac, col_prey_ac = col_prey_ac, col_edge_type = col_edge_type, 
                  filename_prefix = filename_prefix, vizmap_filename = vizmap_filename,
				  directed_edges = directed_edges ) )
}


RClassToCytoscapeValueClass <- function( className )
{
    switch ( className,
             numeric = 'Double',
             integer = 'Integer',
             character = 'String',
             logical = 'Boolean' )
}

setMethod( "cytoWrite", signature( .Object = "PPIDataFrameCytoscapeFilesWriter", writable = 'CytoscapeNodeAttributeDescriptor' ), 
    function( .Object, writable ) {
        # start node attribute file
        print( sprintf( "Node attribute '%s', Class=%s", writable@label, writable@valueClass ) )
        fileAndName <- openCytoscapeAttributeFile( writable, .Object@filename_prefix )

        writeNodeAttrib <- function( nodeAc, head_rows, tail_rows ) {
            nodeValue <- writable@nodeValueFunc( nodeAc, head_rows, tail_rows )
            if ( !is.null( nodeValue ) ) {
                cat( nodeAc, "=", nodeValue, file = fileAndName$file, sep='' )
                cat( '\n', file = fileAndName$file, sep='' )
            }
        }
        node_acs <- union( unique( .Object@dataframe[,.Object@col_bait_ac] ), 
            unique( .Object@dataframe[,.Object@col_prey_ac] ) )
        for ( node_ac in node_acs ) {
            # compose pairlist to pass to nodeAttribDescriptors (and calculate info per hyperedge)
            head_rows <- .Object@dataframe[ .Object@dataframe[,.Object@col_bait_ac] == node_ac, ]
            tail_rows <- .Object@dataframe[ .Object@dataframe[,.Object@col_prey_ac] == node_ac, ]
            writeNodeAttrib( node_ac, head_rows, tail_rows )
        }
        close( fileAndName$file )
        
        .Object@description@nodeAttributeFiles <- append( .Object@description@nodeAttributeFiles, fileAndName$filename )
        
        return ( .Object )
    }
)

setMethod( "cytoWrite", signature( .Object = "PPIDataFrameCytoscapeFilesWriter", writable = 'CytoscapeEdgeAttributeDescriptor' ), 
    function( .Object, writable ) {
        # update valueClass 
        if ( writable@valueClass == '' ) writable@valueClass = RClassToCytoscapeValueClass( class( .Object@dataframe[,writable@label] ) )
        print( sprintf( "Edge attribute '%s', Class=%s", writable@label, writable@valueClass ) )
        # start edge attribute files
        fileAndName <- openCytoscapeAttributeFile( writable, .Object@filename_prefix, ext = 'eda' )
        if ( .Object@directed_edges ) {
            edges.df <- .Object@dataframe
        }
        else {
            # remove self-edges and edges to alphabetically preceding baits
            baits_ac <- unique( as.character( .Object@dataframe[, .Object@col_bait_ac ] ) )
            edges.df <- .Object@dataframe[ !(.Object@col_prey_ac %in% baits_ac) | ( as.character( .Object@col_prey_ac ) > as.character( .Object@col_bait_ac ) ), ]
        }
        ddply( edges.df, c( .Object@col_bait_ac, .Object@col_prey_ac ), function( edge_rows ) {
                head_ac <- edge_rows[,.Object@col_bait_ac]
                tail_ac <- edge_rows[,.Object@col_prey_ac]
                edgeValue <- writable@edgeValueFunc( head_ac, tail_ac, edge_rows )
                if ( !is.null( edgeValue ) && !is.na( edgeValue ) ) {
                    cat( head_ac, " (pp) ", tail_ac, "=", edgeValue, file = fileAndName$file, sep='' )
                    cat( '\n', file = fileAndName$file, sep='' )
                }
                return (data.frame())
            }
        )
        close( fileAndName$file )
        
        .Object@description@edgeAttributeFiles <- append( .Object@description@edgeAttributeFiles,
            fileAndName$filename )
        return ( .Object )
    }
)

PPIDataFrameCytoscapeNodeDataSetAttribute <- function( info_field, description, dataset = NULL, key_column = NULL )
{
    CytoscapeNodeAttributeDescriptor( info_field, description,
        valueClass = RClassToCytoscapeValueClass( class(dataset[,info_field]) ), 
        nodeValue = function( node_ac, head_rows, tail_rows ) {
            if ( is.null( dataset ) ) {
                if ( nrow( tail_rows ) > 0 ) {
                    return( tail_rows[[1,info_field]] )
                }
                if ( nrow( head_rows ) > 0 ) {
                    return( head_rows[[1,info_field]] )
                } 
            } else {
                mask <- dataset[,key_column] == node_ac & !is.na(dataset[,info_field])
                if ( any(mask) ) {
                    return ( dataset[ mask, info_field][[1]] )
                } else {
                    return ( NULL );
                }
            }
        } )
}

PPIDataFrameCytoscapeEdgeDataSetAttribute <- function( info_field, description )
{
    return ( new( 'CytoscapeEdgeAttributeDescriptor',
              label=info_field, description=description,
              edgeValueFunc = function( outNode, inNode, rows ) return ( rows[,info_field] )
        ) )
}
