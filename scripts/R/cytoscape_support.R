# Starting Cytoscape with all node/edge attributes and Viz mapping automatically
# loaded.
# 
# Author: astukalov
###############################################################################

# description of the network, node/edge attributes and style for Cytoscape 
setClass( "CytoscapeViewDescription",
        representation(
                networkFile = "character",
                nodeAttributeFiles = "list",
                edgeAttributeFiles = "list",
                vizMapFile = "character"
        )
)

setGeneric( "invoke",
        function( .Object, ... ) standardGeneric( "invoke" )
)

setMethod( "invoke", signature( .Object = "CytoscapeViewDescription" ),
        function( .Object, cytoscape_path = cytoscape_path ) {
            options <- paste( "--network", .Object@networkFile )
            if ( length( .Object@nodeAttributeFiles ) > 0 ) {
                        options <- paste( options, paste(
                        lapply( .Object@nodeAttributeFiles, 
                                function( nodeAttrFile ) paste("--node-attrs", shQuote(nodeAttrFile) ) ), collapse = ' ' )
                )
            }
            if ( length( .Object@edgeAttributeFiles ) > 0 ) {
                options <- paste( options, paste( 
                        lapply( .Object@edgeAttributeFiles, 
                                function( edgeAttrFile ) paste("--edge-attrs", shQuote(edgeAttrFile) ) ), collapse = ' ' ) 
                )
            }
            if ( length( .Object@vizMapFile ) > 0 ) {
                options <- paste( options, "--vizmap", .Object@vizMapFile )
            }
            options <- paste( options, "--plugin", file.path( cytoscape_path, 'plugins' ) ) 
            cmd_line <- paste( file.path( cytoscape_path, "cytoscape.sh" ), options ) 
            print( cmd_line )
            system( cmd_line, wait = FALSE )
        }
)

# common class for node and edge attribute descriptors 
setClass( "CytoscapeAttributeDescriptor",
    representation( label = "character", 
        description = "character",
        valueClass = "character"
    )
)

setClass( "CytoscapeNodeAttributeDescriptor",
    contains = "CytoscapeAttributeDescriptor",
    representation( 
        nodeValueFunc = "function",
        contextFunc = "function"
    )
)

setMethod( 'initialize', "CytoscapeNodeAttributeDescriptor",
    function( .Object, label = "unknown", description = "<Unknown Attribute>", valueClass = '', 
              contextFunc = function( context ) return ( NULL ),
              nodeValueFunc = function( node, context, ... ) return ( NULL )
    ){
        .Object@label <- label
        .Object@description <- description
        .Object@valueClass <- valueClass
        .Object@contextFunc <- contextFunc
        .Object@nodeValueFunc <- nodeValueFunc
        return ( .Object )
    }
)

CytoscapeNodeAttributeDescriptor <- function( label, description, valueClass = '', 
    contextFunc = function( context ) return ( NULL ),
    nodeValueFunc = function( node, context, ... ) return ( NULL )
){
    return ( new( 'CytoscapeNodeAttributeDescriptor', 
            label=label, description=description, valueClass=valueClass, 
            contextFunc = contextFunc, nodeValueFunc = nodeValueFunc
        ) )
}

setClass( "CytoscapeEdgeAttributeDescriptor",
    contains = "CytoscapeAttributeDescriptor",
    representation( 
        edgeValueFunc = "function"
    )
)

setMethod( 'initialize', "CytoscapeEdgeAttributeDescriptor",
    function( .Object, label = "unknown", description = "<Unknown Attribute>", valueClass = '', 
        edgeValueFunc = function( outNode, inNode ) return ( NULL )
    )
    {
        .Object@label <- label
        .Object@description <- description
        .Object@valueClass <- valueClass
        .Object@edgeValueFunc <- edgeValueFunc
        return ( .Object )
    }
)

CytoscapeEdgeAttributeDescriptor <- function( label, description, valueClass = '', 
    edgeValueFunc = function( outNode, inNode ) return ( NULL )
){
    return ( new( 'CytoscapeEdgeAttributeDescriptor', 
                  label=label, description=description, valueClass=valueClass, 
                  edgeValueFunc = edgeValueFunc
        ) )
}

openCytoscapeAttributeFile <- function( attrDescr, filename_prefix, ext = 'noa', path = '' )
{
    title <- attrDescr@description
    if ( length( title ) == 0 ) {
        title <- attrDescr@label
    }
    attr_filename <- paste( filename_prefix, attrDescr@label, ext, sep='.' )
    if ( nchar( path ) > 0 ) attr_filename <- file.path( path, attr_filename )
    attr_file <- file( attr_filename, "w" )
    # start the file
    cat( title, file = attr_file )
    if ( nchar( attrDescr@valueClass ) > 0 ) {
        cat( " (class=", attrDescr@valueClass, ")", file = attr_file, sep='' )
    }
    cat( '\n', file = attr_file )
    return ( list( file = attr_file, filename = attr_filename ) )
}

setGeneric( "prepareHyperedgeNodes",
    function( .Object, hyperedge, ... ) standardGeneric( "prepareHyperedgeNodes" )
)

setMethod( "prepareHyperedgeNodes", signature( .Object = "CytoscapeNodeAttributeDescriptor", hyperedge = "Hyperedge" ),
    function( .Object, hyperedge, ... ) {
        return ( NULL )
    }
)

setGeneric( "nodeValue",
        function( .Object, node, isHead = true, hyperedgeInfo = NA, ... ) standardGeneric( "nodeValue" )
)

setClass( "CytoscapeFilesWriter",
    representation = representation(
        filename_prefix = "character",
        description = "CytoscapeViewDescription"
    )
)


setGeneric( "cytoDescribe", function( .Object ) standardGeneric( "cytoDescribe" ) )

#' 
#' @param .Object 
#' @returnType CytoscapeViewDescriptor
#' @return view descriptor object, to be passed into invoke() method
#' @author astukalov
#' @export
setMethod( 'cytoDescribe', signature( .Object = "CytoscapeFilesWriter" ),
    function ( .Object ) {
        return ( .Object@description )
    }
)

setGeneric( "cytoWrite", function( .Object, writable, ... ) standardGeneric( "cytoWrite" ) )

cytoWriteProteinAcc <- function( cytoWriter )
{
    cytoWrite( CytoscapeNodeAttributeDescriptor( "label", "Label", 
                    nodeValue = function( node, isHead, hyperedgeInfo ) return( label( node ) ) ) )
    cytoWrite( CytoscapeNodeAttributeDescriptor( "uniprot", "UniProt Accession", 
        nodeValueFunc = function( node, isHead, hyperedgeInfo ) return( node@uniprot_ac ) ) )
}

CytoscapeNodeDataSetAttribute <- function( node_data_set, node_id_field, info_field, description )
{
    filter_expr <- parse(text=paste('subset(node_data_set,',node_id_field,'==node@id)[1,\'', info_field,'\']', sep='' ))
    CytoscapeNodeAttributeDescriptor( info_field, description, 
        nodeValue = function( node, isHead, hyperedgeInfo ) {
            return( eval(filter_expr) ) } )
}

CytoscapeEdgeDataSetAttribute <- function( edge_data_set, out_node_id_field, in_node_id_field, info_field, description )
{
    filter_expr <- parse(text=paste('as.character(subset(edge_data_set,',out_node_id_field,'==outNode & ',in_node_id_field,'==inNode)$', info_field,')', sep='' ))
    CytoscapeEdgeAttributeDescriptor( info_field, description, 
        edgeValue = function( outNode, inNode ) {
            return( eval(filter_expr) ) } )
}
