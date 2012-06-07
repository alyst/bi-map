require( xlsx )
require( gplots )

source( file.path( bimap_scripts_path, "BIMAP_plot.R" ) )

#' Generate Excel sheet for biclustering
#' @param bimap.props BI-MAP biclustering
#' @param bimap.data input data for BI-MAP
#' @param blocks.pal palette to use for on-blocks
#' @param cells.pal palette to use for non-empty cells in off-blocks
#' @param grid.col grid color
#' @param bait.border.col grid color for detected baits
#' @param col.width width of columns
#' @param workbook XLSX workbook object, if NULL the new workbook is created 
#' @param sheet.name name of the BI-MAP results sheet
#' @param ... parameters for BIMAP.prepare_plot() 
#' @returnType 
#' @return XLSX object
#' @author Alexey Stukalov
#' @see BIMAP.plot_prepare(), BIMAP.plot()
BIMAP.create_xlsx <- function( 
    bimap.props, bimap.data,
    blocks.pal = colorRamp( c("blue","cyan","yellow") ),
    cells.pal = colorRamp( c( "pink","yellow") ),
    grid.col = 'black',
    bait.border.col = 'red',
    col.width = 5,
    workbook = NULL,
    sheet.name = 'BI-MAP results',
    ...
){
    message( 'Preparing BI-MAP output...' )
    args = list(...)
    bimap.plot_internal <- do.call( 'BIMAP.plot_prepare', c( list( bimap.props, bimap.data ), args ) )

    putBorder <- function( block, rowIndex, colIndex, border ) {
        CB.setBorder( block,
                      rowIndex = rep.int( rowIndex, length(colIndex) ),
                      colIndex = rep( colIndex, each = length(rowIndex) ),
                      border )
    }

    message( 'Generating BI-MAP XLSX...' )
    bimap.workbook <- if ( is.null( workbook ) ) createWorkbook( type = 'xlsx' ) else workbook
    wsh <- createSheet( bimap.workbook, sheet.name )
    with( bimap.plot_internal, {
        #print( head( proteins ) )
        #print( head( samples ) )
        prot_cols <- c( 'protein_ac', 'name', 'description' )
        prot_cols <- c( intersect( prot_cols, colnames( proteins ) ),
                        setdiff( colnames( proteins ), c( prot_cols, 'proteins.cluster', 'short_label' ) ) )
        message( 'Protein info to output: ', paste( prot_cols, ' ' ) )
        samp_cols <- c( 'bait_short_label', 'bait_ac', 'sample', 'msrun' )
        samp_cols <- c( intersect( samp_cols, colnames( samples ) ),
                        setdiff( colnames( samples ), c( samp_cols, 'samples.cluster', 'short_label', 'col_id' ) ) )
        message( 'Samples info to output: ', paste( samp_cols, ' ' ) )
        row_offset <- length( samp_cols )
        col_offset <- length( prot_cols )

        createFreezePane( wsh, 1 + row_offset, 1 + col_offset )

        proteins$order <- 1:nrow(proteins)
        samples$order <- 1:nrow(samples)
        proteins <- ddply( proteins, c( 'proteins.cluster' ), function( proteins_clu ) {
            proteins_clu$clu.order <- 1:nrow( proteins_clu )
            return ( proteins_clu )
        } )
        proteins <- proteins[ order( proteins$order ), ]
        samples <- ddply( samples, c( 'samples.cluster' ), function( samples_clu ) {
            samples_clu$clu.order <- 1:nrow( samples_clu )
            return ( samples_clu )
        } )
        samples <- samples[ order( samples$order ), ]

        header_style <- CellStyle( bimap.workbook ) + Fill( 'lightblue' ) +
                        Font( bimap.workbook, isItalic = TRUE, isBold = TRUE )
        info_style <- CellStyle( bimap.workbook ) +
                       Font( bimap.workbook, isItalic = TRUE )

        message( 'Samples data writing...' )
        sampleInfoBlock <- addDataFrame( samples[,samp_cols], wsh,
                      startRow = 1, startCol = col_offset,
                      col.names = TRUE, row.names = FALSE,
                      byrow = TRUE, colnamesStyle = header_style,
                      colStyle = info_style )
        # cluster borders
        putBorder( sampleInfoBlock,
                   colIndex = 1 + subset( samples, clu.order == 1 & order > 1 )$order,
                   rowIndex = 1:length(samp_cols),
                   Border( color = grid.col, position = c( 'LEFT' ), pen = c( 'BORDER_THIN' ) ) )
        # outer frame
        putBorder( sampleInfoBlock,
                   Border( color = grid.col, position = c( 'LEFT' ), pen = c( 'BORDER_THICK' ) ),
                   colIndex = 2, rowIndex =1:length(samp_cols) )
        putBorder( sampleInfoBlock,
                   Border( color = grid.col, position = c( 'RIGHT' ), pen = c( 'BORDER_THICK' ) ),
                   colIndex = 1 + nrow(samples), rowIndex = 1:length(samp_cols) )
        putBorder( sampleInfoBlock,
                   Border( color = grid.col, position = c( 'TOP' ), pen = c( 'BORDER_THICK' ) ),
                   colIndex = 1 + 1:nrow(samples), rowIndex = 1 )
        putBorder( sampleInfoBlock,
                   Border( color = grid.col, position = c( 'BOTTOM' ), pen = c( 'BORDER_THICK' ) ),
                   colIndex = 1 + 1:nrow(samples), rowIndex = length(samp_cols) )
        message( 'Merging sample info columns' )
        for ( icol in 1:length(samp_cols) ) {
            ddply( samples, c( "samples.cluster", samp_cols[[icol]] ), function( merged_cells ) {
                addMergedRegion( wsh, icol, icol,
                                 col_offset + min( merged_cells$order ),
                                 col_offset + max( merged_cells$order ) )
            } )
        }
        setColumnWidth( wsh, col_offset + 1:nrow(samples), col.width )

        message( 'Proteins data writing...' )
        proteinInfoBlock <- addDataFrame( proteins[,prot_cols], wsh,
                      startRow = row_offset, startCol = 1,
                      col.names = TRUE, row.names = FALSE,
                      byrow = FALSE,
                      colStyle = info_style, colnamesStyle = header_style )
        # cluster borders
        putBorder( proteinInfoBlock,
                   rowIndex = 1 + subset( proteins, clu.order == 1 & order > 1 )$order,
                   colIndex = 1 + 1:length(samp_cols),
                   Border( color = grid.col, position = c( 'TOP' ), pen = c( 'BORDER_THIN' ) ) )
        # outer frame
        putBorder( proteinInfoBlock,
                   Border( color = grid.col, position = c( 'TOP' ), pen = c( 'BORDER_THICK' ) ),
                   rowIndex = 2, colIndex = 1:length(prot_cols) )
        putBorder( proteinInfoBlock,
                   Border( color = grid.col, position = c( 'BOTTOM' ), pen = c( 'BORDER_THICK' ) ),
                   rowIndex = 1 + nrow(proteins), colIndex = 1:length(prot_cols) )
        putBorder( proteinInfoBlock,
                   Border( color = grid.col, position = c( 'LEFT' ), pen = c( 'BORDER_THICK' ) ),
                   rowIndex = 1 + 1:nrow(proteins), colIndex = 1 )
        putBorder( proteinInfoBlock,
                   Border( color = grid.col, position = c( 'RIGHT' ), pen = c( 'BORDER_THICK' ) ),
                   rowIndex = 1 + 1:nrow(proteins), colIndex = length(prot_cols) )

        if ( !is.null( cells.on.matrix ) ) {
            message( 'MS-data writing...' )
            cells.matrix <- cells.on.matrix
            cells.matrix[ is.na( cells.on.matrix) ] <- cells.off.matrix[ is.na( cells.on.matrix ) ]
        } else {
            cells.matrix <- matrix( nrow = nrow( proteins ), ncol = ncol( samples ) )
        }
        biclustersBlock <- addDataFrame( cells.matrix, wsh,
                          startRow = 1 + row_offset,
                          startCol = 1 + col_offset,
                          col.names = FALSE, row.names = FALSE )
        message( 'Biclusters grid...' )
        putBorder( biclustersBlock,
                   rowIndex = subset( proteins, clu.order == 1 & order > 1 )$order,
                   colIndex = 1:ncol( cells.matrix ),
                   Border( color = grid.col, position = c( 'TOP' ), pen = c( 'BORDER_THIN' ) ) )
        putBorder( biclustersBlock,
                   rowIndex = 1:nrow( cells.matrix ),
                   colIndex = subset( samples, clu.order == 1 & order > 1 )$order,
                   Border( color = grid.col, position = c( 'LEFT' ), pen = c( 'BORDER_THIN' ) ) )

        message( 'On-Blocks color assignment...' )
        min.signal <- min( blocks$signal.mean )
        max.signal <- max( blocks$signal.mean )
        blocks <- within( blocks, {
            q <- ( signal.mean - min.signal ) / ( max.signal - min.signal )
            fill.color <- sapply( q, function( k ) do.call( rgb, as.list( blocks.pal( k )[1,]/255 ) ) )
        } )
        
        # TODO: clear border for cells below or to the right of bait cells to preserve the color

        message( 'On-Blocks styling...' ) 
        for ( i in 1:nrow(blocks) ) {
            block <- blocks[i,]
            block_rows <- subset( proteins, proteins.cluster == block$proteins.cluster )$order
            block_cols <- subset( samples, samples.cluster == block$samples.cluster )$order
            if ( !is.na( block$fill.color ) ) {
                CB.setFill( biclustersBlock,
                            Fill( foregroundColor = blocks[ i, 'fill.color' ] ),
                            rowIndex = rep.int( block_rows, length(block_cols) ),
                            colIndex = rep( block_cols, each=length(block_rows)) )
            }
            createCellComment( CB.getCell( biclustersBlock, block_rows[1], block_cols[1] ),
                               string = paste( format( block$signal.mean, digits=3 ), "\u0b1",
                                         format( block$signal.sd, digits=3 ), sep="" ) )
        }

        if ( any( cells.off.matrix ) ) {
            message( 'Off-blocks cells styling...' )
            colnames(cells.off.matrix) <- NULL
            rownames(cells.off.matrix) <- NULL
            off.cells <- adply( cells.off.matrix, c(1,2), function(x) data.frame( data = x ) )
            off.cells <- subset( off.cells, !is.na( data ) )
            colnames(off.cells) <- c( 'order.protein', 'order.sample', 'data' )
            off.cells$order.protein <- as.integer( off.cells$order.protein )
            off.cells$order.sample <- as.integer( off.cells$order.sample )
            min.signal <- min( off.cells$data, na.rm = TRUE )
            max.signal <- max( off.cells$data, na.rm = TRUE ) + 1E-7
            off.cells$q <- ( off.cells$data - min.signal ) / ( max.signal - min.signal )
            off.cells$fill.color <- sapply( off.cells$q, function( k ) do.call( rgb, as.list( cells.pal( k )[1,]/255 ) ) )
            for ( i in 1:nrow(off.cells) ) {
                CB.setFill( biclustersBlock, Fill( foregroundColor = off.cells[ i, 'fill.color' ] ),
                            rowIndex = off.cells[i,'order.protein'],
                            colIndex = off.cells[i,'order.sample'] )
            }
        }

        message( 'Baits...' )
        bait_cells <- merge( subset( proteins, protein_ac %in% samples$bait_ac ),
                             samples, by.x = c( 'protein_ac' ), by.y = c( 'bait_ac' ),
                             suffixes = c( '.protein', '.sample' ), all = FALSE )
        bait_border = Border( color = bait.border.col, position = c( 'TOP', 'LEFT', 'RIGHT', 'BOTTOM' ), pen = 'BORDER_THICK' )
        CB.setBorder( biclustersBlock, bait_border,
                      rowIndex = bait_cells$order.protein,
                      colIndex = bait_cells$order.sample )

        message( "Outer frame..." )
        putBorder( biclustersBlock,
                   rowIndex = nrow( cells.matrix ), colIndex = 1:ncol( cells.matrix ),
                   Border( color = grid.col, position = c( 'BOTTOM' ), pen = c( 'BORDER_THICK' ) ) )
        putBorder( biclustersBlock,
                   rowIndex = 1:nrow( cells.matrix ), colIndex = ncol( cells.matrix ),
                   Border( color = grid.col, position = c( 'RIGHT' ), pen = c( 'BORDER_THICK' ) ) )
    })
    message( "BI-MAP XLSX generation done" )
    return ( bimap.workbook )
}
