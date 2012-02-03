require( xlsx )
require( gplots )

source( file.path( bimap_scripts_path, "BIMAP_plot.R" ) )

BIMAP.create_xlsx <- function( 
        bimap.props, bimap.data,
        blocks.pal = colorRamp( c("blue","cyan","yellow") ),
        cells.pal = colorRamp( c( "pink","yellow") ),
        grid.col = 'black',
        bait.border.col = 'red',
        col.width = 5,
        ...
){
    message( 'Preparing for BI-MAP plotting...' )
    args = list(...)
    bimap.plot_internal <- do.call( 'BIMAP.plot_prepare', c( list( bimap.props, bimap.data ), args ) )

    message( 'Generating BI-MAP XLSX...' )
    bimap.workbook <- createWorkbook( type = 'xlsx' )
    wsh <- createSheet( bimap.workbook, title )
    with( bimap.plot_internal, {
        prot_cols <- intersect( c( 'protein_ac', 'short_label', 'description' ),
                                colnames(proteins) )
        samp_cols <- intersect( c( 'bait_ac', 'short_label', 'col_id' ),
                                colnames(samples) )
        row_offset <- length( samp_cols )
        col_offset <- length( prot_cols )

        createFreezePane( wsh, 1 + row_offset, 1 + col_offset )

        proteins$order <- 1:nrow(proteins)
        samples$order <- 1:nrow(samples)
        proteins <- ddply( proteins, c( 'proteins.cluster' ), function( proteins_clu ) {
            proteins_clu$clu.order <- 1:nrow( proteins_clu )
            return ( proteins_clu )
        } )
        samples <- ddply( samples, c( 'samples.cluster' ), function( samples_clu ) {
            samples_clu$clu.order <- 1:nrow( samples_clu )
            return ( samples_clu )
        } )

        message( 'Generating intermediate cells representation' )
        xl.cells <- merge( proteins, samples, by = c(), all = TRUE, suffixes = c( '.row', '.col' ) )

        message( 'Blocks color assignment...' )
        min.signal <- min( blocks$signal.mean )
        max.signal <- max( blocks$signal.mean )
        blocks <- within( blocks, {
            q <- ( signal.mean - min.signal ) / ( max.signal - min.signal )
            fill.color <- sapply( q, function( k ) do.call( rgb, as.list( blocks.pal( k )[1,]/255 ) ) )
        } )
        xl.cells <- merge( xl.cells, blocks, by = c( 'proteins.cluster', 'samples.cluster' ),
                           all.x = TRUE )
        xl.cells <- within( xl.cells, {
            is_bait <- protein_ac == bait_ac
            block_on <- !is.na( signal.mean )
            has_data <- FALSE
            off_data <- FALSE
            border.left <- clu.order.col == 1 | is_bait
            border.top <- clu.order.row == 1 | is_bait
            border.right <- is_bait
            border.bottom <- is_bait
            border.color <- ifelse( is_bait, bait.border.col,
                            ifelse( border.left | border.top | border.right| border.bottom, grid.col, NA ) )
        } )
        # TODO: clear border for cells below or to the right of bait cells to preserve the color
        # xl.cells.bait <- subset( xl.cells, is_bait )
        

        if ( !is.null( cells.on.matrix ) ) {
            message( 'MS-data writing...' )
            cells.matrix <- cells.on.matrix
            cells.matrix[ is.na( cells.on.matrix) ] <- cells.off.matrix[ is.na( cells.on.matrix ) ]
            xl.cells$data <- apply( xl.cells, 1, function(x) cells.matrix[ x[['protein_ac']], x[['col_id']] ] )
            xl.cells$has_data <- !is.na( xl.cells$data )
            xl.cells$off_data <- xl.cells$has_data & !xl.cells$block_on
            addDataFrame( cells.matrix, wsh,
                          startRow = 1 + row_offset,
                          startCol = 1 + col_offset,
                          col.names = FALSE, row.names = FALSE )
            if ( any( xl.cells$off_data ) ) {
                message( 'Off-blocks cells color assignment' )
                min.signal <- min( xl.cells[ xl.cells$off_data, 'data' ], na.rm = TRUE )
                max.signal <- max( xl.cells[ xl.cells$off_data, 'data' ], na.rm = TRUE ) + 1E-7
                xl.cells[ xl.cells$off_data, 'q.off' ] <- ( xl.cells[ xl.cells$off_data, 'data' ] - min.signal ) /
                                  ( max.signal - min.signal )
                xl.cells[ xl.cells$off_data, 'fill.color' ] <- sapply( xl.cells[ xl.cells$off_data, 'q.off' ], function( k ) do.call( rgb, as.list( cells.pal( k )[1,]/255 ) ) )
            }
        }
        header_style <- CellStyle( bimap.workbook ) + Fill( 'lightblue' ) +
                        Font( bimap.workbook, isItalic = TRUE, isBold = TRUE )
        info_style <- CellStyle( bimap.workbook ) +
                       Font( bimap.workbook, isItalic = TRUE )

        message( 'Samples data writing...' )
        header_rows <- createRow( wsh, rowIndex = 1:length(samp_cols) )
        mapply( setCellValue, createCell( header_rows, colIndex = col_offset ), samp_cols )
        lapply( getCells( header_rows, colIndex = col_offset ), setCellStyle, header_style )
        for ( i in 1:length(samp_cols) ) {
            col_row <- getRows( wsh, i )
            mapply( setCellValue, createCell( col_row, colIndex = col_offset + samples$order ),
                                  samples[,samp_cols[[i]]])
        }
        grid_rows <- getRows( wsh, 1:length(samp_cols) )
        lapply( getCells( grid_rows, colIndex = subset( samples, clu.order > 1 )$order + col_offset ),
                setCellStyle, info_style )
        lapply( getCells( grid_rows, colIndex = 1 + col_offset ),
                setCellStyle, info_style + Border( color = grid.col,
                                     position = c( 'LEFT' ), pen = c( 'BORDER_THICK' ) ) )
        lapply( getCells( grid_rows, colIndex = subset( samples, clu.order == 1 & order > 1 )$order + col_offset ),
                setCellStyle, info_style + Border( color = grid.col,
                                     position = c( 'LEFT' ), pen = c( 'BORDER_THIN' ) ) )
        message( 'Protein data writing...' )
        rows <- getRows( wsh, proteins$order + row_offset )
        header_row <- getRows( wsh, rowIndex = row_offset )
        mapply( setCellValue, createCell( header_row, colIndex = 1:length(prot_cols) ), prot_cols )
        lapply( getCells( header_row, colIndex = 1:length(prot_cols) ), setCellStyle, header_style )
        for ( i in 1:length(prot_cols) ) {
            mapply( setCellValue, createCell( rows, colIndex = i ),
                                  proteins[,prot_cols[[i]]] )
        }
        lapply( getCells( getRows( wsh, subset( proteins, clu.order > 1 )$order + row_offset ),
                          colIndex = 1:length(prot_cols) ),
                setCellStyle, info_style )
        lapply( getCells( getRows( wsh, 1 + row_offset ),
                          colIndex = 1:length(prot_cols) ),
                setCellStyle, info_style + Border( color = grid.col,
                                     position = c( 'TOP' ), pen = c( 'BORDER_THICK' ) ) )
        lapply( getCells( getRows( wsh, subset( proteins, clu.order == 1 & order > 1 )$order + row_offset ),
                          colIndex = 1:length(prot_cols) ),
                setCellStyle, info_style + Border( color = grid.col,
                                     position = c( 'TOP' ), pen = c( 'BORDER_THIN' ) ) )

        message( 'Applying cell styles...' )
        #print( head( xl.cells ))
        ddply( xl.cells, c( 'is_bait', 'off_data', 'fill.color', 'border.color',
                            'border.left', 'border.right', 'border.top', 'border.bottom' ),
            function( sub.cells ) {
                # construct style
                style_desc <- sub.cells[ 1, ]
                cs <- CellStyle( bimap.workbook )
                if ( !is.na( style_desc$fill.color ) ) {
                    cs <- cs + Fill( foregroundColor = sub.cells[ 1, 'fill.color' ] )
                }
                cs <- cs + Font( bimap.workbook, isBold = style_desc$is_bait,
                                                 isItalic = style_desc$off_data )
                border_pos <- c()
                if ( style_desc$border.left ) border_pos <- c( border_pos, c( 'LEFT' ) )
                if ( style_desc$border.right ) border_pos <- c( border_pos, c( 'RIGHT' ) )
                if ( style_desc$border.top ) border_pos <- c( border_pos, c( 'TOP' ) )
                if ( style_desc$border.bottom ) border_pos <- c( border_pos, c( 'BOTTOM' ) )
                if ( length( border_pos) > 0 ) {
                  cs <- cs + Border( color = style_desc$border.color,
                                     position = border_pos, pen = c( 'BORDER_THIN' ) )
                }
                # apply style
                prots.pos <- unique( sub.cells$order.row )
                for ( prot.pos in prots.pos ) {
                    prot.row <- getRows( wsh, rowIndex= prot.pos + row_offset )
                    cells <- getCells( prot.row, colIndex= subset( sub.cells, order.row == prot.pos )$order.col + col_offset )
                    lapply( cells, setCellStyle, cs )
                }
            
            return ( NULL )
        } )
        setColumnWidth( wsh, col_offset + 1:nrow(samples), col.width )

        message( 'On-Blocks signal comments...' ) 
        for ( i in 1:nrow(blocks) ) {
            block <- blocks[i,]
            first.cell <- subset( xl.cells, proteins.cluster == block$proteins.cluster &
                                            samples.cluster == block$samples.cluster &
                                            clu.order.row == 1 & clu.order.col == 1 )
            signal_txt <- paste( format( block$signal.mean, digits=3 ), "\u0b1",
                                 format( block$signal.sd, digits=3 ), sep="" )
            first_row <- getRows( wsh, rowIndex = first.cell[ 1, 'order.row' ] + row_offset )
            lapply( getCells( first_row, colIndex = first.cell[ 1, 'order.col' ] + col_offset ),
                    createCellComment, string = signal_txt )
        }
        message( 'Merging sample info columns' )
        for ( icol in 1:length(samp_cols) ) {
            ddply( samples, c( "samples.cluster", samp_cols[[icol]] ), function( merged_cells ) {
                addMergedRegion( wsh, icol, icol,
                                 col_offset + min( merged_cells$order ),
                                 col_offset + max( merged_cells$order ) )
            } )
        }

        message( "Outer frame..." )
        lapply( createCell( createRow( wsh, row_offset + nrow( proteins ) + 1 ),
                            colIndex = 1:(col_offset+nrow(samples) ) ),
                setCellStyle, CellStyle( bimap.workbook ) + Border( grid.col, position = c( 'TOP' ), pen = c( 'BORDER_THICK' ) ) )
        lapply( createCell( getRows( wsh, 1:( row_offset + nrow( proteins ) ) ), 
                            colIndex = col_offset+nrow(samples) + 1 ),
                setCellStyle, CellStyle( bimap.workbook ) + Border( grid.col, position = c( 'LEFT' ), pen = c( 'BORDER_THICK' ) ) )
    })
    return ( bimap.workbook )
}
