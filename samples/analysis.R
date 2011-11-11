# Sample AP-MS data analysis script
# 
# Author: Alexey Stukalov (astukalov@cemm.oeaw.ac.at)
###############################################################################

# base path to R scripts and BI-MAP (update to actual location)
base_path <- "/home/<user>/projects/bi-map"

# path to BI-MAP scripts
bimap_scripts_path <- file.path( base_path, "scripts/R" )
# path to data files
data_path <- file.path( base_path, "samples/tip49" )
# path to output files
results_path <- file.path( base_path, "results" )

# define path to RBIMAP library (update to actual location)
RBIMAP.libpath <- file.path( base_path, "build/release/src/rcpp" )

# load R BIMAP functions
source( file.path( bimap_scripts_path, "BIMAP.R" ) )
source( file.path( bimap_scripts_path, "BIMAP_plot.R" ) )
source( file.path( bimap_scripts_path, "BIMAP_graphml.R" ) )

# load AP-MS data
bimap.data <- BIMAP.msdata.load( data_path, format = 'CSV' )

# load calculated BI-MAP MCMC walk
# (run samples/run_bimap_tip49.sh to get run)
bimap.walk <- BIMAP.mcmcwalk.load( file.path( results_path, 'tip49_walk.xml.bz2' ),
    0.02, 0.02 )

# or execute BI-MAP

# generate a dataframe of top n clustering models sorting according to their posterior probability
bimap.best_clusterings.df <- BIMAP.mcmcwalk_biclusterings_order( bimap.walk, n = 10 )
# extract the model with max posterior probability
bimap.best_clustering <- BIMAP.mcmcwalk.extract_biclustering( bimap.walk, bimap.best_clusterings.df[1,'clustering.serial'] )

# render the matrix view of BI-MAP model to PDF file
pdf( file= file.path( results_path, 'tip49_blocks.pdf' ),
     title="TIP49a/b", width=34, height=60 )
BIMAP.plot( bimap.best_clustering,
    bimap.data,
    #proteins.info = apms.proteins_info,
    sample_name_col = 'condition',
    protein_name_col = 'short_id',
    protein_description_col = 'description',
    show.abundance.labels = FALSE,
    show.protein_ac = TRUE, show.sample_id = TRUE, show.msruns = TRUE,
    show.measurements = TRUE,
    aspect = 3.0
)
dev.off()

