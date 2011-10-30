#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "ParametersReader.h"

namespace po = boost::program_options;


BIMAPIOParams::BIMAPIOParams()
    : minCrossClusRefCount( 1 )
    , minObjectsPtnRefCount( 1 )
    , minProbesPtnRefCount( 1 )
    , csvColumnSeparator( '\t' )
{
}

bool BIMAPParamsRead(
    int argc, char* argv[],
    ChessboardBiclusteringHyperPriors&     hyperpriors,
    GibbsSamplerParams&             gibbsParams,
    TurbineCascadeParams&           cascadeParams,
    CellSignalParams&               signalParams,
    PrecomputedDataParams&          precomputedDataParams,
    ChessboardBiclusteringPriors&   priors,
    BIMAPSampleCollectorParams&     collectorParams,
    BIMAPIOParams&                  ioParams
){
    // Declare the supported options
    po::positional_options_description p;
    p.add( "input_file", -1 );

    po::options_description basic_cmdline_options_desc( "Program options" );

    basic_cmdline_options_desc.add_options()
        ( "help", "options description" )
        ( "config_file", po::value< std::string >(),
            "configuration (INI) file" )
        ( "input_file", po::value< std::string >( &ioParams.dataFilename ),
            "input data <serialized OPAData>" )

        ( "proteins_file", po::value< std::string >( &ioParams.proteinsFilename ),
            "proteins table (AC, sequence length) in CSV (tab-separated) format, required if not input_file provided" )
        ( "exp_design_file", po::value< std::string >( &ioParams.expDesignFilename ),
            "experimental design table (bait AC, sample ID, MS run ID, MS run multiplier) in CSV (tab-separated) format, required if not input_file provided" )
        ( "measurements_file", po::value< std::string >( &ioParams.expDesignFilename ),
            "MS measurements table (MS run ID, protein AC, Spectral Counts, [Peptide Counts]) in CSV (tab-separated) format, required if not input_file provided" )
        ( "csv_column_separator", po::value<char>( &ioParams.csvColumnSeparator )->default_value( ioParams.csvColumnSeparator ),
            "CSV files column separator" )

        ( "output_file", po::value< std::string >( &ioParams.outputFilename ),
            "output filename <serialized BIMAP walk>" );

    po::options_description cascade_options_desc( "EE sampler cascade parameters" );
    cascade_options_desc.add_options()
        ( "cascade.levels", po::value< size_t >( &cascadeParams.levelsCount )->default_value( cascadeParams.levelsCount ), 
            "number of temperature levels in the cascade" )
        ( "cascade.turbines", po::value< size_t >( &cascadeParams.turbinesCount )->default_value( cascadeParams.turbinesCount ),
            "number of turbines in the cascade (by default equals to MPI communicator rank)" )
        ( "cascade.turbines_per_level", po::value< std::string >(), 
            "comma-separated number of turbines per level" )
        ( "cascade.turbine_capacity", po::value<size_t>( &cascadeParams.turbineParams.maxParticles )->default_value( cascadeParams.turbineParams.maxParticles ),
            "maximum number of particles stored in turbine's disk" )
        ( "cascade.temperature_multiplier", po::value<energy_type>( &cascadeParams.temperatureMultiplier )->default_value( cascadeParams.temperatureMultiplier ),
            "temperature ratio between subsequent hotter and colder turbines" )
        ( "cascade.burnin_iterations", po::value<size_t>( &cascadeParams.burnInIterations )->default_value( cascadeParams.burnInIterations ),
            "number of cascade burn-in iterations (after all turbines burned-in)" )
        ( "cascade.ladder_adjust_period", po::value<size_t>( &cascadeParams.ladderAdjustPeriod )->default_value( cascadeParams.ladderAdjustPeriod ),
            "period (# of leading turbine's iterations) the energy ladder is adjusted during burn-in" )
        ( "cascade.ee_jump_rate", po::value<prob_t>( &cascadeParams.turbineParams.eeJumpRate )->default_value( cascadeParams.turbineParams.eeJumpRate ),
            "rate of equi-energy jumps" )
        ( "cascade.ee_jump_energy_tolerance", po::value<energy_type>( &cascadeParams.turbineParams.eeJumpEnergyTolerance )->default_value( cascadeParams.turbineParams.eeJumpEnergyTolerance ),
            "fraction of total current energy range, that defines energy ring of possibe equi-energy jump destinations" )
        ( "cascade.ee_max_energy_quantile", po::value<prob_t>( &cascadeParams.turbineParams.eeMaxEnergyQuantile )->default_value( cascadeParams.turbineParams.eeMaxEnergyQuantile ),
            "quantile of energy distribution used to define the energy range" )
        ( "cascade.ee_ladder_energy_quantile", po::value<prob_t>( &cascadeParams.turbineParams.eeLadderEnergyQuantile )->default_value( cascadeParams.turbineParams.eeLadderEnergyQuantile ),
            "quantile of energy distribution used to define energy threshold for step of energy ladder (using energies distributions on the previous step)" )
        ( "cascade.generate_rate", po::value<prob_t>( &cascadeParams.turbineParams.generateRate )->default_value( cascadeParams.turbineParams.generateRate ),
            "rate of particles generation to be sent to hotter turbines" )
        ( "cascade.detention_iterations", po::value<size_t>( &cascadeParams.turbineParams.detentionIterations )->default_value( cascadeParams.turbineParams.detentionIterations ),
            "# of iterations the generated particle need to do, before it's put to the energy disk" )
        ( "cascade.turbine_snapshot_period", po::value< size_t >( &cascadeParams.turbineParams.particleSnapshotPeriod )->default_value( cascadeParams.turbineParams.particleSnapshotPeriod ),
            "period (# of iterations) when current turbine's particle is stored to turbine's disk" )
        ( "cascade.turbine_status_broadcast_period", po::value< size_t >( &cascadeParams.broadcastStatusPeriod )->default_value( cascadeParams.broadcastStatusPeriod ),
            "period (# of iterations) that each turbine broadcasts its status to all nodes" )

        ( "samples", po::value< size_t >( &collectorParams.walkSamples )->default_value( collectorParams.walkSamples ),
            "number of samples to collect" )
        ( "storage_period", po::value< size_t >( &cascadeParams.turbineParams.particleSamplingPeriod )->default_value( cascadeParams.turbineParams.particleSamplingPeriod ),
            "period (# of iterations) before new sample is collected" )
        ( "priors_storage_period", po::value< size_t >( &collectorParams.priorsStoragePeriod )->default_value( collectorParams.priorsStoragePeriod ),
            "period (# of cluster samples collected) before new sample of prior parameters is collected" );

    po::options_description prior_params_desc( "Prior model distribution parameters" );
    prior_params_desc.add_options()
        ( "prior.obj_clustering_concentration", po::value<prob_t>( &priors.objectClustering.concentration )->default_value( priors.objectClustering.concentration ),
            "'concentration' parameter for Pitman-Yor prior of objects to clusters distribution" )
        ( "prior.obj_clustering_discount", po::value<prob_t>( &priors.objectClustering.discount )->default_value( priors.objectClustering.discount ),
            "'discount' parameter for Pitman-Yor prior of objects to clusters distribution" )

        ( "prior.probe_clustering_concentration", po::value<prob_t>( &priors.probeClustering.concentration )->default_value( priors.probeClustering.concentration ),
            "'concentration' parameter for Pitman-Yor prior of probes to clusters distribution" )
        ( "prior.probe_clustering_discount", po::value<prob_t>( &priors.probeClustering.discount )->default_value( priors.probeClustering.discount ),
            "'discount' parameter for Pitman-Yor prior of probes to clusters distribution" )

        ( "prior.ccluster_enabled", po::value<prob_t>( &priors.cellEnablementProb )->default_value( priors.cellEnablementProb ),
            "prior probability that cross-cluster is enabled" )
        ( "prior.objects_clusters_off_probes_cluster_rate", po::value<prob_t>( &priors.objectsCluOffProbesCluRate )->default_value( priors.objectsCluOffProbesCluRate ),
            "failure rate of geometrical distribution, which is prior expectation of the number of objects clusters related to probes cluster via probes baits" )
        ( "prior.probes_clusters_off_objects_cluster_rate", po::value<prob_t>( &priors.probesCluOffObjectsCluRate )->default_value( priors.probesCluOffObjectsCluRate ),
            "failure rate of geometrical distribution, which is prior expectation of the number of probes clusters related to objects cluster via probes baits" )
        ( "prior.objects_clusters_per_probes_cluster_rate", po::value<prob_t>(),
            "success rate of geometrical distribution, which is prior expectation of the number of objects clusters related to probes cluster via probes baits" )
        ( "prior.probes_clusters_per_objects_cluster_rate", po::value<prob_t>(),
            "success rate of geometrical distribution, which is prior expectation of the number of probes clusters related to objects cluster via probes baits" )
        ( "prior.false_hits", po::value<prob_t>( &priors.noise.failures )->default_value( priors.noise.failures ),
            "prior counts of false protein discoveries (in disabled cross-clusters)" )
        ( "prior.true_misses", po::value<prob_t>( &priors.noise.successes )->default_value( priors.noise.successes ),
            "prior counts of true non-detections of protein (in disabled cross-clusters)" )
        ( "prior.object_multiple_rate", po::value<prob_t>( &priors.objectMultipleRate )->default_value( priors.objectMultipleRate ),
            "failure rate of object multiple (geometrical) distribution" )
            ;

    po::options_description hyperprior_params_desc( "Hyperprior model parameters" );
    hyperprior_params_desc.add_options()
        ( "hyperprior.signal_mean_mean", po::value<signal_t>( &hyperpriors.signalHyperprior.meanMean )->default_value( hyperpriors.signalHyperprior.meanMean ),
            "prior for signal's normal distribution -- mean parameter for normal-scaled inverse Gamma" )
        ( "hyperprior.signal_mean_var_scale", po::value<signal_t>( &hyperpriors.signalHyperprior.meanVarScale )->default_value( hyperpriors.signalHyperprior.meanVarScale ),
            "prior for signal's normal distribution -- mean variance scale parameter for normal-scaled inverse Gamma" )
        ( "hyperprior.signal_var_scale", po::value<signal_t>( &hyperpriors.signalHyperprior.varDistrib.scale )->default_value( hyperpriors.signalHyperprior.varDistrib.scale ),
            "prior for signal's normal distribution -- scale of variance for normal-scaled inverse Gamma" )
        ( "hyperprior.signal_var_shape", po::value<signal_t>( &hyperpriors.signalHyperprior.varDistrib.shape )->default_value( hyperpriors.signalHyperprior.varDistrib.shape ),
            "prior for signal's normal distribution -- shape of variance for normal-scaled inverse Gamma" )
        ;

    po::options_description signal_params_desc( "Signal parameters" );
    signal_params_desc.add_options()
        ( "signal.seq_length_factor", po::value<log_prob_t>( &signalParams.sequenceLengthFactor )->default_value( signalParams.sequenceLengthFactor ),
            "the power sequence length is raised to (signal is divided by L^p for signal normalization)" )
        ( "signal.signal_counts_shape", po::value<log_prob_t>( &signalParams.scShape )->default_value( signalParams.scShape ),
            "[0,1), shape parameter of Lagrangian Poisson TRUE-POSITIVE spectral counts distribution (0 = classic Poisson)" )
        ( "signal.zero_inflation_factor", po::value<log_prob_t>( &signalParams.zeroInflationFactor )->default_value( signalParams.zeroInflationFactor ),
            "PDF(factor*mode()) is added to PDF(0) (0 = disable)" )
        ( "signal.preliminary_signal_prior_shape", po::value<double>( &signalParams.preliminarySignalPrior.shape )->default_value( signalParams.preliminarySignalPrior.shape ),
            "shape (alpha) parameter of Gamma preliminary signals prior" )
        ( "signal.preliminary_signal_prior_scale", po::value<double>( &signalParams.preliminarySignalPrior.scale )->default_value( signalParams.preliminarySignalPrior.scale ),
            "scale (beta) parameter of Gamma preliminary signals prior" )
        ;

    po::options_description precomputed_params_desc( "Precomputed data parameters" );
    precomputed_params_desc.add_options()
        ( "precomputed.object_freq_threshold", po::value<double>( &precomputedDataParams.objectFreqThreshold )->default_value( precomputedDataParams.objectFreqThreshold ),
            "maximal frequence of object occurences in probes to be still used for probes neighbourhood checking" )
        ( "precomputed.probe_freq_threshold", po::value<double>( &precomputedDataParams.probeFreqThreshold )->default_value( precomputedDataParams.probeFreqThreshold ),
            "maximal frequence of object occurences in probes to be still used for probes neighbourhood checking" )
    ;

    po::options_description gibbs_sampler_params_desc( "Gibbs sampler parameters" );
    gibbs_sampler_params_desc.add_options()
        ( "sampler.signal_rate", po::value<prob_t>( &gibbsParams.signalRate )->default_value( gibbsParams.signalRate ),
            "rate of cross-cluster signal sampling steps" )
        ( "sampler.object_multiple_rate", po::value<prob_t>( &gibbsParams.objectMultipleRate )->default_value( gibbsParams.objectMultipleRate ),
            "rate of object multiple (stoichiometry) sampling steps" )
            ( "sampler.ccluster_flip_rate", po::value<prob_t>( &gibbsParams.crossClusterFlipRate )->default_value( gibbsParams.crossClusterFlipRate ),
            "rate of cross-cluster probe flipping sampling steps" )
        ( "sampler.object_membership_rate", po::value<prob_t>( &gibbsParams.objectMembershipRate )->default_value( gibbsParams.objectMembershipRate ),
            "rate of single object membership changing sampling steps" )
        ( "sampler.probe_membership_rate", po::value<prob_t>( &gibbsParams.probeMembershipRate )->default_value( gibbsParams.probeMembershipRate ),
            "rate of single probe membership changing sampling steps" )
        ( "sampler.object_split_merge_rate", po::value<prob_t>( &gibbsParams.objectsSplitMergeRate )->default_value( gibbsParams.objectsSplitMergeRate ),
            "rate of objects clusters split/merge sampling steps" )
        ( "sampler.probe_split_merge_rate", po::value<prob_t>( &gibbsParams.probesSplitMergeRate )->default_value( gibbsParams.probesSplitMergeRate ),
            "rate of objects clusters split/merge sampling steps" )
        ( "sampler.ccluster_resamples", po::value<size_t>( &gibbsParams.crossClusterResamples )->default_value( gibbsParams.crossClusterResamples ),
            "number of signal sampling steps after cross-cluster flipped" )

        ( "sampler.object_cluster_existing_samples", po::value<size_t>( &gibbsParams.objectClusterParams.existingClusterSamples )->default_value( gibbsParams.objectClusterParams.existingClusterSamples ),
            "number of signal sampling steps (for each cross-cluster) after objects cluster contents changed" )
        ( "sampler.object_cluster_new_samples", po::value<size_t>( &gibbsParams.objectClusterParams.newClusterSamples )->default_value( gibbsParams.objectClusterParams.newClusterSamples ),
            "number of signal sampling steps (for each cross-cluster) for new objects cluster" )
        ( "sampler.object_membership_samples", po::value<size_t>( &gibbsParams.objectClusterParams.existingClusterSamples )->default_value( gibbsParams.objectClusterParams.existingClusterSamples ),
            "number of samples for assigning object to existing cluster" )
        ( "sampler.object_new_cluster_samples", po::value<size_t>( &gibbsParams.objectClusterParams.newClusterSamples )->default_value( gibbsParams.objectClusterParams.newClusterSamples ),
            "number of samples for assigning object to a new (singleton) cluster" )
        ( "sampler.object_cluster_split_samples", po::value<size_t>( &gibbsParams.objectsSplitMergeParams.splitLaunchSamples )->default_value( gibbsParams.objectsSplitMergeParams.splitLaunchSamples ),
            "number of signal sampling steps (for each cross-cluster) for split object cluster launch probe" )
        ( "sampler.object_cluster_merge_samples", po::value<size_t>( &gibbsParams.objectsSplitMergeParams.mergeLaunchSamples )->default_value( gibbsParams.objectsSplitMergeParams.mergeLaunchSamples ),
            "number of signal sampling steps (for each cross-cluster) for merged object clusters launch probe" )

        ( "sampler.probe_cluster_existing_samples", po::value<size_t>( &gibbsParams.probeClusterParams.existingClusterSamples )->default_value( gibbsParams.probeClusterParams.existingClusterSamples ),
            "number of signal sampling steps (for each cross-cluster) after probes cluster contents changed" )
        ( "sampler.probe_cluster_new_samples", po::value<size_t>( &gibbsParams.probeClusterParams.newClusterSamples )->default_value( gibbsParams.probeClusterParams.newClusterSamples ),
            "number of signal sampling steps (for each cross-cluster) for new probes cluster" )
        ( "sampler.probe_membership_samples", po::value<size_t>( &gibbsParams.probeClusterParams.existingClusterSamples )->default_value( gibbsParams.probeClusterParams.existingClusterSamples ),
            "number of samples for assigning probe to existing cluster" )
        ( "sampler.probe_new_cluster_samples", po::value<size_t>( &gibbsParams.probeClusterParams.newClusterSamples )->default_value( gibbsParams.probeClusterParams.newClusterSamples ),
            "number of samples for assigning probe to a new (singleton) cluster" )
        ( "sampler.probe_cluster_split_samples", po::value<size_t>( &gibbsParams.probesSplitMergeParams.splitLaunchSamples )->default_value( gibbsParams.probesSplitMergeParams.splitLaunchSamples ),
            "number of signal sampling steps (for each cross-cluster) for split probe cluster launch probe" )
        ( "sampler.probe_cluster_merge_samples", po::value<size_t>( &gibbsParams.probesSplitMergeParams.mergeLaunchSamples )->default_value( gibbsParams.probesSplitMergeParams.mergeLaunchSamples ),
            "number of signal sampling steps (for each cross-cluster) for merged probe clusters launch probe" )

        ( "sampler.mean_object_rank", po::value<size_t>( &gibbsParams.meanObjectRank )->default_value( gibbsParams.meanObjectRank ),
            "fuzziness of the ranks to use for split/merge and object membership moves (not the very best, but rank-d best)" )
        ( "sampler.mean_probe_rank", po::value<size_t>( &gibbsParams.meanProbeRank )->default_value( gibbsParams.meanProbeRank ),
            "fuzziness of the ranks to use for split/merge and probe membership moves (not the very best, but rank-d best)" )

        ( "sampler.priors_update_period", po::value<size_t>( &gibbsParams.priorsUpdatePeriod )->default_value( gibbsParams.priorsUpdatePeriod ),
            "period (# of cclustering sampling iterations) of priors updating sampling step" )
        ;

    po::options_description log_params_desc( "Logging parameters" );
    gibbs_sampler_params_desc.add_options()
        ( "log.samples_reporting_period", po::value<prob_t>( &gibbsParams.signalRate )->default_value( gibbsParams.signalRate ),
            "period (# of samples collected) to output info about collected samples" )
        ( "log.samples_max_reporting_delay", po::value<prob_t>( &gibbsParams.crossClusterFlipRate )->default_value( gibbsParams.crossClusterFlipRate ),
            "maximal time without reporting about collected samples (in case it takes too long)" )
        ;

    po::options_description output_params_desc( "Output parameters" );
    gibbs_sampler_params_desc.add_options()
        ( "output.min_cross_clus_ref_counts", po::value<size_t>( &ioParams.minCrossClusRefCount )->default_value( ioParams.minCrossClusRefCount ),
            "minimal number of times cross clustering should appear in random walk, less frequent are removed" )
        ( "output.min_objects_ptn_ref_counts", po::value<size_t>( &ioParams.minObjectsPtnRefCount )->default_value( ioParams.minObjectsPtnRefCount ),
            "minimal number of times objects partition should appear in random walk, less frequent are removed" )
        ( "output.min_probes_ptn_ref_counts", po::value<size_t>( &ioParams.minProbesPtnRefCount )->default_value( ioParams.minProbesPtnRefCount ),
            "minimal number of times probes partition should appear in random walk, less frequent are removed" )
        ;

    po::options_description all_cmdline_options_desc = basic_cmdline_options_desc;
    all_cmdline_options_desc
        .add( cascade_options_desc )
        .add( prior_params_desc )
        .add( hyperprior_params_desc )
        .add( signal_params_desc )
        .add( precomputed_params_desc )
        .add( gibbs_sampler_params_desc )
        .add( log_params_desc )
        .add( output_params_desc )
    ;

    po::variables_map opt_map;
    po::store( po::command_line_parser( argc, argv )
                .options( all_cmdline_options_desc )
                .positional( p ).run(), opt_map );

    if ( opt_map.count( "help" ) ) {
        std::cout << all_cmdline_options_desc;
        return ( false );
    }
    po::options_description all_config_options_desc( "Configuration file parameters" );
    all_config_options_desc
        .add( cascade_options_desc )
        .add( prior_params_desc )
        .add( hyperprior_params_desc )
        .add( signal_params_desc )
        .add( precomputed_params_desc )
        .add( gibbs_sampler_params_desc )
        .add( log_params_desc )
        .add( output_params_desc )
    ;

    if ( opt_map.count( "config_file" ) ) {
        po::store( po::parse_config_file<char>( opt_map["config_file"].as<std::string>().c_str(),
                                                all_config_options_desc ),
                    opt_map );
    }
    po::notify( opt_map );

    // use _per_, if available and _off_ not provided
    if ( opt_map.count( "prior.probes_clusters_off_objects_cluster_rate" ) == 0
         && opt_map.count( "prior.probes_clusters_per_objects_cluster_rate" ) > 0 ) {
        priors.probesCluOffObjectsCluRate = 1.0 - opt_map[ "prior.probes_clusters_per_objects_cluster_rate" ].as<prob_t>();
    }
    if ( opt_map.count( "prior.objects_clusters_off_probes_cluster_rate" ) == 0
         && opt_map.count( "prior.objects_clusters_per_probes_cluster_rate" ) > 0 ) {
        priors.objectsCluOffProbesCluRate = 1.0 - opt_map[ "prior.objects_clusters_per_probes_cluster_rate" ].as<prob_t>();
    }

    // check output path existence
    boost::filesystem::path out_path( opt_map[ "output_file" ].as<std::string>() );
    if ( !boost::filesystem::exists( out_path.parent_path() ) ) {
        THROW_EXCEPTION( std::invalid_argument, "Output folder doesn't exist: " << out_path.parent_path() );
    }
    if ( !boost::filesystem::is_directory( out_path.parent_path() ) ) {
        THROW_EXCEPTION( std::invalid_argument, "Output path is not a folder: " << out_path.parent_path() );
    }

    // update internals of the options
    priors.updateCachedDistributions();

    signalParams.preliminarySignalPrior = GammaDistribution(
            signalParams.preliminarySignalPrior.shape,
            signalParams.preliminarySignalPrior.scale
    );

    if ( opt_map.count( "input_file" ) == 0
        && (  opt_map.count( "proteins_file" ) == 0
           || opt_map.count( "exp_design_file" ) == 0
           || opt_map.count( "measurements_file" ) == 0
    ) ){
        THROW_RUNTIME_ERROR( "No input file(s) defined" );
    }

    if ( opt_map.count( "cascade.turbines_per_level" ) > 0 ) {
        std::stringstream ss( opt_map["cascade.turbines_per_level"].as<std::string>() );
        cascadeParams.turbinesPerLevel.clear();
        cascadeParams.turbinesCount = 0;
        int number;
        while (ss >> number) {
            if ( number > 0 ) cascadeParams.turbinesPerLevel.push_back( (size_t)number );
            else THROW_RUNTIME_ERROR( number << " is not positive" );
            if (ss.peek() == ',') ss.ignore();
            cascadeParams.turbinesCount += number;
        }
        cascadeParams.levelsCount = cascadeParams.turbinesPerLevel.size();
    }
#if 0 /// @todo print parameters used
    LOG_INFO( "Parameters used:" );
    for ( po::variables_map::const_iterator it = opt_map.begin(); it != opt_map.end(); ++it ) {
        LOG_INFO( it->first << "=" << it->second.as<std::string>() );
    }
#endif
    return ( true );
}

