#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "ParametersReader.h"

namespace po = boost::program_options;

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
          "configuration (.INI) file" )
        ( "input_file", po::value< std::string >( &ioParams.dataFilename ),
          "input data <serialized OPAData>" )

        ( "proteins_file", po::value< std::string >( &ioParams.proteinsFilename ),
          "proteins table (AC, sequence length) in CSV (tab-separated) format.\nRequired iff no --input_file provided" )
        ( "exp_design_file", po::value< std::string >( &ioParams.expDesignFilename ),
          "experimental design table (bait AC, sample ID, MS run ID, MS run multiplier) in CSV (tab-separated) format.\nRequired iff no --input_file provided" )
        ( "measurements_file", po::value< std::string >( &ioParams.measurementsFilename ),
          "MS measurements table (MS run ID, protein AC, Spectral Counts, [Peptide Counts]) in CSV (tab-separated) format.\nRequired iff no --input_file provided" )
        ( "csv_column_separator", po::value<char>( &ioParams.csvColumnSeparator )->default_value( ioParams.csvColumnSeparator ),
          "CSV files column separator (TAB used by default)" )

        ( "map_baits_to_preys", po::value<bool>( &ioParams.mapBaitsToObjects )->default_value( ioParams.mapBaitsToObjects ),
          "if true, baits are regarded as proteins and experimental design likelihood component is calculated" )
        ( "objects_universe_size", po::value<size_t>( &ioParams.objectsUniverseSize )->default_value( ioParams.objectsUniverseSize ),
          "total number of objects in the universe (proteome size)" )

        ( "output_file", po::value< std::string >( &ioParams.outputFilename ),
          "output filename <serialized BI-MAP MCMC walk>" );

    po::options_description eesampler_options_desc( "EE sampler parameters" );
    eesampler_options_desc.add_options()
        ( "eesampler.levels", po::value< size_t >( &cascadeParams.levelsCount )->default_value( cascadeParams.levelsCount ), 
          "number of temperature levels" )
        ( "eesampler.simple_samplers", po::value< size_t >( &cascadeParams.turbinesCount )->default_value( cascadeParams.turbinesCount ),
          "number of simple (Gibbs) samplers (by default equals to MPI communicator rank)" )
        ( "eesampler.simple_samplers_per_level", po::value< std::string >(), 
          "numbers (comma-separated) of simple (Gibbs) samplers per level" )
        ( "eesampler.sampler_storage_capacity", po::value<size_t>( &cascadeParams.turbineParams.maxParticles )->default_value( cascadeParams.turbineParams.maxParticles ),
          "maximum number of intermediate models cached by each simple sampler" )
        ( "eesampler.temperature_multiplier", po::value<energy_type>( &cascadeParams.temperatureMultiplier )->default_value( cascadeParams.temperatureMultiplier ),
          "temperature ratio between neighbouring levels" )
        ( "eesampler.burnin_iterations", po::value<size_t>( &cascadeParams.burnInIterations )->default_value( cascadeParams.burnInIterations ),
          "number of burn-in iterations (# of iterations of the leading simple turbine)" )
        ( "eesampler.ladder_adjust_period", po::value<size_t>( &cascadeParams.ladderAdjustPeriod )->default_value( cascadeParams.ladderAdjustPeriod ),
          "period (# of leading simple sampler iterations) after which the energy levels of the whole ladder are adjusted during burn-in" )
        ( "eesampler.ee_jump_rate", po::value<prob_t>( &cascadeParams.turbineParams.eeJumpRate )->default_value( cascadeParams.turbineParams.eeJumpRate ),
          "rate of equi-energy jumps" )
        ( "eesampler.ee_jump_energy_tolerance", po::value<energy_type>( &cascadeParams.turbineParams.eeJumpEnergyTolerance )->default_value( cascadeParams.turbineParams.eeJumpEnergyTolerance ),
          "relative difference of current and considered destination energies, a fraction of current energy range (E_max - E_min)" )
        ( "eesampler.max_energy_quantile", po::value<prob_t>( &cascadeParams.turbineParams.eeMaxEnergyQuantile )->default_value( cascadeParams.turbineParams.eeMaxEnergyQuantile ),
          "quantile of energies distribution that is used to define E_max" )
        ( "eesampler.step_energy_quantile", po::value<prob_t>( &cascadeParams.turbineParams.eeLadderEnergyQuantile )->default_value( cascadeParams.turbineParams.eeLadderEnergyQuantile ),
          "quantile of the energy distribution from the lower energy step that defines energy threshold E_step for the next step of energy ladder" )
        ( "eesampler.generate_rate", po::value<prob_t>( &cascadeParams.turbineParams.generateRate )->default_value( cascadeParams.turbineParams.generateRate ),
          "rate of particles generation that are sent to hotter samplers" )
        ( "eesampler.detention_iterations", po::value<size_t>( &cascadeParams.turbineParams.detentionIterations )->default_value( cascadeParams.turbineParams.detentionIterations ),
         "# of iterations a generated particle needs to stay at hotter turbine, before it could be picked by EE jump" )
        ( "eesampler.caching_period", po::value< size_t >( &cascadeParams.turbineParams.particleSnapshotPeriod )->default_value( cascadeParams.turbineParams.particleSnapshotPeriod ),
          "period (# of iterations) of putting current simple sampler's model in cache" )
        ( "eesampler.broadcasting_period", po::value< size_t >( &cascadeParams.broadcastStatusPeriod )->default_value( cascadeParams.broadcastStatusPeriod ),
          "period (# of iterations) each simple sampler broadcasts its status to all nodes" )

        // basic sampler parameters
        ( "samples", po::value< size_t >( &collectorParams.walkSamples )->default_value( collectorParams.walkSamples ),
            "number of samples to collect" )
        ( "storage_period", po::value< size_t >( &cascadeParams.turbineParams.particleSamplingPeriod )->default_value( cascadeParams.turbineParams.particleSamplingPeriod ),
          "period (# of iterations) of collecting new model sample" )
        ( "priors_storage_period", po::value< size_t >( &collectorParams.priorsStoragePeriod )->default_value( collectorParams.priorsStoragePeriod ),
          "period (# of cluster samples collected) of collecting new prior parameters sample" );

    po::options_description prior_params_desc( "Prior model distribution parameters" );
    prior_params_desc.add_options()
        ( "prior.obj_clustering_concentration", po::value<prob_t>( &priors.objectClustering.concentration )->default_value( priors.objectClustering.concentration ),
          "'concentration' parameter for Pitman-Yor prior of objects clustering" )
        ( "prior.obj_clustering_discount", po::value<prob_t>( &priors.objectClustering.discount )->default_value( priors.objectClustering.discount ),
          "'discount' parameter for Pitman-Yor prior of objects clustering" )

        ( "prior.probe_clustering_concentration", po::value<prob_t>( &priors.probeClustering.concentration )->default_value( priors.probeClustering.concentration ),
          "'concentration' parameter for Pitman-Yor prior of probes clustering" )
        ( "prior.probe_clustering_discount", po::value<prob_t>( &priors.probeClustering.discount )->default_value( priors.probeClustering.discount ),
          "'discount' parameter for Pitman-Yor prior of probes clustering" )

        ( "prior.on_block_probability", po::value<prob_t>( &priors.cellEnablementProb )->default_value( priors.cellEnablementProb ),
          "prior probability that block is in \"on\" state" )
        ( "prior.objects_clusters_off_probes_cluster_rate", po::value<prob_t>( &priors.objectsCluOffProbesCluRate )->default_value( priors.objectsCluOffProbesCluRate ),
          "prior expectation of the number of additional objects clusters related to probes cluster via probes baits.\nFailure rate of geometrical distribution" )
        ( "prior.probes_clusters_off_objects_cluster_rate", po::value<prob_t>( &priors.probesCluOffObjectsCluRate )->default_value( priors.probesCluOffObjectsCluRate ),
          "prior expectation of the number of additional probes clusters related to objects cluster via probes baits.\nFailure rate of geometrical distribution" )

        ( "prior.false_hits", po::value<prob_t>( &priors.noise.failures )->default_value( priors.noise.failures ),
          "prior counts of false protein identifications (in off-blocks)" )
        ( "prior.true_misses", po::value<prob_t>( &priors.noise.successes )->default_value( priors.noise.successes ),
          "prior counts of true protein \"non-identifications\" (in off-blocks)" )

        ( "prior.object_multiple_rate", po::value<prob_t>( &priors.objectMultipleRate )->default_value( priors.objectMultipleRate ),
          "expected number of additional copies of object in its cluster. Failure rate of geometrical distribution" )
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
          "p, the power sequence length (L) is raised to (signal is divided by L^p for signal normalization)" )
        ( "signal.signal_counts_shape", po::value<log_prob_t>( &signalParams.scShape )->default_value( signalParams.scShape ),
          "[0,1), shape parameter of Lagrangian Poisson TRUE-POSITIVE spectral counts distribution (0 = classic Poisson)" )
        ( "signal.zero_inflation_factor", po::value<log_prob_t>( &signalParams.zeroInflationFactor )->default_value( signalParams.zeroInflationFactor ),
          "PDF(factor*mode()) is added to PDF(0) (0 = disable). Experimental" )
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
        ( "gibbs.signal_rate", po::value<prob_t>( &gibbsParams.signalRate )->default_value( gibbsParams.signalRate ),
          "rate of on-block abundance sampling steps" )
        ( "gibbs.object_multiple_rate", po::value<prob_t>( &gibbsParams.objectMultipleRate )->default_value( gibbsParams.objectMultipleRate ),
          "rate of object multiple (stoichiometry) sampling steps" )
        ( "gibbs.block_flip_rate", po::value<prob_t>( &gibbsParams.blockFlipRate )->default_value( gibbsParams.blockFlipRate ),
          "rate of on/off switching of random matrix blocks" )
        ( "gibbs.object_membership_rate", po::value<prob_t>( &gibbsParams.objectMembershipRate )->default_value( gibbsParams.objectMembershipRate ),
          "rate of changing cluster of single object" )
        ( "gibbs.probe_membership_rate", po::value<prob_t>( &gibbsParams.probeMembershipRate )->default_value( gibbsParams.probeMembershipRate ),
          "rate of changing cluster of single probe" )
        ( "gibbs.object_split_merge_rate", po::value<prob_t>( &gibbsParams.objectsSplitMergeRate )->default_value( gibbsParams.objectsSplitMergeRate ),
          "rate of objects clusters split/merge sampling steps" )
        ( "gibbs.probe_split_merge_rate", po::value<prob_t>( &gibbsParams.probesSplitMergeRate )->default_value( gibbsParams.probesSplitMergeRate ),
          "rate of objects clusters split/merge sampling steps" )
        ( "gibbs.block_resamples", po::value<size_t>( &gibbsParams.blockResamples )->default_value( gibbsParams.blockResamples ),
          "number of mandatory block abundance sampling steps after it was switched to on-state" )

        // object clustering steps parameters
        ( "gibbs.object_membership_samples", po::value<size_t>( &gibbsParams.objectClusterParams.existingClusterSamples )->default_value( gibbsParams.objectClusterParams.existingClusterSamples ),
          "number of sub-samples for block parameters after the assignment of object to an existing cluster" )
        ( "gibbs.object_new_cluster_samples", po::value<size_t>( &gibbsParams.objectClusterParams.newClusterSamples )->default_value( gibbsParams.objectClusterParams.newClusterSamples ),
          "number of sub-samples for block parameters after the assignment of object to a new cluster (singleton)" )
        ( "gibbs.object_cluster_split_samples", po::value<size_t>( &gibbsParams.objectsSplitMergeParams.splitLaunchSamples )->default_value( gibbsParams.objectsSplitMergeParams.splitLaunchSamples ),
          "number of sub-samples for block parameters after the split of object cluster" )
        ( "gibbs.object_cluster_merge_samples", po::value<size_t>( &gibbsParams.objectsSplitMergeParams.mergeLaunchSamples )->default_value( gibbsParams.objectsSplitMergeParams.mergeLaunchSamples ),
          "number of sub-samples for block parameters after the merge of object clusters" )

        // probe clustering steps parameters
        ( "gibbs.probe_membership_samples", po::value<size_t>( &gibbsParams.probeClusterParams.existingClusterSamples )->default_value( gibbsParams.probeClusterParams.existingClusterSamples ),
          "number of sub-samples for block parameters after the assignment of probe to an existing cluster" )
        ( "gibbs.probe_new_cluster_samples", po::value<size_t>( &gibbsParams.probeClusterParams.newClusterSamples )->default_value( gibbsParams.probeClusterParams.newClusterSamples ),
          "number of sub-samples for block parameters after the assignment of probe to a new cluster (singleton)" )
        ( "gibbs.probe_cluster_split_samples", po::value<size_t>( &gibbsParams.probesSplitMergeParams.splitLaunchSamples )->default_value( gibbsParams.probesSplitMergeParams.splitLaunchSamples ),
          "number of sub-samples for block parameters after the split of a probe cluster" )
        ( "gibbs.probe_cluster_merge_samples", po::value<size_t>( &gibbsParams.probesSplitMergeParams.mergeLaunchSamples )->default_value( gibbsParams.probesSplitMergeParams.mergeLaunchSamples ),
          "number of sub-samples for block parameters after the merge of probe clusters" )

        ( "gibbs.mean_object_rank", po::value<size_t>( &gibbsParams.meanObjectRank )->default_value( gibbsParams.meanObjectRank ),
          "rank (as defined by Pearson correlation) of an object considered as a candidate for split/merge and membership moves" )
        ( "gibbs.mean_probe_rank", po::value<size_t>( &gibbsParams.meanProbeRank )->default_value( gibbsParams.meanProbeRank ),
          "rank (as defined by Pearson correlation) of a probe considered as a candidate for split/merge and membership moves" )

        ( "gibbs.priors_update_period", po::value<size_t>( &gibbsParams.priorsUpdatePeriod )->default_value( gibbsParams.priorsUpdatePeriod ),
            "period (# of cclustering sampling iterations) of priors updating sampling step" )
        ;

    po::options_description log_params_desc( "Logging parameters" );
    log_params_desc.add_options()
        ( "log.samples_reporting_period", po::value<prob_t>( &gibbsParams.signalRate )->default_value( gibbsParams.signalRate ),
          "period (# of samples collected) to output info about collected samples" )
        ( "log.samples_max_reporting_delay", po::value<prob_t>( &gibbsParams.blockFlipRate )->default_value( gibbsParams.blockFlipRate ),
          "maximal time without reporting about collected samples (in case it takes too long)" )
        ;

    po::options_description output_params_desc( "Output parameters" );
    output_params_desc.add_options()
        ( "output.min_biclustering_ref_counts", po::value<size_t>( &ioParams.minCrossClusRefCount )->default_value( ioParams.minCrossClusRefCount ),
          "minimal number of times chessboard biclustering should appear in the random walk. Less frequent biclustering are removed from the output" )
        ( "output.min_objects_clustering_ref_counts", po::value<size_t>( &ioParams.minObjectsPtnRefCount )->default_value( ioParams.minObjectsPtnRefCount ),
          "minimal number of times objects clustering should appear in the random walk. Less frequent are removed" )
        ( "output.min_probes_clustering_ref_counts", po::value<size_t>( &ioParams.minProbesPtnRefCount )->default_value( ioParams.minProbesPtnRefCount ),
          "minimal number of times probes clustering should appear in the random walk. Less frequent are removed" )
        ;

    po::options_description all_cmdline_options_desc = basic_cmdline_options_desc;
    all_cmdline_options_desc
        .add( eesampler_options_desc )
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
        .add( eesampler_options_desc )
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

    // check output path existence
    if ( opt_map.count("output_file") == 0 ) {
        THROW_EXCEPTION( std::invalid_argument, "Output file not specified" );
    }
    boost::filesystem::path out_path( opt_map[ "output_file" ].as<std::string>() );
    if ( !out_path.parent_path().empty()
     && !boost::filesystem::exists( out_path.parent_path() )
    ){
        THROW_EXCEPTION( std::invalid_argument, "Output folder doesn't exist: " << out_path.parent_path() );
    }
    if ( !out_path.parent_path().empty()
     && !boost::filesystem::is_directory( out_path.parent_path() )
    ){
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

    if ( opt_map.count( "eesampler.simple_samplers_per_level" ) > 0 ) {
        std::stringstream ss( opt_map["eesampler.simple_samplers_per_level"].as<std::string>() );
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

