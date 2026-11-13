function DOpts = davidson_opts_from_Opt(Opt)

DOpts = struct();

DOpts.select = 'smallest';
DOpts.target_type = 'omega';

DOpts.tol = Opt.eig_tol;
DOpts.maxit = Opt.max_outer_iterations;

DOpts.initDim = Opt.subspace.init_dim;
DOpts.maxSubspace = Opt.subspace.max_dim;
DOpts.restartDim = Opt.subspace.restart_dim;
DOpts.expandCandidates = Opt.subspace.expand_pairs;
DOpts.extraCandidates = Opt.subspace.extra_candidates;
DOpts.addCross = Opt.add_cross;

DOpts.randomSeed = Opt.random_seed;
DOpts.verbose = true;

DOpts.gamma = Opt.gamma; 

DOpts.linear.solver = Opt.linear.solver;
DOpts.linear.tol = Opt.linear.tol;
DOpts.linear.maxit = Opt.linear.maxit;
DOpts.linear.useDirect = false;
DOpts.linear.useIC = false;
DOpts.linear.allowFallback = Opt.linear.allow_fallback;

DOpts.efficiencyStats = Opt.stats.enable;
DOpts.trackMatvec = Opt.stats.track_matvec;
DOpts.saveHistory = Opt.stats.save_history;
DOpts.keepLinearLog = Opt.stats.keep_inner_log;
DOpts.historyCount = Opt.stats.history_count;
DOpts.referenceLambda = Opt.reference_lambda;
DOpts.cacheBasisActions = Opt.cacheBasisActions; 
DOpts.detailedTiming = Opt.stats.detailed_timing;
DOpts.trackPreconditionerApply = true;  

end