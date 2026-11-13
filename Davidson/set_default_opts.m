function opts = set_default_opts(opts, nev)
    opts = setfield_if_missing(opts, 'select', 'smallest');
    opts = setfield_if_missing(opts, 'target_type', 'omega');
    opts = setfield_if_missing(opts, 'tol', 1e-8);
    opts = setfield_if_missing(opts, 'maxit', 80);
    opts = setfield_if_missing(opts, 'initDim', max(nev+3, 6));
    opts = setfield_if_missing(opts, 'maxSubspace', 50);
    opts = setfield_if_missing(opts, 'restartDim', max(nev+6, 12));
    opts = setfield_if_missing(opts, 'expandCandidates', min(2, max(1, nev)));
    opts = setfield_if_missing(opts, 'extraCandidates', max(4, nev));
    opts = setfield_if_missing(opts, 'addCross', true);
    opts = setfield_if_missing(opts, 'verbose', true);
    opts = setfield_if_missing(opts, 'realTol', 1e-8);
    opts = setfield_if_missing(opts, 'dupTol', 1e-7);
    opts = setfield_if_missing(opts, 'orthTol', 1e-12);
    opts = setfield_if_missing(opts, 'maxReorth', 5);
    opts = setfield_if_missing(opts, 'minOmega', 1e-12);
    opts = setfield_if_missing(opts, 'randomSeed', []);

    if ~isfield(opts, 'linear') || isempty(opts.linear)
        opts.linear = struct();
    end
    opts.linear = setfield_if_missing(opts.linear, 'solver', 'minres');
    opts.linear = setfield_if_missing(opts.linear, 'tol', 1e-2);
    opts.linear = setfield_if_missing(opts.linear, 'maxit', 300);
    opts.linear = setfield_if_missing(opts.linear, 'directThreshold', 1000);
    opts.linear = setfield_if_missing(opts.linear, 'useDirect', false);
    opts.linear = setfield_if_missing(opts.linear, 'useIC', false);
    opts.linear = setfield_if_missing(opts.linear, 'icholDroptol', 1e-3);
end