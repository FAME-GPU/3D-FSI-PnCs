function opts = set_local_defaults(opts,nev)
if ~isfield(opts,'efficiencyStats') || isempty(opts.efficiencyStats), opts.efficiencyStats = true; end
if ~isfield(opts,'trackMatvec') || isempty(opts.trackMatvec), opts.trackMatvec = opts.efficiencyStats; end
if ~isfield(opts,'saveHistory') || isempty(opts.saveHistory), opts.saveHistory = true; end
if ~isfield(opts,'keepLinearLog') || isempty(opts.keepLinearLog), opts.keepLinearLog = false; end
if ~isfield(opts,'historyCount') || isempty(opts.historyCount), opts.historyCount = max(10,nev); end
if ~isfield(opts,'referenceLambda'), opts.referenceLambda = []; end
if ~isfield(opts,'cacheBasisActions') || isempty(opts.cacheBasisActions), opts.cacheBasisActions = true; end
if ~isfield(opts,'linear') || isempty(opts.linear), opts.linear = struct(); end
if ~isfield(opts.linear,'allowFallback') || isempty(opts.linear.allowFallback), opts.linear.allowFallback = false; end

if ~isfield(opts,'gamma') || isempty(opts.gamma), opts.gamma = struct(); end
if ~isfield(opts.gamma,'enabled') || isempty(opts.gamma.enabled), opts.gamma.enabled = false; end
if ~isfield(opts.gamma,'lockSolid') || isempty(opts.gamma.lockSolid), opts.gamma.lockSolid = false; end
if ~isfield(opts.gamma,'lockFluid') || isempty(opts.gamma.lockFluid), opts.gamma.lockFluid = false; end
if ~isfield(opts.gamma,'eliminateZero') || isempty(opts.gamma.eliminateZero), opts.gamma.eliminateZero = false; end
if ~isfield(opts.gamma,'useLift') || isempty(opts.gamma.useLift), opts.gamma.useLift = false; end
if ~isfield(opts.gamma,'nullTol') || isempty(opts.gamma.nullTol), opts.gamma.nullTol = 1e-8; end
if ~isfield(opts.gamma,'couplingTol') || isempty(opts.gamma.couplingTol), opts.gamma.couplingTol = 1e-8; end
if ~isfield(opts.gamma,'solidDofLayout') || isempty(opts.gamma.solidDofLayout), opts.gamma.solidDofLayout = 'auto'; end
if ~isfield(opts.gamma,'liftTau_s') || isempty(opts.gamma.liftTau_s), opts.gamma.liftTau_s = 1; end
if ~isfield(opts.gamma,'liftTau_f') || isempty(opts.gamma.liftTau_f), opts.gamma.liftTau_f = 1; end
if ~isfield(opts.gamma,'keepNullVectors') || isempty(opts.gamma.keepNullVectors), opts.gamma.keepNullVectors = false; end
end