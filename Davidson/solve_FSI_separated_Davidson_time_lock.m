function [lambda, omega, U_s, U_f, P, info] = solve_FSI_separated_Davidson_time_lock(K_s, K_f, M_s, M_f, C, nev, target_omega, opts)

if nargin < 7 || isempty(target_omega)
    target_omega = 0;
end
if nargin < 8
    opts = struct();
end

opts = set_default_opts(opts,nev);
opts = set_local_defaults(opts,nev);
opts = set_lock_defaults_lock(opts);
solver_timer = tic;

detailed_timing = isfield(opts,'detailedTiming') && ~isempty(opts.detailedTiming) && logical(opts.detailedTiming);

track_prec_apply = true;
if isfield(opts,'trackPreconditionerApply') && ~isempty(opts.trackPreconditionerApply)
    track_prec_apply = logical(opts.trackPreconditionerApply);
end

timing_reduced_qep = 0;
timing_candidate_evaluation = 0;
timing_m_orthogonalization = 0;  
timing_append_basis_cache = 0;
timing_build_basis_cache = 0;

prec_apply_solid_time = 0;
prec_apply_fluid_time = 0;
prec_apply_solid_calls = 0;
prec_apply_fluid_calls = 0;
prec_s_original = [];
prec_f_original = [];

if track_prec_apply
    if isfield(opts,'prec_s') && ~isempty(opts.prec_s)
        prec_s_original = opts.prec_s;
        opts.prec_s = @timed_prec_s;
    end
    if isfield(opts,'prec_f') && ~isempty(opts.prec_f)
        prec_f_original = opts.prec_f;
        opts.prec_f = @timed_prec_f;
    end
end

if strcmpi(opts.target_type,'lambda')
    target_theta = sqrt(max(real(target_omega),0));
else
    target_theta = target_omega;
end

ns = size(K_s,1);
nf = size(K_f,1);
if size(C,1) ~= nf || size(C,2) ~= ns
    error('C must be n_f-by-n_s, consistent with C*u_s and C''*u_f.');
end

if ~isempty(opts.randomSeed)
    rng(opts.randomSeed);
end

is_complex_problem = ~isreal(K_s) || ~isreal(K_f) || ~isreal(M_s) || ...
                     ~isreal(M_f) || ~isreal(C);

is_gamma = isfield(opts,'gamma') && opts.gamma.enabled; 
Zs_null = [];
Zf_null = [];
ns0 = 0;
nf0 = 0;
gamma_info = struct('enabled',false);
if isfield(opts,'gamma') && isfield(opts.gamma,'isExactGamma')
    gamma_info.is_exact_gamma = logical(opts.gamma.isExactGamma);
else
    gamma_info.is_exact_gamma = false;
end

if is_gamma
    [Zs_null,Zf_null,gamma_info] = build_gamma_nullspaces(K_s,K_f,M_s,M_f,C,opts);
    ns0 = size(Zs_null,2);
    nf0 = size(Zf_null,2);

    opts.gamma.Zs = Zs_null;
    opts.gamma.Zf = Zf_null;
    opts.gamma.nullDimS = ns0;
    opts.gamma.nullDimF = nf0;

    nzero = ns0+nf0;
    gamma_info.zero_omega = zeros(nzero,1);
    gamma_info.zero_lambda = zeros(nzero,1);
    gamma_info.elimination_used = false;
    gamma_info.C00_relative_reduced_last = NaN;
    gamma_info.output_convention = ['nev counts TOTAL eigenpairs; Davidson computes ', ...
                                    'nev-nZeroLocked positive Ritz pairs and prepends ', ...
                                    'the selected exact Gamma zero eigenpairs.'];
    gamma_info.zero_mode_labels = [repmat({'solid_translation'},ns0,1); ...
                                   repmat({'fluid_constant'},nf0,1)];
    if opts.gamma.keepNullVectors
        gamma_info.Zs = Zs_null;
        gamma_info.Zf = Zf_null;
    end
else
    nzero = 0;
end

nev_positive = nev-nzero;
if nev_positive < 0
    error('nev=%d is smaller than the selected Gamma zero-mode count %d.',nev,nzero);
elseif nev_positive == 0
    error(['This implementation requires at least one positive eigenpair in addition to ', ...
           'the selected Gamma zero modes. Increase nev above %d.'],nzero);
end
if is_gamma
    gamma_info.total_requested = nev;
    gamma_info.positive_requested = nev_positive;
end

cap_s = opts.maxSubspace+ns0;
cap_f = opts.maxSubspace+nf0;

initial_space_timer = tic;
Zs0 = opts.Vs0;
Zf0 = opts.Vf0;

if is_gamma
    
    if ns0 > 0
        [V_s,locked_added_s] = append_m_orth_Scale_Invariant([],Zs_null,M_s,opts);
    else
        V_s = zeros(ns,0,'like',Zs0);
        locked_added_s = 0;
    end
    if nf0 > 0
        [V_f,locked_added_f] = append_m_orth_Scale_Invariant([],Zf_null,M_f,opts);
    else
        V_f = zeros(nf,0,'like',Zf0);
        locked_added_f = 0;
    end
    if locked_added_s ~= ns0 || locked_added_f ~= nf0
        error('Failed to insert the selected Gamma nullspaces into the trial spaces.');
    end
    [V_s,added_s] = append_m_orth_Scale_Invariant(V_s,Zs0,M_s,opts);
    [V_f,added_f] = append_m_orth_Scale_Invariant(V_f,Zf0,M_f,opts);

    if added_s == 0
        [V_s,added_s] = add_random_until_working(V_s,M_s,ns0,ns,is_complex_problem,opts,cap_s);
    end
    if added_f == 0
        [V_f,added_f] = add_random_until_working(V_f,M_f,nf0,nf,is_complex_problem,opts,cap_f);
    end
else
    [V_s,added_s] = append_m_orth_Scale_Invariant([],Zs0,M_s,opts);
    [V_f,added_f] = append_m_orth_Scale_Invariant([],Zf0,M_f,opts);
end

initial_subspace_time = toc(initial_space_timer);
if added_s == 0 || added_f == 0
    error('Initial working-subspace construction failed. Check M_s/M_f and initial vectors.');
end

cache = [];
if opts.cacheBasisActions
    if detailed_timing, tdetail = tic; end
    cache = build_basis_cache(K_s,K_f,M_s,M_f,C,V_s,V_f);
    if detailed_timing, timing_build_basis_cache = timing_build_basis_cache + toc(tdetail); end
end

if is_gamma
    if opts.cacheBasisActions
        opts.gamma.MZs = cache.MsVs(:,1:ns0);
        opts.gamma.MZf = cache.MfVf(:,1:nf0);
    else
        if ns0 > 0, opts.gamma.MZs = M_s*Zs_null; else, opts.gamma.MZs = zeros(ns,0,'like',Zs_null); end
        if nf0 > 0, opts.gamma.MZf = M_f*Zf_null; else, opts.gamma.MZf = zeros(nf,0,'like',Zf_null); end
    end
end

info = struct();
info.converged = false;
info.success = false;
info.message = '';
info.restart_count = 0;
info.outer_iterations = 0;
info.initial_subspace_time = initial_subspace_time;
info.history = empty_history();
info.linear = struct('domain',{},'theta',{},'rhsnorm',{}, ...
                     'flag',{},'relres',{},'iter',{},'solver',{}, ...
                     'time',{},'matvec',{},'matvec_time',{}, ...
                     'fallback_used',{},'exception',{});
info.gamma = gamma_info;
info.ew_count_time = zeros(nev,1);
info.for_itmax_time = 0;

info.timing = struct('enabled',detailed_timing, ...
                     'reduced_qep',0, ...
                     'candidate_evaluation',0, ...
                     'm_orthogonalization',0, ...
                     'append_basis_cache',0, ...
                     'build_basis_cache',0);

locked = init_positive_locks_lock(ns,nf,is_complex_problem);
info.locking = struct('enabled',opts.locking.enabled, ...
                      'thetaTol',opts.locking.thetaTol, ...
                      'vectorTol',opts.locking.vectorTol, ...
                      'count',0,'theta_lock_order',[], ...
                      'residual_lock_order',[],'lock_iteration',[], ...
                      'final_output_locked_mask',[], ...
                      'final_output_lock_iteration',[]);
last_active_cands = struct([]);

t_for_itmax = tic;
t_ew_count_time = tic;
for it = 1:opts.maxit
    info.outer_iterations = it;

    if detailed_timing, tdetail = tic; end
    red = solve_reduced_qep(K_s,K_f,M_s,M_f,C,V_s,V_f,target_theta,opts,cache);
    if detailed_timing, timing_reduced_qep = timing_reduced_qep + toc(tdetail); end
    if is_gamma
        info.gamma.elimination_used = info.gamma.elimination_used || red.gamma_elimination_used;
        info.gamma.C00_relative_reduced_last = red.C00_relative;
    end
    if isempty(red.theta)
        error(['The reduced QEP produced no admissible positive real Ritz values. ', ...
               'At Gamma, check that the initial working spaces contain non-null vectors.']);
    end

    maxCheck = min(numel(red.theta),max([nev_positive+opts.extraCandidates, ...
                                        opts.restartDim,opts.expandCandidates, ...
                                        opts.historyCount]));
    if detailed_timing, tdetail = tic; end
    cands = evaluate_candidates(red,K_s,K_f,M_s,M_f,C,V_s,V_f,maxCheck,cache);
    if detailed_timing, timing_candidate_evaluation = timing_candidate_evaluation + toc(tdetail); end

    [is_locked_replica,lock_id,~] = match_locked_candidates_lock(cands,locked,opts);

    if opts.locking.enabled && ~isempty(locked.theta)
        nseen = numel(unique(lock_id(lock_id>0)));
        if nseen < numel(locked.theta) && opts.verbose
            fprintf('[SepDav] warning: matched %d/%d stored positive locks in current Ritz set.\n', ...
                    nseen,numel(locked.theta));
        end
    end

    active_idx = find(~is_locked_replica);
    nlocked_before = numel(locked.theta);
    nneed_before = max(nev_positive-nlocked_before,0);
    target_active_idx = active_idx(1:min(nneed_before,numel(active_idx)));

    if opts.locking.enabled && ~isempty(target_active_idx)
        target_res = [cands(target_active_idx).res];
        new_lock_idx = target_active_idx(target_res <= opts.tol);
    else
        new_lock_idx = [];
    end

    if ~isempty(new_lock_idx)
        [locked,new_lock_count] = append_new_locks_lock(locked,cands,new_lock_idx,M_s,M_f,it);
        elapsed_lock = toc(t_ew_count_time);
        for kk = nlocked_before+1:numel(locked.theta)
            if kk <= numel(info.ew_count_time)
                info.ew_count_time(kk) = elapsed_lock;
            end
        end
        if opts.verbose
            for kk = 1:new_lock_count
                lid = nlocked_before+kk;
                fprintf('[SepDav] LOCK +%d/%d at it=%d | theta=% .12e | lambda=% .12e | res=%.3e\n', ...
                        lid,nev_positive,it,locked.theta(lid),locked.theta(lid)^2,locked.res(lid));
            end
        end

        [is_locked_replica,lock_id,~] = match_locked_candidates_lock(cands,locked,opts);
        active_idx = find(~is_locked_replica);
    end

    conv_count = numel(locked.theta);
    nneed = max(nev_positive-conv_count,0);
    last_active_cands = cands(active_idx);

    if opts.saveHistory
        info.history(end+1) = make_history_entry(it,V_s,V_f,cands, ...
                                                 conv_count,nev_positive,opts);
    end

    if opts.verbose
        
        display_idx = [];
        display_pool = active_idx(1:min(nneed,numel(active_idx)));
        if ~isempty(display_pool)
            display_res = [cands(display_pool).res];
            kk_display = find(display_res > opts.tol,1,'first');
            if ~isempty(kk_display)
                display_idx = display_pool(kk_display);
            end
        end

        if isempty(display_idx)
            if is_gamma
                fprintf(['[SepDav-Gamma-lock] it=%3d | workdim=(%d,%d) | null=(%d,%d) | ', ...
                         'lambda= N/A | res= N/A | locked=%d/%d\n'], ...
                        it,size(V_s,2)-ns0,size(V_f,2)-nf0,ns0,nf0, ...
                        conv_count,nev_positive);
            else
                fprintf(['[SepDav] it=%3d | dim=(%d,%d) | ', ...
                         'lambda= N/A | res= N/A | locked=%d/%d\n'], ...
                        it,size(V_s,2),size(V_f,2),conv_count,nev_positive);
            end
        else
            if is_gamma
                fprintf(['[SepDav-Gamma-lock] it=%3d | workdim=(%d,%d) | null=(%d,%d) | ', ...
                         'lambda=% .10e | res=%8.2e | locked=%d/%d\n'], ...
                        it,size(V_s,2)-ns0,size(V_f,2)-nf0,ns0,nf0, ...
                        cands(display_idx).theta^2,cands(display_idx).res,conv_count,nev_positive);
            else
                fprintf(['[SepDav] it=%3d | dim=(%d,%d) | ', ...
                         'lambda=% .10e | res=%8.2e | locked=%d/%d\n'], ...
                        it,size(V_s,2),size(V_f,2), ...
                        cands(display_idx).theta^2,cands(display_idx).res,conv_count,nev_positive);
            end
        end
    end

    if conv_count >= nev_positive
        info.converged = true;
        info.message = 'Converged with persistent positive Ritz-pair locking.';
        break;
    end

    target_active_idx = active_idx(1:min(nneed,numel(active_idx)));
    if isempty(target_active_idx)
        bad = [];
    else
        target_bad = target_active_idx([cands(target_active_idx).res] > opts.tol);
        bad = target_bad(:).';
    end

    if numel(bad) < opts.expandCandidates
        extra_start = min(nneed,numel(active_idx))+1;
        if extra_start <= numel(active_idx)
            extra_pool = active_idx(extra_start:end);
            extra_bad = extra_pool([cands(extra_pool).res] > opts.tol);
            take = min(opts.expandCandidates-numel(bad),numel(extra_bad));
            bad = [bad, extra_bad(1:take)]; 
        end
    end
    if isempty(bad) && ~isempty(active_idx)
        bad = active_idx(1:min(opts.expandCandidates,numel(active_idx)));
    end
    bad = bad(1:min(opts.expandCandidates,numel(bad)));

    nnew_per_side = numel(bad)*(1+double(opts.addCross));
    need_restart = (size(V_s,2)+nnew_per_side > cap_s) || ...
                   (size(V_f,2)+nnew_per_side > cap_f);

    if need_restart
        
        nlock_pos = numel(locked.theta);
        active_keep_count = max(opts.restartDim-nlock_pos,0);
        active_keep_idx = active_idx(1:min(active_keep_count,numel(active_idx)));
        restart_active = cands(active_keep_idx);

        if nlock_pos > opts.restartDim && opts.verbose
            fprintf(['[SepDav] note: positive locks (%d) exceed restartDim=%d; ', ...
                     'all locks are retained and no active Ritz vector is kept by quota.\n'], ...
                    nlock_pos,opts.restartDim);
        end

        if is_gamma
            [V_s,V_f] = restart_from_candidates_gamma_lock(locked,restart_active,Zs_null,Zf_null, ...
                                                                           M_s,M_f,opts,ns,nf, ...
                                                                           is_complex_problem,cap_s,cap_f);
        else
            [V_s,V_f] = restart_from_candidates_lock(locked,restart_active,M_s,M_f,opts, ...
                                                                    ns,nf,is_complex_problem);
        end

        if size(V_s,2)+nnew_per_side > cap_s || size(V_f,2)+nnew_per_side > cap_f
            error(['Persistent locks leave insufficient subspace room for the planned correction block. ', ...
                   'After restart dim=(%d,%d), planned new=%d/side, cap=(%d,%d). ', ...
                   'Increase maxSubspace or reduce expandCandidates/restartDim.'], ...
                  size(V_s,2),size(V_f,2),nnew_per_side,cap_s,cap_f);
        end

        if opts.cacheBasisActions
            if detailed_timing, tdetail = tic; end
            cache = build_basis_cache(K_s,K_f,M_s,M_f,C,V_s,V_f);
            if detailed_timing, timing_build_basis_cache = timing_build_basis_cache + toc(tdetail); end
        end

        info.restart_count = info.restart_count+1;
        if opts.verbose
            if is_gamma
                fprintf(['[SepDav-Gamma-lock] thick restart -> workdim=(%d,%d), null=(%d,%d), ', ...
                         'positive locks=%d, active kept=%d.\n'], ...
                        size(V_s,2)-ns0,size(V_f,2)-nf0,ns0,nf0,nlock_pos,numel(restart_active));
            else
                fprintf('[SepDav] thick restart -> dim=(%d,%d), positive locks=%d, active kept=%d.\n', ...
                        size(V_s,2),size(V_f,2),nlock_pos,numel(restart_active));
            end
        end
    end

    old_ms_for_cache = size(V_s,2);
    old_mf_for_cache = size(V_f,2);

    Zs = [];
    Zf = [];
    for kk = 1:numel(bad)
        j = bad(kk);
        th = cands(j).theta;

        [zs,st] = solve_shifted_system_time(K_s,M_s,th,cands(j).r_s,'solid',opts);
        info.linear(end+1) = st; 
        [zf,st] = solve_shifted_system_time(K_f,M_f,th,cands(j).r_f,'fluid',opts);
        info.linear(end+1) = st; 

        Zs = [Zs,zs]; 
        Zf = [Zf,zf]; 

        if opts.addCross
            rhs_s_cross = C' * cands(j).u_f;
            rhs_f_cross = C  * cands(j).u_s;
            [zs_c,st] = solve_shifted_system_time(K_s,M_s,th,rhs_s_cross,'solid',opts);
            info.linear(end+1) = st; 
            [zf_c,st] = solve_shifted_system_time(K_f,M_f,th,rhs_f_cross,'fluid',opts);
            info.linear(end+1) = st; 
            Zs = [Zs,zs_c]; 
            Zf = [Zf,zf_c]; 
        end
    end

    [V_s,added_s] = append_m_orth_Scale_Invariant(V_s,Zs,M_s,opts);
    [V_f,added_f] = append_m_orth_Scale_Invariant(V_f,Zf,M_f,opts);

    if added_s == 0 && size(V_s,2) < cap_s
        z = random_block(ns,1,is_complex_problem);
        [V_s,~] = append_m_orth_Scale_Invariant(V_s,z,M_s,opts);
    end
    if added_f == 0 && size(V_f,2) < cap_f
        z = random_block(nf,1,is_complex_problem);
        [V_f,~] = append_m_orth_Scale_Invariant(V_f,z,M_f,opts);
    end

    if opts.cacheBasisActions && ...
            (size(V_s,2) > old_ms_for_cache || size(V_f,2) > old_mf_for_cache)
        if detailed_timing, tdetail = tic; end
        cache = append_basis_cache(cache,K_s,K_f,M_s,M_f,C,V_s,V_f, ...
                                   old_ms_for_cache,old_mf_for_cache);
        if detailed_timing, timing_append_basis_cache = timing_append_basis_cache + toc(tdetail); end
    end

    if size(V_s,2) > cap_s || size(V_f,2) > cap_f
        error(['Subspace cap exceeded after restart/expansion: dim=(%d,%d), ', ...
               'allowed=(%d,%d).'],size(V_s,2),size(V_f,2),cap_s,cap_f);
    end

end
info.for_itmax_time = toc(t_for_itmax);

if isempty(last_active_cands) && numel(locked.theta) < nev_positive
    if detailed_timing, tdetail = tic; end
    red = solve_reduced_qep(K_s,K_f,M_s,M_f,C,V_s,V_f,target_theta,opts,cache);
    if detailed_timing, timing_reduced_qep = timing_reduced_qep + toc(tdetail); end
    if detailed_timing, tdetail = tic; end
    cands_last = evaluate_candidates(red,K_s,K_f,M_s,M_f,C,V_s,V_f, ...
                                     min(numel(red.theta),nev_positive+opts.extraCandidates),cache);
    if detailed_timing, timing_candidate_evaluation = timing_candidate_evaluation + toc(tdetail); end
    [is_locked_replica,~,~] = match_locked_candidates_lock(cands_last,locked,opts);
    last_active_cands = cands_last(~is_locked_replica);
end

[omega_pos,U_s_pos,U_f_pos,qep_res_pos,final_locked_mask,final_lock_iter] = ...
    assemble_final_positive_lock(locked,last_active_cands,nev_positive,target_theta,opts, ...
                                 ns,nf,is_complex_problem);
nkeep = numel(omega_pos);
lambda_pos = omega_pos.^2;
P_pos = bsxfun(@times,U_f_pos,reshape(omega_pos,1,[]));

if ~info.converged && isempty(info.message)
    if nkeep == 0
        info.message = sprintf('Stopped after maxit=%d with no admissible positive Ritz pairs.',opts.maxit);
    else
        info.message = sprintf('Stopped after maxit=%d with %d/%d positive pairs locked. Best returned QEP residual = %.3e.', ...
                               opts.maxit,numel(locked.theta),nev_positive,min(qep_res_pos));
    end
end

if nkeep > 0
    [gep_res_pos,Rmax] = compute_FSI_GEP_residual(K_s,K_f,M_s,M_f,C, ...
                                                  lambda_pos,U_s_pos,P_pos);
else
    gep_res_pos = [];
    Rmax = NaN;
end

if nzero > 0
    U_s_zero = [Zs_null, zeros(ns,nf0,'like',U_s_pos)];
    U_f_zero = [zeros(nf,ns0,'like',U_f_pos), Zf_null];
    omega_zero = zeros(nzero,1);
    lambda_zero = zeros(nzero,1);
    P_zero = zeros(nf,nzero,'like',U_f_pos);
else
    U_s_zero = zeros(ns,0,'like',U_s_pos);
    U_f_zero = zeros(nf,0,'like',U_f_pos);
    omega_zero = zeros(0,1);
    lambda_zero = zeros(0,1);
    P_zero = zeros(nf,0,'like',U_f_pos);
end

omega = [omega_zero;omega_pos];
lambda = [lambda_zero;lambda_pos];
U_s = [U_s_zero,U_s_pos];
U_f = [U_f_zero,U_f_pos];
P = [P_zero,P_pos];

nout = numel(omega);
if nout > 0
    Gm = U_s'*M_s*U_s+U_f'*M_f*U_f;
    info.mass_norm_error = max(abs(diag(Gm)-1));
    Gc = U_s'*C'*U_f+U_f'*C*U_s;
    OmegaSum = omega(:)+omega(:).';
    Gorth = Gc+OmegaSum.*Gm;
    Gorth_off = Gorth-diag(diag(Gorth));
    info.qep_orth_error = norm(Gorth_off,2);
else
    info.mass_norm_error = NaN;
    info.qep_orth_error = NaN;
end

info.final_QEP_residuals = qep_res_pos;
info.final_GEP_residuals = gep_res_pos;
info.final_R_max = Rmax;
info.final_omega = omega_pos;
info.final_lambda = lambda_pos;
info.final_QEP_residuals_positive = qep_res_pos;
info.final_GEP_residuals_positive = gep_res_pos;
info.final_QEP_residuals_with_zeros = [zeros(nzero,1);qep_res_pos];
info.final_GEP_residuals_with_zeros = [NaN(nzero,1);gep_res_pos(:)];
info.returned_omega = omega;
info.returned_lambda = lambda;
info.locking.count = numel(locked.theta);
info.locking.theta_lock_order = locked.theta;
info.locking.residual_lock_order = locked.res;
info.locking.lock_iteration = locked.iter;
info.locking.final_output_locked_mask = final_locked_mask;
info.locking.final_output_lock_iteration = final_lock_iter;
if is_gamma
    info.gamma.zero_indices = (1:nzero).';
    info.gamma.positive_indices = (nzero+(1:nkeep)).';
    info.gamma.returned_zero_mode_count = nzero;
end

info.timing.enabled = detailed_timing;
info.timing.reduced_qep = timing_reduced_qep;
info.timing.candidate_evaluation = timing_candidate_evaluation;
info.timing.m_orthogonalization = timing_m_orthogonalization;
info.timing.append_basis_cache = timing_append_basis_cache;
info.timing.build_basis_cache = timing_build_basis_cache;

info.efficiency = build_efficiency_summary(info.linear);
info.efficiency.preconditioner_apply.solid = struct('calls',prec_apply_solid_calls,'time',prec_apply_solid_time);
info.efficiency.preconditioner_apply.fluid = struct('calls',prec_apply_fluid_calls,'time',prec_apply_fluid_time);
info.efficiency.total_preconditioner_apply_time = prec_apply_solid_time + prec_apply_fluid_time;
info.efficiency.total_time = toc(solver_timer);
info.efficiency.initial_subspace_time = initial_subspace_time;
info.efficiency.timing = info.timing;

info.inner_failure_count = info.efficiency.inner_failure_count;
info.inner_maxit_count = info.efficiency.inner_maxit_count;
info.inner_hard_failure_count = info.efficiency.inner_hard_failure_count;
info.fallback_used = info.efficiency.fallback_used;
info.success = info.converged && info.inner_hard_failure_count == 0 && ~info.fallback_used;
info.opts = strip_large_opts(opts);

if opts.verbose
    if is_gamma
        fprintf('\nGamma nullspace diagnostics:\n');
        fprintf('  lock solid / fluid        = %d / %d\n',info.gamma.lock_solid,info.gamma.lock_fluid);
        fprintf('  locked zero modes         = %d (solid %d, fluid %d)\n', ...
                nzero,ns0,nf0);
        if ns0 > 0
            fprintf('  solid null residual       = %.3e\n',info.gamma.solid_null_residual);
        else
            fprintf('  solid null residual       = N/A (not locked)\n');
        end
        if nf0 > 0
            fprintf('  fluid null residual       = %.3e\n',info.gamma.fluid_null_residual);
        else
            fprintf('  fluid null residual       = N/A (not locked)\n');
        end
        if ns0 > 0 && nf0 > 0
            fprintf('  relative ||Zf''*C*Zs||    = %.3e\n',info.gamma.C00_relative);
        else
            fprintf('  relative ||Zf''*C*Zs||    = N/A (one-sided lock)\n');
        end
        fprintf('  reduced zero elimination  = %d\n',info.gamma.elimination_used);
        fprintf('  total / positive returned = %d / %d\n',numel(omega),nkeep);
        fprintf('  correction lift tau_s/f   = %.3e / %.3e\n',opts.gamma.liftTau_s,opts.gamma.liftTau_f);
    end
    fprintf('Mass normalization error = %.3e\n',info.mass_norm_error);
    fprintf('QEP orthogonality error   = %.3e\n',info.qep_orth_error);
    fprintf('Unified original-GEP Rmax = %.3e\n',info.final_R_max);
    if opts.efficiencyStats
        print_efficiency_summary(info.efficiency);
    end
end

if ~opts.keepLinearLog
    info.linear = [];
end

    function y = timed_prec_s(x)
        t0 = tic;
        y = prec_s_original(x);
        prec_apply_solid_time = prec_apply_solid_time + toc(t0);
        prec_apply_solid_calls = prec_apply_solid_calls + 1;
    end

    function y = timed_prec_f(x)
        t0 = tic;
        y = prec_f_original(x);
        prec_apply_fluid_time = prec_apply_fluid_time + toc(t0);
        prec_apply_fluid_calls = prec_apply_fluid_calls + 1;
    end
end

function opts = set_lock_defaults_lock(opts)

if ~isfield(opts,'locking') || ~isstruct(opts.locking)
    opts.locking = struct();
end
if ~isfield(opts.locking,'enabled') || isempty(opts.locking.enabled)
    opts.locking.enabled = true;
end
if ~isfield(opts.locking,'thetaTol') || isempty(opts.locking.thetaTol)
    opts.locking.thetaTol = 1e-6;
end
if ~isfield(opts.locking,'vectorTol') || isempty(opts.locking.vectorTol)
    opts.locking.vectorTol = 1e-3;
end
end

function locked = init_positive_locks_lock(ns,nf,is_complex_problem)
locked = struct();
locked.theta = zeros(0,1);
locked.res = zeros(0,1);
locked.iter = zeros(0,1);
locked.U_s = zeros(ns,0);
locked.U_f = zeros(nf,0);
locked.MsU_s = zeros(ns,0);
locked.MfU_f = zeros(nf,0);
if is_complex_problem
    locked.U_s = complex(locked.U_s);
    locked.U_f = complex(locked.U_f);
    locked.MsU_s = complex(locked.MsU_s);
    locked.MfU_f = complex(locked.MfU_f);
end
end

function [locked,nadded] = append_new_locks_lock(locked,cands,idx,M_s,M_f,it)
nadded = 0;
for kk = 1:numel(idx)
    c = cands(idx(kk));
    locked.theta(end+1,1) = c.theta;
    locked.res(end+1,1) = c.res;
    locked.iter(end+1,1) = it;
    locked.U_s(:,end+1) = c.u_s;
    locked.U_f(:,end+1) = c.u_f;
    
    locked.MsU_s(:,end+1) = M_s*c.u_s;
    locked.MfU_f(:,end+1) = M_f*c.u_f;
    nadded = nadded+1;
end
end

function [omega_pos,U_s_pos,U_f_pos,qep_res_pos,locked_mask,lock_iter] = ...
    assemble_final_positive_lock(locked,active_cands,nev_positive,target_theta,opts,ns,nf,is_complex_problem)

nlocked = numel(locked.theta);
nfill = min(max(nev_positive-nlocked,0),numel(active_cands));
ntotal = min(nev_positive,nlocked+nfill);

omega_pos = zeros(ntotal,1);
qep_res_pos = zeros(ntotal,1);
U_s_pos = zeros(ns,ntotal);
U_f_pos = zeros(nf,ntotal);
locked_mask = false(ntotal,1);
lock_iter = NaN(ntotal,1);
if is_complex_problem
    U_s_pos = complex(U_s_pos);
    U_f_pos = complex(U_f_pos);
end

if nlocked > 0
    omega_pos(1:nlocked) = locked.theta;
    qep_res_pos(1:nlocked) = locked.res;
    U_s_pos(:,1:nlocked) = locked.U_s;
    U_f_pos(:,1:nlocked) = locked.U_f;
    locked_mask(1:nlocked) = true;
    lock_iter(1:nlocked) = locked.iter;
end
for j = 1:nfill
    k = nlocked+j;
    omega_pos(k) = active_cands(j).theta;
    qep_res_pos(k) = active_cands(j).res;
    U_s_pos(:,k) = active_cands(j).u_s;
    U_f_pos(:,k) = active_cands(j).u_f;
end

ord = output_order_lock(omega_pos,target_theta,opts);
omega_pos = omega_pos(ord);
qep_res_pos = qep_res_pos(ord);
U_s_pos = U_s_pos(:,ord);
U_f_pos = U_f_pos(:,ord);
locked_mask = locked_mask(ord);
lock_iter = lock_iter(ord);
end

function ord = output_order_lock(theta,target_theta,opts)
switch lower(opts.select)
    case 'smallest'
        score = theta;
    case 'closest'
        score = abs(theta-target_theta);
    case 'greater_than'
        score = abs(theta-target_theta);
        score(theta <= target_theta) = inf;
    otherwise
        error('Unknown opts.select = %s.',opts.select);
end
[~,ord] = sort(score,'ascend');
end

function [V,added] = add_random_until_working(V,M,nlocked,n,is_complex_problem,opts,cap)
added = 0;
trials = 0;
while size(V,2) <= nlocked && size(V,2) < cap && trials < 10
    z = random_block(n,1,is_complex_problem);
    [V,a] = append_m_orth_Scale_Invariant(V,z,M,opts);
    added = added+a;
    trials = trials+1;
end
end

function [V_s,V_f] = restart_from_candidates_gamma_lock(locked,active_cands,Zs_null,Zf_null,M_s,M_f,opts, ...
                                                              ns,nf,is_complex_problem,cap_s,cap_f)

if ~isempty(Zs_null)
    [V_s,~] = append_m_orth_Scale_Invariant([],Zs_null,M_s,opts);
else
    V_s = zeros(ns,0);
    if is_complex_problem, V_s = complex(V_s); end
end
if ~isempty(Zf_null)
    [V_f,~] = append_m_orth_Scale_Invariant([],Zf_null,M_f,opts);
else
    V_f = zeros(nf,0);
    if is_complex_problem, V_f = complex(V_f); end
end

if ~isempty(locked.theta)
    [V_s,~] = append_m_orth_Scale_Invariant(V_s,locked.U_s,M_s,opts);
    [V_f,~] = append_m_orth_Scale_Invariant(V_f,locked.U_f,M_f,opts);
end
if ~isempty(active_cands)
    Us = [active_cands.u_s];
    Uf = [active_cands.u_f];
    [V_s,~] = append_m_orth_Scale_Invariant(V_s,Us,M_s,opts);
    [V_f,~] = append_m_orth_Scale_Invariant(V_f,Uf,M_f,opts);
end

target_work_s = min([opts.initDim,opts.restartDim,cap_s-size(Zs_null,2)]);
target_work_f = min([opts.initDim,opts.restartDim,cap_f-size(Zf_null,2)]);

trials = 0;
while size(V_s,2)-size(Zs_null,2) < target_work_s && size(V_s,2) < cap_s && trials < 50
    z = random_block(ns,1,is_complex_problem);
    [V_s,~] = append_m_orth_Scale_Invariant(V_s,z,M_s,opts);
    trials = trials+1;
end

trials = 0;
while size(V_f,2)-size(Zf_null,2) < target_work_f && size(V_f,2) < cap_f && trials < 50
    z = random_block(nf,1,is_complex_problem);
    [V_f,~] = append_m_orth_Scale_Invariant(V_f,z,M_f,opts);
    trials = trials+1;
end
end

function Z = random_block(n,m,is_complex_problem)
Z = randn(n,m);
if is_complex_problem
    Z = Z+1i*randn(n,m);
end
end

function cache = build_basis_cache(K_s,K_f,M_s,M_f,C,V_s,V_f)

cache = struct();
cache.KsVs = K_s*V_s;
cache.MsVs = M_s*V_s;
cache.CVs  = C*V_s;
cache.KfVf = K_f*V_f;
cache.MfVf = M_f*V_f;
cache.CtVf = C'*V_f;
end

function cache = append_basis_cache(cache,K_s,K_f,M_s,M_f,C,V_s,V_f,old_ms,old_mf)

new_ms = size(V_s,2);
new_mf = size(V_f,2);

if isempty(cache) || ~isstruct(cache) || ...
        ~isfield(cache,'KsVs') || size(cache.KsVs,2) ~= old_ms || ...
        ~isfield(cache,'KfVf') || size(cache.KfVf,2) ~= old_mf
    cache = build_basis_cache(K_s,K_f,M_s,M_f,C,V_s,V_f);
    return;
end

if new_ms > old_ms
    dVs = V_s(:,old_ms+1:new_ms);
    cache.KsVs = [cache.KsVs, K_s*dVs];
    cache.MsVs = [cache.MsVs, M_s*dVs];
    cache.CVs  = [cache.CVs,  C*dVs];
end
if new_mf > old_mf
    dVf = V_f(:,old_mf+1:new_mf);
    cache.KfVf = [cache.KfVf, K_f*dVf];
    cache.MfVf = [cache.MfVf, M_f*dVf];
    cache.CtVf = [cache.CtVf, C'*dVf];
end
end

function h = empty_history()
h = struct('iter',{},'dim_s',{},'dim_f',{}, ...
           'ritz_theta',{},'ritz_lambda',{},'ritz_residuals',{}, ...
           'num_converged',{},'max_rel_eig_error',{}, ...
           'max_scaled_residual',{});
end

function h = make_history_entry(it,V_s,V_f,cands,conv_count,nev,opts)
nstore = min(opts.historyCount,numel(cands));
theta = NaN(opts.historyCount,1);
lambda = NaN(opts.historyCount,1);
res = NaN(opts.historyCount,1);
if nstore > 0
    theta(1:nstore) = [cands(1:nstore).theta].';
    lambda(1:nstore) = theta(1:nstore).^2;
    res(1:nstore) = [cands(1:nstore).res].';
end

h.iter = it;
h.dim_s = size(V_s,2);
h.dim_f = size(V_f,2);
h.ritz_theta = theta;
h.ritz_lambda = lambda;
h.ritz_residuals = res;
h.num_converged = conv_count;

if numel(opts.referenceLambda) >= nev && nstore >= nev
    ref = opts.referenceLambda(1:nev);
    cur = lambda(1:nev);
    denom = max(abs(ref(:)),eps);
    h.max_rel_eig_error = max(abs(cur-ref(:))./denom);
else
    h.max_rel_eig_error = NaN;
end

if nstore > 0
    h.max_scaled_residual = max(res(1:min(nev,nstore)),[],'omitnan');
else
    h.max_scaled_residual = NaN;
end
end

function eff = build_efficiency_summary(linear)
eff = struct();
eff.linear.solid = summarize_linear(linear,'solid');
eff.linear.fluid = summarize_linear(linear,'fluid');
eff.matvec.solid = summarize_matvec(linear,'solid');
eff.matvec.fluid = summarize_matvec(linear,'fluid');
eff.total_shifted_solve_calls = eff.linear.solid.calls+eff.linear.fluid.calls;
eff.total_krylov_iterations = eff.linear.solid.iter_total+eff.linear.fluid.iter_total;
eff.total_linear_solve_time = sum_finite([eff.linear.solid.time,eff.linear.fluid.time]);
eff.total_operator_applications = sum_finite([eff.matvec.solid.count,eff.matvec.fluid.count]);
eff.total_matvec_time = sum_finite([eff.matvec.solid.time,eff.matvec.fluid.time]);

if isempty(linear)
    eff.inner_failure_count = 0;
    eff.inner_maxit_count = 0;
    eff.inner_hard_failure_count = 0;
    eff.fallback_used = false;
else
    flags = [linear.flag].';
    finite_flags = flags(isfinite(flags));
    
    eff.inner_failure_count = sum(finite_flags ~= 0);
    eff.inner_maxit_count = sum(finite_flags == 1);
    eff.inner_hard_failure_count = sum(finite_flags ~= 0 & finite_flags ~= 1);
    eff.fallback_used = any([linear.fallback_used]);
end
end

function s = summarize_linear(linear,domain)
idx = select_domain(linear,domain);
if ~any(idx)
    s = struct('calls',0,'iter_total',0,'time',0,'failed_calls',0, ...
               'maxit_calls',0,'hard_failed_calls',0);
    return;
end
L = linear(idx);
it = [L.iter].';
tm = [L.time].';
flags = [L.flag].';
s.calls = numel(L);
s.iter_total = sum(it(isfinite(it)));
s.time = sum(tm(isfinite(tm)));
s.failed_calls = sum(isfinite(flags) & flags ~= 0);
s.maxit_calls = sum(isfinite(flags) & flags == 1);
s.hard_failed_calls = sum(isfinite(flags) & flags ~= 0 & flags ~= 1);
end

function s = summarize_matvec(linear,domain)
idx = select_domain(linear,domain);
if ~any(idx)
    s = struct('count',0,'time',0);
    return;
end
L = linear(idx);
mv = [L.matvec].';
mt = [L.matvec_time].';
if any(isfinite(mv)), s.count = sum(mv(isfinite(mv))); else, s.count = NaN; end
if any(isfinite(mt)), s.time = sum(mt(isfinite(mt))); else, s.time = NaN; end
end

function idx = select_domain(linear,domain)
if isempty(linear)
    idx = false(0,1);
else
    idx = strcmpi({linear.domain}.',domain);
end
end

function value = sum_finite(x)
x = x(isfinite(x));
if isempty(x), value = NaN; else, value = sum(x); end
end

function print_efficiency_summary(eff)
fprintf('\n============= SepDav V2 efficiency =================\n');
fprintf('Davidson solver wall time    : %.6f s\n',eff.total_time);
fprintf('Solid shifted solves         : %d calls, %g iterations, %.6f s\n', ...
        eff.linear.solid.calls,eff.linear.solid.iter_total,eff.linear.solid.time);
fprintf('Fluid shifted solves         : %d calls, %g iterations, %.6f s\n', ...
        eff.linear.fluid.calls,eff.linear.fluid.iter_total,eff.linear.fluid.time);
fprintf('Total operator applications  : %g\n',eff.total_operator_applications);
fprintf('Total Krylov matvec time     : %.6f s\n',eff.total_matvec_time);
fprintf('Inner nonconv/maxit/hard     : %d / %d / %d\n', ...
        eff.inner_failure_count,eff.inner_maxit_count,eff.inner_hard_failure_count);
fprintf('Fallback used                : %d\n',eff.fallback_used);
if isfield(eff,'preconditioner_apply')
    fprintf('Solid prec apply             : %d calls, %.6f s\n', ...
            eff.preconditioner_apply.solid.calls,eff.preconditioner_apply.solid.time);
    fprintf('Fluid prec apply             : %d calls, %.6f s\n', ...
            eff.preconditioner_apply.fluid.calls,eff.preconditioner_apply.fluid.time);
end
if isfield(eff,'timing')
    fprintf('Reduced QEP time             : %.6f s\n',eff.timing.reduced_qep);
    fprintf('Candidate evaluation time    : %.6f s\n',eff.timing.candidate_evaluation);
    fprintf('M-orthogonalization time     : N/A (not measured)\n');
    fprintf('Append basis cache time      : %.6f s\n',eff.timing.append_basis_cache);
    fprintf('Build basis cache time       : %.6f s\n',eff.timing.build_basis_cache);
end
fprintf('=====================================================\n\n');
end

function small = strip_large_opts(opts)
small = opts;
large_fields = {'Vs0','Vf0','prec_s','prec_f','shiftSolve_s','shiftSolve_f'};
for k = 1:numel(large_fields)
    if isfield(small,large_fields{k})
        small = rmfield(small,large_fields{k});
    end
end
if isfield(small,'gamma')
    gamma_large = {'Zs','Zf','MZs','MZf'};
    for k = 1:numel(gamma_large)
        if isfield(small.gamma,gamma_large{k})
            small.gamma = rmfield(small.gamma,gamma_large{k});
        end
    end
end
end
