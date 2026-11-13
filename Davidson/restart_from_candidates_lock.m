function [V_s_new, V_f_new] = restart_from_candidates_lock(locked, active_cands, M_s, M_f, opts, ns, nf, is_complex_problem)

    nlocked = numel(locked.theta);
    nactive = numel(active_cands);

    Zs = zeros(ns, nlocked+nactive);
    Zf = zeros(nf, nlocked+nactive);
    if is_complex_problem
        Zs = complex(Zs); Zf = complex(Zf);
    end

    if nlocked > 0
        Zs(:,1:nlocked) = locked.U_s;
        Zf(:,1:nlocked) = locked.U_f;
    end
    for j = 1:nactive
        Zs(:,nlocked+j) = active_cands(j).u_s;
        Zf(:,nlocked+j) = active_cands(j).u_f;
    end

    [V_s_new, ~] = append_m_orth_Scale_Invariant([], Zs, M_s, opts);
    [V_f_new, ~] = append_m_orth_Scale_Invariant([], Zf, M_f, opts);

    target_dim = min(opts.initDim, opts.restartDim);
    trials = 0;
    while size(V_s_new,2) < target_dim && trials < 50
        z = randn(ns,1); if is_complex_problem, z = z + 1i*randn(ns,1); end
        [V_s_new, add] = append_m_orth_Scale_Invariant(V_s_new, z, M_s, opts);
        trials = trials+1;
        if add == 0 && trials >= 50, break; end
    end

    trials = 0;
    while size(V_f_new,2) < target_dim && trials < 50
        z = randn(nf,1); if is_complex_problem, z = z + 1i*randn(nf,1); end
        [V_f_new, add] = append_m_orth_Scale_Invariant(V_f_new, z, M_f, opts);
        trials = trials+1;
        if add == 0 && trials >= 50, break; end
    end
end
