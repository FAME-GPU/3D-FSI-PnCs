function [z, st] = solve_shifted_system_time(K, M, theta, b, domain, opts)

A = K - theta^2 * M;
track_stats = get_opt(opts,'efficiencyStats',true);
track_matvec = get_opt(opts,'trackMatvec',track_stats);
allow_fallback = false;
if isfield(opts,'linear') && isfield(opts.linear,'allowFallback')
    allow_fallback = logical(opts.linear.allowFallback);
end

gamma_lift = false;
Z = [];
MZ = [];
tau = 0;
alpha = 0;
if isfield(opts,'gamma') && isfield(opts.gamma,'enabled') && opts.gamma.enabled && ...
        get_opt(opts.gamma,'useLift',true)
    if strcmpi(domain,'solid')
        if isfield(opts.gamma,'Zs'), Z = opts.gamma.Zs; end
        if isfield(opts.gamma,'MZs'), MZ = opts.gamma.MZs; end
        tau = get_opt(opts.gamma,'liftTau_s',1);
    elseif strcmpi(domain,'fluid')
        if isfield(opts.gamma,'Zf'), Z = opts.gamma.Zf; end
        if isfield(opts.gamma,'MZf'), MZ = opts.gamma.MZf; end
        tau = get_opt(opts.gamma,'liftTau_f',1);
    end

    if ~isempty(Z)
        if ~(isscalar(tau) && isfinite(tau) && tau > 0)
            error('Gamma lift parameter tau must be a finite positive scalar for %s.',domain);
        end
        
        if isempty(MZ)
            MZ = M*Z;
        end
        if size(MZ,1) ~= size(M,1) || size(MZ,2) ~= size(Z,2)
            error('Cached Gamma MZ has inconsistent dimensions for %s.',domain);
        end
        
        b = b - MZ*(Z'*b);
        alpha = tau + abs(theta)^2;
        gamma_lift = true;
    end
end

rhsnorm = norm(b);
st = struct('domain',domain,'theta',theta,'rhsnorm',rhsnorm, ...
            'flag',NaN,'relres',NaN,'iter',NaN, ...
            'solver',opts.linear.solver,'time',NaN, ...
            'matvec',NaN,'matvec_time',NaN, ...
            'fallback_used',false,'exception','');

if rhsnorm == 0 || ~isfinite(rhsnorm)
    z = zeros(size(b),'like',b);
    st.flag = 0;
    st.relres = 0;
    st.iter = 0;
    st.time = 0;
    st.matvec = 0;
    st.matvec_time = 0;
    return;
end

solver = lower(opts.linear.solver);
n = size(K,1);

if opts.linear.useDirect || strcmp(solver,'backslash')
    if track_stats, tsolve = tic; end
    if gamma_lift
        
        r = size(Z,2);
        Aug = [A, sparse(MZ); sparse(MZ'), sparse(r,r)];
        rhs_aug = [b; zeros(r,size(b,2),'like',b)];
        sol = Aug \ rhs_aug;
        z = sol(1:n,:);
        z = project_M_complement(z,Z,MZ);
        relres = norm(apply_operator(z)-b)/rhsnorm;
    else
        z = A \ b;
        relres = norm(A*z-b)/rhsnorm;
    end
    if track_stats, st.time = toc(tsolve); else, st.time = NaN; end
    st.flag = 0;
    st.relres = relres;
    st.iter = 0;
    st.matvec = 0;
    st.matvec_time = 0;
    st.solver = 'backslash';
    return;
end

Mfun = [];
if strcmpi(domain,'solid') && isfield(opts,'prec_s') && ~isempty(opts.prec_s)
    Mfun = opts.prec_s;
elseif strcmpi(domain,'fluid') && isfield(opts,'prec_f') && ~isempty(opts.prec_f)
    Mfun = opts.prec_f;
elseif opts.linear.useIC
    try
        setup.type = 'ict';
        setup.droptol = opts.linear.icholDroptol;
        if gamma_lift
            aux_shift = max(abs(theta)^2,tau);
        else
            aux_shift = abs(theta)^2;
        end
        L = ichol(K + aux_shift*M,setup);
        Mfun = @(x) L' \ (L \ x);
    catch
        Mfun = [];
    end
end

tol = opts.linear.tol;
maxit = opts.linear.maxit;

if track_matvec
    mv_count = 0;
    mv_time = 0;
    Aop = @counted_A_mv;
else
    mv_count = NaN;
    mv_time = NaN;
    if gamma_lift
        Aop = @apply_operator;
    else
        Aop = A;
    end
end

if track_stats, tsolve = tic; end
try
    switch solver
        case 'pcg'
            if isempty(Mfun)
                [z,flag,relres,iter] = pcg(Aop,b,tol,maxit);
            else
                [z,flag,relres,iter] = pcg(Aop,b,tol,maxit,Mfun);
            end
        case 'gmres'
            if isempty(Mfun)
                [z,flag,relres,iter] = gmres(Aop,b,[],tol,maxit);
            else
                [z,flag,relres,iter] = gmres(Aop,b,[],tol,maxit,Mfun);
            end
        otherwise
            
            if isempty(Mfun)
                [z,flag,relres,iter] = minres(Aop,b,tol,maxit);
            else
                [z,flag,relres,iter] = minres(Aop,b,tol,maxit,Mfun);
            end
    end

    if gamma_lift
        z = project_M_complement(z,Z,MZ);
        
        relres = norm(apply_operator(z)-b)/rhsnorm;
    end
catch ME
    st.exception = ME.message;
    if allow_fallback
        warning('Iterative shifted solve failed in %s block (%s). Using explicit fallback because allowFallback=true.', ...
                domain,ME.message);
        if gamma_lift
            r = size(Z,2);
            Aug = [A, sparse(MZ); sparse(MZ'), sparse(r,r)];
            rhs_aug = [b; zeros(r,size(b,2),'like',b)];
            sol = Aug \ rhs_aug;
            z = sol(1:n,:);
            z = project_M_complement(z,Z,MZ);
            relres = norm(apply_operator(z)-b)/rhsnorm;
        else
            z = A \ b;
            relres = norm(A*z-b)/rhsnorm;
        end
        flag = 0;
        iter = 0;
        solver = 'fallback_backslash';
        st.fallback_used = true;
    else
        warning('Iterative shifted solve failed in %s block (%s). No direct fallback is used.', ...
                domain,ME.message);
        z = zeros(size(b),'like',b);
        flag = -99;
        relres = Inf;
        iter = 0;
        solver = 'failed_iterative';
    end
end
if track_stats, st.time = toc(tsolve); else, st.time = NaN; end

st.flag = flag;
st.relres = relres;
st.iter = krylov_iter_count(iter,solver,n);
st.matvec = mv_count;
st.matvec_time = mv_time;
st.solver = solver;

    function y = apply_operator(x)
        y = A*x;
        if gamma_lift
            
            y = y + alpha*MZ*(MZ'*x);
        end
    end

    function y = counted_A_mv(x)
        tmv = tic;
        y = apply_operator(x);
        mv_time = mv_time + toc(tmv);
        mv_count = mv_count + size(x,2);
    end
end

function z = project_M_complement(z,Z,MZ)

z = z - Z*(MZ'*z);
end

function niter = krylov_iter_count(iter,solver,n)
if isempty(iter)
    niter = 0;
elseif isscalar(iter)
    niter = iter;
elseif strcmpi(solver,'gmres')
    restart = min(10,n);
    niter = max(0,iter(1)-1)*restart + iter(2);
else
    niter = iter(end);
end
end

function value = get_opt(s,name,default_value)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default_value;
end
end
