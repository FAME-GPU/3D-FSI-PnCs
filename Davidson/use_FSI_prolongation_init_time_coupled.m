function [Vs0, Vf0, initial_info] = use_FSI_prolongation_init_time_coupled( ...
    mesh_c, mesh_h, phys, wave_vector, T_s_h, T_f_h, M_s_h, M_f_h, opts)

    if nargin < 8
        error(['M_s_h and M_f_h (fine reduced mass matrices) are required. ', ...
               'Call with 8 inputs plus optional opts.']);
    end
    if nargin < 9 || isempty(opts), opts = struct(); end
    opts = set_defaults(opts);
    check_positive_integer(opts.nEig, 'opts.nEig');
    if ~isnumeric(wave_vector) || numel(wave_vector) ~= 3
        error('wave_vector must contain three components.');
    end
    if ~isfield(phys,'meteria') || ~isfield(phys,'lattice')
        error('phys must be the output of Physical_settings(model).');
    end
    if isempty(M_s_h) || isempty(M_f_h)
        error('Both M_s_h and M_f_h are required for fine-grid M-orthonormalization.');
    end
    if size(M_s_h,1) ~= size(M_s_h,2) || size(M_f_h,1) ~= size(M_f_h,2)
        error('M_s_h and M_f_h must be square matrices.');
    end

    t_all = tic;
    timing = struct();
    t = tic;
    mesh_c = read_mesh(mesh_c);
    mesh_h = read_mesh(mesh_h);
    timing.mesh_load = toc(t);

    t = tic;
    mat_c = Assembly_mass_and_stiffness_matrix(mesh_c, phys.meteria);
    [T_s_c,T_f_c] = Generate_Bloch_T_Matrix_fast( ...
        mesh_c, wave_vector, phys.lattice.lattice_vector);
    K_s = hermitian_part(T_s_c' * mat_c.K_s_tight * T_s_c);
    K_f = hermitian_part(T_f_c' * mat_c.K_f_tight * T_f_c);
    M_s = hermitian_part(T_s_c' * mat_c.M_s_tight * T_s_c);
    M_f = hermitian_part(T_f_c' * mat_c.M_f_tight * T_f_c);
    C   = T_f_c' * mat_c.C_tight * T_s_c;
    clear mat_c
    timing.coarse_assembly = toc(t);

    ns = size(K_s,1);
    nf = size(K_f,1);
    n = ns + nf;
    if n <= opts.nEig+1
        error('Coarse problem has n=%d DOFs; nEig=%d is too large.', n, opts.nEig);
    end
    if opts.verbose
        fprintf('[Coupled->Davidson] coarse DOFs: solid=%d, fluid=%d; nEig=%d\n', ...
            ns, nf, opts.nEig);
    end

    t = tic;
    P_s = Prec_setup_Ks_ict(K_s, opts.solidDroptol, opts.solidDiagcomp, true);
    P_f = Prec_setup_Ks_ict(K_f, opts.fluidDroptol, opts.fluidDiagcomp, true);
    timing.preconditioner = toc(t);

    inner = struct('operator_applications',0, ...
        'solid_pcg_calls',0,'fluid_pcg_calls',0, ...
        'solid_pcg_iterations',0,'fluid_pcg_iterations',0, ...
        'solid_pcg_failures',0,'fluid_pcg_failures',0, ...
        'solid_pcg_time',0,'fluid_pcg_time',0);

    eopts = struct('isreal',false,'issym',false,'tol',opts.eigsTol, ...
                   'maxit',opts.maxit,'disp',0);
    eopts.p = choose_p(n, opts.nEig, opts.p);
    t = tic;
    [V,D,eigs_flag] = eigs(@apply_AinvB, n, opts.nEig, 'largestabs', eopts);
    timing.coarse_eigs = toc(t);

    [lambda_c, V_c] = positive_real_modes(1./diag(D), V, opts.realTol);
    if numel(lambda_c) < opts.nEig
        error(['Only %d / %d eigenpairs have positive essentially-real eigenvalues. ', ...
               'No complete Davidson initialization was constructed; ', ...
               'check eigsTol, innerTol and the coarse problem.'], ...
               numel(lambda_c), opts.nEig);
    end
    if eigs_flag ~= 0
        warning('Coarse eigs flag=%d; inspect initial_info.coarse_residuals.', eigs_flag);
    end
    phi_s_c = V_c(1:ns, 1:opts.nEig);
    phi_f_c = V_c(ns+1:end, 1:opts.nEig);

    t = tic;
    coarse_residuals = zeros(opts.nEig,1);
    for j = 1:opts.nEig
        us = phi_s_c(:,j);
        pf = phi_f_c(:,j);
        lam = lambda_c(j);
        ax = [K_s*us-C'*pf; K_f*pf];
        bx = [M_s*us; C*us+M_f*pf];
        coarse_residuals(j) = norm(ax-lam*bx) / ...
            (norm(ax) + abs(lam)*norm(bx) + eps);
    end
    timing.coarse_residual = toc(t);

    t = tic;
    if isempty(T_s_h) || isempty(T_f_h)
        [T_s_new,T_f_new] = Generate_Bloch_T_Matrix_fast( ...
            mesh_h, wave_vector, phys.lattice.lattice_vector);
        if isempty(T_s_h), T_s_h = T_s_new; end
        if isempty(T_f_h), T_f_h = T_f_new; end
    end
    timing.fine_bloch = toc(t);
    if size(T_s_h,2) ~= size(M_s_h,1) || size(T_f_h,2) ~= size(M_f_h,1)
        error(['Fine Bloch maps and fine mass matrices have incompatible dimensions. ', ...
               'Use the same wave_vector and fine mesh as in Davidson.']);
    end

    popts = opts.prolongOpts;
    popts.inputSpace = 'reduced';
    popts.massOrth = false; 
    popts.computeSolid = true;
    popts.computeFluid = true;
    popts.verbose = opts.verbose && get_option(popts,'verbose',false);
    t = tic;
    [U_h,P_h,prolong_info] = build_FSI_prolongation_init( ...
        mesh_c, mesh_h, T_s_c, T_f_c, T_s_h, T_f_h, ...
        phi_s_c, phi_f_c, [], [], popts);
    timing.prolongation = toc(t);

    if size(U_h,2) ~= opts.nEig || size(P_h,2) ~= opts.nEig
        error('Prolongation did not return nEig columns for BOTH subsystems.');
    end
    
    block_norm_s = sqrt(sum(abs(U_h).^2,1)).';
    block_norm_f = sqrt(sum(abs(P_h).^2,1)).';

    t = tic;
    [Vs0, solid_orth] = orthonormalize_mass(U_h, M_s_h, opts.orthTol, 'solid');
    timing.solid_orth = toc(t);
    t = tic;
    [Vf0, fluid_orth] = orthonormalize_mass(P_h, M_f_h, opts.orthTol, 'fluid');
    timing.fluid_orth = toc(t);

    if size(Vs0,2) ~= opts.nEig || size(Vf0,2) ~= opts.nEig
        error('Each Davidson initial basis must have exactly nEig columns.');
    end

    initial_info = struct();
    initial_info.method = 'coarse_full_coupled_eigs_to_separated_Davidson';
    initial_info.options = opts;
    initial_info.coarse_dofs = [ns,nf];
    initial_info.fine_dofs = [size(T_s_h,2),size(T_f_h,2)];
    initial_info.nEig_requested = opts.nEig;
    initial_info.nSolid_returned = size(Vs0,2);
    initial_info.nFluid_returned = size(Vf0,2);
    initial_info.coarse_eigs_flag = eigs_flag;
    initial_info.coarse_lambda = lambda_c;
    initial_info.coarse_residuals = coarse_residuals;
    initial_info.coarse_max_residual = max(coarse_residuals);
    initial_info.prolongated_solid_column_norms = block_norm_s;
    initial_info.prolongated_fluid_column_norms = block_norm_f;
    initial_info.orthogonalization.solid = solid_orth;
    initial_info.orthogonalization.fluid = fluid_orth;
    initial_info.inner = inner;
    initial_info.inner.total_pcg_calls = inner.solid_pcg_calls + inner.fluid_pcg_calls;
    initial_info.inner.total_pcg_iterations = ...
        inner.solid_pcg_iterations + inner.fluid_pcg_iterations;
    initial_info.prolongation = prolong_info;
    initial_info.timings = timing;
    initial_info.total_time = toc(t_all); 
    initial_info.save_file = char(opts.saveFile);

    if ~isempty(opts.saveFile)
        output_dir = fileparts(char(opts.saveFile));
        if ~isempty(output_dir) && ~isfolder(output_dir), mkdir(output_dir); end
        save(char(opts.saveFile), 'Vs0', 'Vf0', 'initial_info', '-v7.3');
    end
    if opts.verbose
        fprintf(['[Coupled->Davidson] Vs0: %d x %d, Vf0: %d x %d; ', ...
            'M-orth errors: %.2e / %.2e; total %.3f s\n'], ...
            size(Vs0,1),size(Vs0,2),size(Vf0,1),size(Vf0,2), ...
            solid_orth.error,fluid_orth.error,initial_info.total_time);
        if ~isempty(opts.saveFile)
            fprintf('[Coupled->Davidson] saved: %s\n',char(opts.saveFile));
        end
    end

    function y = apply_AinvB(x)
        y = zeros(size(x),'like',x);
        for col = 1:size(x,2)
            u = x(1:ns,col);
            p = x(ns+1:end,col);

            rhs_p = C*u + M_f*p;
            tc = tic;
            [yp,flag_f,~,it_f] = pcg(K_f,rhs_p,opts.innerTol, ...
                                        opts.innerMaxit,P_f.apply);
            inner.fluid_pcg_time = inner.fluid_pcg_time + toc(tc);
            inner.fluid_pcg_calls = inner.fluid_pcg_calls + 1;
            inner.fluid_pcg_iterations = inner.fluid_pcg_iterations + count_iter(it_f);
            inner.fluid_pcg_failures = inner.fluid_pcg_failures + (flag_f~=0);

            rhs_u = M_s*u + C'*yp;
            tc = tic;
            [yu,flag_s,~,it_s] = pcg(K_s,rhs_u,opts.innerTol, ...
                                        opts.innerMaxit,P_s.apply);
            inner.solid_pcg_time = inner.solid_pcg_time + toc(tc);
            inner.solid_pcg_calls = inner.solid_pcg_calls + 1;
            inner.solid_pcg_iterations = inner.solid_pcg_iterations + count_iter(it_s);
            inner.solid_pcg_failures = inner.solid_pcg_failures + (flag_s~=0);

            y(:,col) = [yu;yp];
            inner.operator_applications = inner.operator_applications + 1;
        end
    end
end

function opts = set_defaults(opts)
    opts = set_default(opts,'nEig',8);
    opts = set_default(opts,'eigsTol',1e-2);
    opts = set_default(opts,'innerTol',1e-3);
    opts = set_default(opts,'innerMaxit',200);
    opts = set_default(opts,'maxit',300);
    opts = set_default(opts,'p',[]);
    opts = set_default(opts,'realTol',1e-5);
    opts = set_default(opts,'orthTol',1e-12);
    opts = set_default(opts,'solidDroptol',1e-3);
    opts = set_default(opts,'solidDiagcomp',0);
    opts = set_default(opts,'fluidDroptol',5e-3);
    opts = set_default(opts,'fluidDiagcomp',1e-5);
    opts = set_default(opts,'verbose',true);
    opts = set_default(opts,'prolongOpts',struct());
    
    if ~isfield(opts,'saveFile')
        opts.saveFile = 'Davidson_init_coupled_coarse.mat';
    end
end

function s = set_default(s,name,val)
    if ~isfield(s,name) || isempty(s.(name)), s.(name) = val; end
end

function val = get_option(s,name,default)
    if isfield(s,name) && ~isempty(s.(name)), val = s.(name);
    else, val = default; end
end

function mesh = read_mesh(mesh)
    if ischar(mesh) || isstring(mesh)
        file = char(mesh);
        if ~isfile(file), error('Mesh file not found: %s',file); end
        data = load(file,'meshData');
        if ~isfield(data,'meshData')
            error('MAT file must contain variable meshData: %s',file);
        end
        mesh = data.meshData;
    end
    if isstruct(mesh) && isfield(mesh,'meshData') && ~isfield(mesh,'nodes')
        mesh = mesh.meshData;
    end
    if ~isstruct(mesh) || ~isfield(mesh,'nodes') || ~isfield(mesh,'elements')
        error('mesh_c/mesh_h must be meshData structs or MAT paths.');
    end
end

function H = hermitian_part(H)
    H = (H+H')*0.5;
end

function check_positive_integer(x,name)
    if ~isnumeric(x) || ~isscalar(x) || ~isfinite(x) || x<1 || x~=floor(x)
        error('%s must be a positive integer.',name);
    end
end

function p = choose_p(n,k,p_user)
    if isempty(p_user), p=min(n,max(20,3*k));
    else, p=p_user; end
    if ~isnumeric(p)||~isscalar(p)||~isfinite(p)||p<=k+1||p>n||p~=floor(p)
        error('Invalid opts.p; must satisfy nEig+1 < p <= coarse total DOFs.');
    end
end

function [lambda,X] = positive_real_modes(lambda,X,rtol)
    lambda = lambda(:);
    valid = isfinite(real(lambda)) & isfinite(imag(lambda)) & real(lambda)>0 & ...
        abs(imag(lambda)) <= rtol.*max(1,abs(real(lambda)));
    lambda = real(lambda(valid));
    X = X(:,valid);
    [lambda,order] = sort(lambda,'ascend');
    X = X(:,order);
end

function [Q,diagnostics] = orthonormalize_mass(X,M,tol,domain)

    k = size(X,2);
    colnorm = sqrt(sum(abs(X).^2,1));
    if any(~isfinite(colnorm)) || any(colnorm<=0)
        error('The prolonged %s block contains a zero/nonfinite column.',domain);
    end
    X = bsxfun(@rdivide,X,colnorm);
    G = hermitian_part(X'*(M*X));
    [W,d] = eig(full(G),'vector');
    d = real(d(:));
    if any(~isfinite(d)) || max(d)<=0
        error('%s M-Gram is not positive definite.',domain);
    end
    threshold = tol*max(d);
    if any(d<=threshold)
        error(['The prolonged %s modes have numerical rank %d < nEig=%d ', ...
            'in the mass inner product (orthTol=%.1e). ', ...
            'Reduce nEig or inspect coupled modes / prolongation.'], ...
            domain,nnz(d>threshold),k,tol);
    end
    Q = (X*W)*diag(1./sqrt(d));
    
    G2 = hermitian_part(Q'*(M*Q));
    [R,cholflag] = chol(full(G2));
    if cholflag ~= 0
        error('Second-pass %s M-Gram is not positive definite.',domain);
    end
    Q = Q/R;
    err = norm(hermitian_part(Q'*(M*Q))-eye(k),'fro');
    if ~isfinite(err) || err > 1e-7
        warning('Large M-orthogonality error in %s initial basis: %.3e.',domain,err);
    end
    diagnostics = struct('rank',k,'error',err,'gram_eigenvalues',d, ...
        'gram_condition',max(d)/min(d), ...
        'column_norms_before_scaling',colnorm(:));
end

function n = count_iter(it)
    if isempty(it), n=0;
    elseif isscalar(it), n=it;
    else, n=it(end); end
end
