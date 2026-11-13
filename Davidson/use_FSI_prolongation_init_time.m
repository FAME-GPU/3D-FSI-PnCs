function [Vs0, Vf0, pinfo] = use_FSI_prolongation_init_time( ...
    mesh_c, phys, Opt, wave_vector, mesh_h, T_s_h, T_f_h, ...
    M_s_h, M_f_h, K_s_h, K_f_h)

if nargin < 10, K_s_h = []; end
if nargin < 11, K_f_h = []; end

solid_mode = 'random';
fluid_mode = 'coarse';
if isfield(Opt,'init')
    if isfield(Opt.init,'solid') && ~isempty(Opt.init.solid), solid_mode = lower(Opt.init.solid); end
    if isfield(Opt.init,'fluid') && ~isempty(Opt.init.fluid), fluid_mode = lower(Opt.init.fluid); end
end
validate_mode(solid_mode,'solid');
validate_mode(fluid_mode,'fluid');

flag_gpu = Opt.flag_gpu;
Vs0 = [];
Vf0 = [];
pinfo = struct();
pinfo.solid_mode = solid_mode;
pinfo.fluid_mode = fluid_mode;
pinfo.Phi_s_time = 0;
pinfo.Phi_f_time = 0;
pinfo.prolongation_time = 0;
pinfo.coarse_setup_time = 0;
pinfo.total_time = 0;
pinfo.efficiency = empty_efficiency();

all_timer = tic;
reset_init_efficiency_stats();

need_coarse_s = strcmp(solid_mode,'coarse');
need_coarse_f = strcmp(fluid_mode,'coarse');
need_fine_s = strcmp(solid_mode,'fine');
need_fine_f = strcmp(fluid_mode,'fine');

Phi_s_c = [];
Phi_f_c = [];
T_s_c = [];
T_f_c = [];
if need_coarse_s || need_coarse_f
    if isempty(mesh_c)
        error('A coarse mesh is required because coarse initialization was requested.');
    end
    tsetup = tic;
    mat_c = Assembly_mass_and_stiffness_matrix(mesh_c, phys.meteria);
    [T_s_c, T_f_c] = Generate_Bloch_T_Matrix_fast(mesh_c, wave_vector, phys.lattice.lattice_vector);
    K_s_c = hermitianize(T_s_c' * mat_c.K_s_tight * T_s_c);
    K_f_c = hermitianize(T_f_c' * mat_c.K_f_tight * T_f_c);
    M_s_c = hermitianize(T_s_c' * mat_c.M_s_tight * T_s_c);
    M_f_c = hermitianize(T_f_c' * mat_c.M_f_tight * T_f_c);
    pinfo.coarse_setup_time = toc(tsetup);

    if need_coarse_f
        t0 = tic;
        Phi_f_c = isolated_low_modes(K_f_c, M_f_c, Opt.subspace.init_dim, 'fluid', Opt, flag_gpu);
        pinfo.Phi_f_time = toc(t0);
    end
    if need_coarse_s
        t0 = tic;
        Phi_s_c = isolated_low_modes(K_s_c, M_s_c, Opt.subspace.init_dim, 'solid', Opt, flag_gpu);
        pinfo.Phi_s_time = toc(t0);
    end
end

if need_coarse_s || need_coarse_f
    popts = struct();
    popts.inputSpace = 'reduced';
    popts.massOrth = true;
    popts.verbose = isfield(Opt.init,'verbose') && Opt.init.verbose;
    popts.computeSolid = need_coarse_s;
    popts.computeFluid = need_coarse_f;

    t0 = tic;
    [Vs_c2f, Vf_c2f, prolong_info] = build_FSI_prolongation_init( ...
        mesh_c, mesh_h, T_s_c, T_f_c, T_s_h, T_f_h, ...
        Phi_s_c, Phi_f_c, M_s_h, M_f_h, popts);
    pinfo.prolongation_time = toc(t0);
    pinfo.prolongation = prolong_info;
    if need_coarse_s, Vs0 = Vs_c2f; end
    if need_coarse_f, Vf0 = Vf_c2f; end
end

if need_fine_s
    if isempty(K_s_h)
        error('K_s_h must be supplied when Opt.init.solid=''fine''.');
    end
    t0 = tic;
    Vs0 = isolated_low_modes(K_s_h, M_s_h, Opt.subspace.init_dim, 'solid', Opt, flag_gpu);
    pinfo.Phi_s_time = toc(t0);
    Vs0 = mass_orthonormalize_local(Vs0, M_s_h, Opt);
end

if need_fine_f
    if isempty(K_f_h)
        error('K_f_h must be supplied when Opt.init.fluid=''fine''.');
    end
    t0 = tic;
    Vf0 = isolated_low_modes(K_f_h, M_f_h, Opt.subspace.init_dim, 'fluid', Opt, flag_gpu);
    pinfo.Phi_f_time = toc(t0);
    Vf0 = mass_orthonormalize_local(Vf0, M_f_h, Opt);
end

pinfo.efficiency = collect_init_efficiency_stats();
pinfo.efficiency.isolated_mode_time = pinfo.Phi_s_time + pinfo.Phi_f_time;
pinfo.total_time = toc(all_timer);

if isfield(Opt.init,'verbose') && Opt.init.verbose
    fprintf('[Init] solid=%s, fluid=%s | Vs0=%d, Vf0=%d | total %.4f s\n', ...
            solid_mode, fluid_mode, size(Vs0,2), size(Vf0,2), pinfo.total_time);
end
end

function Phi = isolated_low_modes(K, M, nvec, domain, Opt, flag_gpu)

n = size(K,1);
if nvec <= 0
    Phi = zeros(n,0);
    return;
end

if flag_gpu
    v0 = rand(n,1,'double','gpuArray') + 1i*rand(n,1,'double','gpuArray');
else
    v0 = rand(n,1) + 1i*rand(n,1);
end
v0 = v0 / norm(v0);

if strcmpi(domain,'fluid')
    
    pre_struct = Prec_setup_Ks_ict(K, Opt.Lanczos.fluid.droptol, Opt.Lanczos.fluid.diagcomp, true);
    pre = pre_struct.apply;
    
else
    pre_struct = Prec_setup_Ks_ict(K, Opt.Lanczos.solid.droptol, Opt.Lanczos.solid.diagcomp, true);
    pre = pre_struct.apply;
end

fncK = @(x) counted_mv(K, x, domain, 'K');
fncM = @(x) counted_mv(M, x, domain, 'M');
fncKinv = @(x) silent_pcg(fncK, x, Opt.Lanczos.LSsolver.tol, ...
                           Opt.Lanczos.LSsolver.itmax, pre, domain);

lanczos_p = Opt.Lanczos.subspace_dim;
lanczos_p = max(lanczos_p, nvec + 4);
lanczos_p = min(lanczos_p, n);

[Phi, ~, ~] = Lanczos(fncM, fncK, fncKinv, n, nvec, ...
                      Opt.Lanczos.eigensolver.tol, Opt.Lanczos.eigensolver.itmax, ...
                      lanczos_p, v0, flag_gpu);
if flag_gpu, wait(gpuDevice); Phi = gather(Phi); end
end

function result = silent_pcg(Ax_func, rhs, tol, itmax, pre, domain)
global FSI_INIT_LIN_STATS;

mv0 = current_K_count(domain);
t0 = tic;
if isempty(pre)
    [result, flag, relres, iter] = pcg(Ax_func, rhs, tol, itmax);
elseif isa(pre,'function_handle')
    [result, flag, relres, iter] = pcg(Ax_func, rhs, tol, itmax, pre);
else
    [result, flag, relres, iter] = pcg(Ax_func, rhs, tol, itmax, pre, pre');
end
dt = toc(t0);
mv1 = current_K_count(domain);

fld = lower(domain);
FSI_INIT_LIN_STATS.(fld).calls = FSI_INIT_LIN_STATS.(fld).calls + 1;
FSI_INIT_LIN_STATS.(fld).iter_total = FSI_INIT_LIN_STATS.(fld).iter_total + iter;
FSI_INIT_LIN_STATS.(fld).time = FSI_INIT_LIN_STATS.(fld).time + dt;
FSI_INIT_LIN_STATS.(fld).matvec_total = FSI_INIT_LIN_STATS.(fld).matvec_total + (mv1-mv0);
FSI_INIT_LIN_STATS.(fld).last_flag = flag;
FSI_INIT_LIN_STATS.(fld).last_relres = relres;
if flag ~= 0
    FSI_INIT_LIN_STATS.(fld).failed_calls = FSI_INIT_LIN_STATS.(fld).failed_calls + 1;
end
end

function y = counted_mv(A, x, domain, kind)
global Kf_count Kf_time Ks_count Ks_time Mf_count Mf_time Ms_count Ms_time;
t0 = tic;
y = A*x;
dt = toc(t0);
ncol = size(x,2);
if strcmpi(domain,'fluid') && strcmpi(kind,'K')
    Kf_count = Kf_count + ncol; Kf_time = Kf_time + dt;
elseif strcmpi(domain,'fluid')
    Mf_count = Mf_count + ncol; Mf_time = Mf_time + dt;
elseif strcmpi(kind,'K')
    Ks_count = Ks_count + ncol; Ks_time = Ks_time + dt;
else
    Ms_count = Ms_count + ncol; Ms_time = Ms_time + dt;
end
end

function n = current_K_count(domain)
global Kf_count Ks_count;
if strcmpi(domain,'fluid'), n = Kf_count; else, n = Ks_count; end
end

function reset_init_efficiency_stats()
global Kf_count Kf_time Ks_count Ks_time Mf_count Mf_time Ms_count Ms_time;
global FSI_INIT_LIN_STATS;
Kf_count=0; Kf_time=0; Ks_count=0; Ks_time=0;
Mf_count=0; Mf_time=0; Ms_count=0; Ms_time=0;
empty = struct('calls',0,'iter_total',0,'time',0,'matvec_total',0, ...
               'failed_calls',0,'last_flag',NaN,'last_relres',NaN);
FSI_INIT_LIN_STATS = struct('fluid',empty,'solid',empty);
end

function eff = collect_init_efficiency_stats()
global Kf_count Kf_time Ks_count Ks_time Mf_count Mf_time Ms_count Ms_time;
global FSI_INIT_LIN_STATS;
eff = struct();
eff.linear = FSI_INIT_LIN_STATS;
eff.linear.total.calls = eff.linear.fluid.calls + eff.linear.solid.calls;
eff.linear.total.iter_total = eff.linear.fluid.iter_total + eff.linear.solid.iter_total;
eff.linear.total.time = eff.linear.fluid.time + eff.linear.solid.time;
eff.linear.total.failed_calls = eff.linear.fluid.failed_calls + eff.linear.solid.failed_calls;
eff.matvec.K_f = struct('count',Kf_count,'time',Kf_time);
eff.matvec.M_f = struct('count',Mf_count,'time',Mf_time);
eff.matvec.K_s = struct('count',Ks_count,'time',Ks_time);
eff.matvec.M_s = struct('count',Ms_count,'time',Ms_time);
eff.matvec.total_count = Kf_count + Mf_count + Ks_count + Ms_count;
eff.matvec.total_time = Kf_time + Mf_time + Ks_time + Ms_time;
end

function eff = empty_efficiency()
empty = struct('calls',0,'iter_total',0,'time',0,'matvec_total',0, ...
               'failed_calls',0,'last_flag',NaN,'last_relres',NaN);
eff = struct('linear',struct('fluid',empty,'solid',empty));
end

function Q = mass_orthonormalize_local(X, M, Opt)
if isempty(X), Q = X; return; end
orthTol = 1e-12;
maxReorth = 5;
if isfield(Opt,'orthTol'), orthTol = Opt.orthTol; end
if isfield(Opt,'maxReorth'), maxReorth = Opt.maxReorth; end
Q = zeros(size(X,1),0,'like',X);
for j = 1:size(X,2)
    v = X(:,j);
    for rep = 1:maxReorth
        if ~isempty(Q), v = v - Q*(Q'*(M*v)); end
    end
    nv = sqrt(abs(v'*(M*v)));
    if isfinite(nv) && nv > orthTol, Q(:,end+1) = v/nv; end 
end
end

function validate_mode(mode, domain)
if ~ismember(lower(mode), {'random','coarse','fine'})
    error('Unknown %s initialization mode: %s.', domain, mode);
end
end
