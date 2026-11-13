function result = main_eigs(model, TOL, num_level, opts, mesh_file, mesh_file_coarse)

addpath('Config'); 
addpath('FEM_Assembly'); 
addpath('Davidson'); 
addpath('data_paper')
addpath('mesh_data')
addpath('fun_eigs')

Opt = Initialize_parameters(model);
Opt.model = model;
if strcmp(model, 'SC')
    Opt.k_index = 10; 
elseif strcmp(model, 'SC_InnerFluid')
    Opt.k_index = 40;
end
Opt.nev = num_level;

Opt.eigs.tol = TOL;
Opt.eigs.inner_tol = Opt.eigs.tol;

phys = Physical_settings(model);

loaded_data = load(mesh_file);
meshData = loaded_data.meshData; loaded_data = [];
solidTag = 2; fluidTag = 1; cellSize = 1.0;
meshData_c = readGmsh(mesh_file_coarse, cellSize, solidTag, fluidTag);

fprintf('\n============================================================\n');
fprintf(' Single-k-point eigs paper baseline\n');
fprintf(' mesh      : %s\n', mesh_file);
fprintf(' k index   : %d\n', Opt.k_index);
fprintf(' nev       : %d\n', Opt.nev);
fprintf(' eigs tol  : %.3e\n', Opt.eigs.tol);
fprintf(' PCG tol   : %.3e\n', Opt.eigs.inner_tol);
fprintf(' PCG maxit : %d\n', Opt.eigs.inner_maxit);
fprintf(' profiling : %d\n', Opt.eigs.profile_efficiency);
fprintf('============================================================\n');

rng(1,'twister');

t_end_to_end = tic;

t_matrix = tic;
matrices = Assembly_mass_and_stiffness_matrix(meshData, phys.meteria);

wave_vector = phys.lattice.wave_vec_array(:,Opt.k_index);

[T_s,T_f] = Generate_Bloch_T_Matrix_fast( ...
meshData, wave_vector, phys.lattice.lattice_vector);

K_s = hermitianize(T_s' * matrices.K_s_tight * T_s);
K_f = hermitianize(T_f' * matrices.K_f_tight * T_f);
M_s = hermitianize(T_s' * matrices.M_s_tight * T_s);
M_f = hermitianize(T_f' * matrices.M_f_tight * T_f);
C   = T_f' * matrices.C_tight * T_s;

matrix_setup_time = toc(t_matrix);

n_s = size(K_s,1);
n_f = size(K_f,1);
n   = n_s + n_f;

fprintf('DOF: n_s=%d, n_f=%d, total=%d\n', n_s, n_f, n);
fprintf('Matrix setup time: %.4f s\n', matrix_setup_time);

stats = struct();
stats.eigs_operator_calls = 0;
stats.kf_pcg_calls = 0;
stats.ks_pcg_calls = 0;
stats.kf_iter_total = 0;
stats.ks_iter_total = 0;
stats.kf_matvec_count = 0;
stats.ks_matvec_count = 0;
stats.kf_matvec_time = 0;
stats.ks_matvec_time = 0;
stats.kf_linear_time = 0;
stats.ks_linear_time = 0;
stats.kf_prec_apply_count = 0;
stats.ks_prec_apply_count = 0;
stats.kf_prec_apply_time = 0;
stats.ks_prec_apply_time = 0;

inner_failure_count = 0;

t_method = tic;

if strcmp(model, 'SC')
    solid_diagcomp = 0;    solid_droptol = 1e-3;
    fluid_diagcomp = 1e-5; fluid_droptol = 5e-3;
elseif strcmp(model, 'SC_InnerFluid')
    solid_diagcomp = 1e-6;    solid_droptol = 1e-4;
    fluid_diagcomp = 1e-5; fluid_droptol = 1e-3;
end
t_prec = tic;
P_f = Prec_setup_Ks_ict(K_f, fluid_droptol, fluid_diagcomp, true);
P_s = Prec_setup_Ks_ict(K_s, solid_droptol, solid_diagcomp, true);
preconditioner_setup_time = toc(t_prec);

opts_eigs = struct();
opts_eigs.isreal = false;
opts_eigs.tol = Opt.eigs.tol;
opts_eigs.disp = 1;

[opts_eigs.v0, initial_info] = init_eigs_v0_coupled_coarse(meshData_c, meshData, phys, wave_vector, T_s, T_f, opts);

V = [];
D = [];
error_message = '';

t_eigs = tic;

[V,D] = eigs(@apply_AinvB, n, Opt.nev, 'largestabs', opts_eigs);
eigs_solver_time = toc(t_eigs);

solver_total_time = toc(t_method); 
fprintf('Preconditioner setup time: %.4f s\n', preconditioner_setup_time);
fprintf('eigs solver time         : %.4f s\n', eigs_solver_time);
fprintf('Solver total time        : %.4f s\n', solver_total_time);

t_post = tic;
lambda = [];
omega = [];
frequency_Hz = [];
residuals = [];
Rmax = NaN;

if ~isempty(D)
mu = diag(D);
lam_all = 1 ./ mu;

imag_tol = 1e-7 .* max(1,abs(lam_all));
good = isfinite(lam_all) & real(lam_all) > 0 & ...
       abs(imag(lam_all)) <= imag_tol;

lam_all = real(lam_all(good));
Vgood = V(:,good);

[lambda,ord] = sort(lam_all,'ascend');
Vgood = Vgood(:,ord);

nkeep = min(Opt.nev,numel(lambda));
lambda = lambda(1:nkeep);
Vgood = Vgood(:,1:nkeep);

U_s = Vgood(1:n_s,:);
P   = Vgood(n_s+1:end,:);

omega = sqrt(lambda);

[residuals,Rmax] = compute_FSI_GEP_residual( ...
    K_s,K_f,M_s,M_f,C,lambda,U_s,P);
end

residual_evaluation_time = toc(t_post);
end_to_end_time = toc(t_end_to_end);
converged = numel(lambda) >= Opt.nev;
success = isempty(error_message) && converged && ...
      inner_failure_count == 0 && isfinite(Rmax) && ...
      Rmax <= Opt.eigs.tol;

metrics = struct();
metrics.method = 'eigs';
metrics.mesh_file = mesh_file;
metrics.k_index = Opt.k_index;
metrics.k_norm = norm(wave_vector);
metrics.n_s = n_s;
metrics.n_f = n_f;
metrics.n_total = n;

metrics.matrix_setup_time = matrix_setup_time;
metrics.preconditioner_setup_time = preconditioner_setup_time;
metrics.eigs_solver_time = eigs_solver_time;
metrics.arnoldi_iterations = stats.eigs_operator_calls;
metrics.eigs_operator_calls = stats.eigs_operator_calls;
metrics.total_time = solver_total_time; 
metrics.residual_evaluation_time = residual_evaluation_time;
metrics.end_to_end_time = end_to_end_time;

metrics.K_f_preconditioner_apply_time = stats.kf_prec_apply_time;
metrics.K_s_preconditioner_apply_time = stats.ks_prec_apply_time;
metrics.total_preconditioner_apply_time = stats.kf_prec_apply_time + stats.ks_prec_apply_time;
metrics.K_f_preconditioner_apply_count = stats.kf_prec_apply_count;
metrics.K_s_preconditioner_apply_count = stats.ks_prec_apply_count;
metrics.total_preconditioner_apply_count = stats.kf_prec_apply_count + stats.ks_prec_apply_count;

metrics.final_R_max = Rmax;
metrics.converged = converged;
metrics.success = success;
metrics.inner_failure_count = inner_failure_count;
metrics.eigs_tolerance = Opt.eigs.tol;
metrics.inner_tolerance = Opt.eigs.inner_tol;
metrics.inner_maxit = Opt.eigs.inner_maxit;
metrics.target_tolerance = Opt.eigs.tol;
metrics.eigs.profile_efficiency = Opt.eigs.profile_efficiency;

if Opt.eigs.profile_efficiency
    efficiency = struct();
    efficiency.arnoldi_iterations = stats.eigs_operator_calls;
    efficiency.eigs_operator_calls = stats.eigs_operator_calls;
    efficiency.K_f_PCG_calls = stats.kf_pcg_calls;
    efficiency.K_s_PCG_calls = stats.ks_pcg_calls;
    efficiency.total_PCG_calls = stats.kf_pcg_calls + stats.ks_pcg_calls;
    efficiency.K_f_PCG_iterations = stats.kf_iter_total;
    efficiency.K_s_PCG_iterations = stats.ks_iter_total;
    efficiency.total_krylov_iterations = stats.kf_iter_total + stats.ks_iter_total;
    efficiency.K_f_matvec_count = stats.kf_matvec_count;
    efficiency.K_s_matvec_count = stats.ks_matvec_count;
    efficiency.total_operator_applications = stats.kf_matvec_count + stats.ks_matvec_count;
    efficiency.K_f_matvec_time = stats.kf_matvec_time;
    efficiency.K_s_matvec_time = stats.ks_matvec_time;
    efficiency.total_matvec_time = stats.kf_matvec_time + stats.ks_matvec_time;
    efficiency.K_f_linear_time = stats.kf_linear_time;
    efficiency.K_s_linear_time = stats.ks_linear_time;
    efficiency.total_linear_solve_time = stats.kf_linear_time + stats.ks_linear_time;
    efficiency.K_f_preconditioner_apply_count = stats.kf_prec_apply_count;
    efficiency.K_s_preconditioner_apply_count = stats.ks_prec_apply_count;
    efficiency.total_preconditioner_apply_count = stats.kf_prec_apply_count + stats.ks_prec_apply_count;
    efficiency.K_f_preconditioner_apply_time = stats.kf_prec_apply_time;
    efficiency.K_s_preconditioner_apply_time = stats.ks_prec_apply_time;
    efficiency.total_preconditioner_apply_time = stats.kf_prec_apply_time + stats.ks_prec_apply_time;
else
    efficiency = struct();
end

result = struct();
result.Opt = Opt;
result.wave_vector = wave_vector;
result.eigenvalues = lambda;
result.angular_frequency = omega;
result.frequency_Hz = frequency_Hz;
result.final_GEP_residuals = residuals;
result.metrics = metrics;
result.efficiency = efficiency;
result.error_message = error_message;
result.initial_info = initial_info;
if strcmp(model, 'SC')
    folderPath = fullfile('data_paper');
end
fileName = sprintf('eigs_%s.mat', model);
fullFilePath = fullfile(folderPath, fileName);
save(fullFilePath, 'result');

fprintf('\n-------------------- Result summary --------------------\n');
fprintf('Valid eigenvalues : %d / %d\n', numel(lambda), Opt.nev);
fprintf('Final GEP Rmax    : %.6e\n', Rmax);
fprintf('Arnoldi steps     : %d\n', stats.eigs_operator_calls);
fprintf('Residual time     : %.4f s\n', residual_evaluation_time);
fprintf('End-to-end time   : %.4f s\n', end_to_end_time);
fprintf('K_f prec apply    : %d calls, %.4f s\n', stats.kf_prec_apply_count, stats.kf_prec_apply_time);
fprintf('K_s prec apply    : %d calls, %.4f s\n', stats.ks_prec_apply_count, stats.ks_prec_apply_time);
fprintf('Inner failures    : %d\n', inner_failure_count);
fprintf('Success           : %d\n', success);

if ~isempty(frequency_Hz)
fprintf('\nEigenfrequencies (Hz):\n');
disp(frequency_Hz);
end

if Opt.eigs.profile_efficiency
fprintf('\nDetailed efficiency profiling:\n');
fprintf('  eigs operator calls      : %d\n', efficiency.eigs_operator_calls);
fprintf('  total PCG calls          : %d\n', efficiency.total_PCG_calls);
fprintf('  total Krylov iterations  : %d\n', efficiency.total_krylov_iterations);
fprintf('  total matvec applications: %d\n', efficiency.total_operator_applications);
fprintf('  total matvec time        : %.4f s\n', efficiency.total_matvec_time);
fprintf('  total prec-apply time    : %.4f s\n', efficiency.total_preconditioner_apply_time);
fprintf('  total linear-solve time  : %.4f s\n', efficiency.total_linear_solve_time);
end
fprintf('--------------------------------------------------------\n');
    function y = apply_AinvB(v)
    ncol = size(v,2);
    y = zeros(size(v),'like',v);
    
    if Opt.eigs.profile_efficiency
        stats.eigs_operator_calls = stats.eigs_operator_calls + ncol;
    end
    
    for jj = 1:ncol
        u = v(1:n_s,jj);
        p = v(n_s+1:end,jj);
        
        rhs_p = C*u + M_f*p;
        if Opt.eigs.profile_efficiency
            t0 = tic;
            [y_p,flag_f,~,iter_f] = pcg( ...
                @Kf_mv, rhs_p, Opt.eigs.inner_tol, ...
                Opt.eigs.inner_maxit, @Pf_apply_timed);
            stats.kf_linear_time = stats.kf_linear_time + toc(t0);
            stats.kf_pcg_calls = stats.kf_pcg_calls + 1;
            stats.kf_iter_total = stats.kf_iter_total + iter_count(iter_f);
        else
            [y_p,flag_f] = pcg( ...
                K_f, rhs_p, Opt.eigs.inner_tol, ...
                Opt.eigs.inner_maxit, @Pf_apply_timed);
        end
        
        if flag_f ~= 0
            inner_failure_count = inner_failure_count + 1;
            warning('K_f PCG exited with flag %d at eigs operator call.', flag_f);
        end
        
        rhs_u = M_s*u + C'*y_p;
        if Opt.eigs.profile_efficiency
            t0 = tic;
            [y_u,flag_s,~,iter_s] = pcg( ...
                @Ks_mv, rhs_u, Opt.eigs.inner_tol, ...
                Opt.eigs.inner_maxit, @Ps_apply_timed);
            stats.ks_linear_time = stats.ks_linear_time + toc(t0);
            stats.ks_pcg_calls = stats.ks_pcg_calls + 1;
            stats.ks_iter_total = stats.ks_iter_total + iter_count(iter_s);
        else
            [y_u,flag_s] = pcg( ...
                K_s, rhs_u, Opt.eigs.inner_tol, ...
                Opt.eigs.inner_maxit, @Ps_apply_timed);
        end
        
        if flag_s ~= 0
            inner_failure_count = inner_failure_count + 1;
            warning('K_s PCG exited with flag %d at eigs operator call.', flag_s);
        end
        
        y(:,jj) = [y_u;y_p];
    end
end

function y = Pf_apply_timed(x)
    t0 = tic;
    y = P_f.apply(x);
    stats.kf_prec_apply_time = stats.kf_prec_apply_time + toc(t0);
    stats.kf_prec_apply_count = stats.kf_prec_apply_count + 1;
end

function y = Ps_apply_timed(x)
    t0 = tic;
    y = P_s.apply(x);
    stats.ks_prec_apply_time = stats.ks_prec_apply_time + toc(t0);
    stats.ks_prec_apply_count = stats.ks_prec_apply_count + 1;
end

function y = Kf_mv(x)
    t0 = tic;
    y = K_f*x;
    stats.kf_matvec_time = stats.kf_matvec_time + toc(t0);
    stats.kf_matvec_count = stats.kf_matvec_count + size(x,2);
end

function y = Ks_mv(x)
    t0 = tic;
    y = K_s*x;
    stats.ks_matvec_time = stats.ks_matvec_time + toc(t0);
    stats.ks_matvec_count = stats.ks_matvec_count + size(x,2);
end

function n = iter_count(iter)
    if isempty(iter)
        n = 0;
    elseif isscalar(iter)
        n = iter;
    else
        n = iter(end);
    end
end

end