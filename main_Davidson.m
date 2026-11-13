clc;clear;close all

addpath('Config'); 
addpath('FEM_Assembly'); 
addpath('Davidson'); 
addpath('data_paper')
addpath('mesh_data')

method_timer = tic;
model = 'SC_Quad';
Opt = Initialize_parameters(model);

if strcmp(Opt.init.fluid, 'coarse_coupled')
    coupled_initial = 1;
else
    coupled_initial = 0;
end
coupled_nEig = Opt.subspace.init_dim;

if strcmp(Opt.model, 'SC')
    solidTag = 2; fluidTag = 1; cellSize = 1.0;
    phys = Physical_settings(Opt.model);
    load('SC_refined.mat')
    meshData_c = readGmsh('SC_coarse.msh', cellSize, solidTag, fluidTag);
elseif strcmp(Opt.model, 'AlWater') 
    solidTag = 2; fluidTag = 1; 
    phys = Physical_settings_AlWater(Opt.model);
    cellSize = phys.geometry.cell_size;
    load('AlWater_refined.mat')
    meshData_c = readGmsh('AlWater_coarse.msh',   cellSize, solidTag, fluidTag);
elseif strcmp(Opt.model, 'SC_Quad')
    solidTag = 1; fluidTag = 2; cellSize = 1.0;
    phys = Physical_settings(Opt.model); 
    load('Quad_refined.mat')
    meshData_c = readGmsh('Quad_coarse.msh', cellSize, solidTag, fluidTag);
end

matrices = Assembly_mass_and_stiffness_matrix(meshData, phys.meteria);

wave_vector = phys.lattice.wave_vec_array(:,Opt.k_index);
[T_s,T_f] = Generate_Bloch_T_Matrix_fast(meshData,wave_vector,phys.lattice.lattice_vector);

K_s = hermitianize(T_s'*matrices.K_s_tight*T_s);
K_f = hermitianize(T_f'*matrices.K_f_tight*T_f);
M_s = hermitianize(T_s'*matrices.M_s_tight*T_s);
M_f = hermitianize(T_f'*matrices.M_f_tight*T_f);
C = T_f'*matrices.C_tight*T_s;
clear matrices;

mesh_summary.nSolid=size(K_s,1); 
mesh_summary.nFluid=size(K_f,1); 
mesh_summary.k_index=Opt.k_index;

rng(Opt.random_seed);

Opts = davidson_opts_from_Opt(Opt);

pinfo = empty_pinfo(Opt);
preconditioner_setup_time = 0;

is_complex_problem = ~isreal(K_s) || ~isreal(K_f) || ...
                     ~isreal(M_s) || ~isreal(M_f) || ~isreal(C);

if strcmpi(Opt.init.solid,'random')
    t_solid_intial = tic;
    Opts.Vs0 = fixed_random_block(size(K_s,1), ...
        Opt.subspace.init_dim, Opt.random_seed+11, is_complex_problem);
    Opts.Vs0 = m_orthonormalize(Opts.Vs0, M_s, 1e-12, 5);
    solid_intial_time = toc(t_solid_intial);
end

if strcmpi(Opt.init.fluid,'random')
    t_fluid_intial = tic;
    Opts.Vf0 = fixed_random_block(size(K_f,1), ...
        Opt.subspace.init_dim, Opt.random_seed+29, is_complex_problem);
    Opts.Vf0 = m_orthonormalize(Opts.Vf0, M_f, 1e-12, 5);
    fluid_intial_time = toc(t_fluid_intial);
end

if Opt.preconditioner.enable
    tprec = tic;

    Ps = Prec_setup_solid_ict(K_s,M_s,...
        Opt.preconditioner.solid.droptol,...
        Opt.preconditioner.solid.lambda,...
        Opt.preconditioner.solid.use_reorder);

    Pf = Prec_setup_fluid_ict(K_f,M_f,...
        Opt.preconditioner.fluid.droptol,...
        Opt.preconditioner.fluid.lambda,...
        Opt.preconditioner.fluid.use_reorder);

    preconditioner_setup_time = toc(tprec);
    Opts.prec_s = Ps.apply;
    Opts.prec_f = Pf.apply;
else
    Opts.prec_s = [];
    Opts.prec_f = [];
end

need_init = ~strcmpi(Opt.init.solid,'random') || ...
            ~strcmpi(Opt.init.fluid,'random');

if need_init
    if coupled_initial
        opts_init = struct();
        opts_init.nEig = coupled_nEig;
        opts_init.solidDroptol = 1e-4;
        opts_init.fluidDroptol = 1e-3;
        opts_init.fluidDiagcomp = 1e-5;
        [Vs0, Vf0, pinfo] = use_FSI_prolongation_init_time_coupled(meshData_c, meshData, phys, wave_vector, T_s, T_f, M_s, M_f, opts_init);
    else
        [Vs0,Vf0,pinfo] = use_FSI_prolongation_init_time(meshData_c, phys,  Opt, wave_vector, meshData,  T_s,T_f,M_s,M_f,K_s,K_f);
    end
    if ~isempty(Vs0), Opts.Vs0 = Vs0; end
    if ~isempty(Vf0), Opts.Vf0 = Vf0; end
end

[lambda,omega,U_s,U_f,P,info] = solve_FSI_separated_Davidson_time_lock(...
    K_s, K_f, M_s, M_f, C, Opt.nev, Opt.target_omega, Opts);

total_time = toc(method_timer);

result = base_result(Opt,mesh_summary);
result.lambda = lambda;
result.omega = omega;
result.metrics = build_metrics(Opt, info, pinfo, preconditioner_setup_time,total_time);
result.history = info.history;
result.initialization = pinfo;
result.info = info;
result.error_message = '';
if strcmpi(Opt.init.solid,'random')
    result.solid_intial_time = solid_intial_time;
end
if strcmpi(Opt.init.fluid,'random')
    result.fluid_intial_time = fluid_intial_time;
end

folderPath = fullfile('data_paper');  
fileName = sprintf('result_%s.mat', model);
fullFilePath = fullfile(folderPath, fileName);
save(fullFilePath, 'result');

function p = empty_pinfo(Opt)
p = struct('solid_mode',Opt.init.solid,'fluid_mode',Opt.init.fluid, ...
           'Phi_s_time',0,'Phi_f_time',0,'prolongation_time',0, ...
           'coarse_setup_time',0,'random_generation_time',0, ...
           'total_time',0,'efficiency',struct());
end

function result = base_result(Opt, mesh_summary)
result = struct();
result.timestamp = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS'));
result.Opt = Opt;
result.mesh = mesh_summary;
end

function Z = fixed_random_block(n,m,seed,is_complex_problem)
rng(seed);
Z = randn(n,m);
if is_complex_problem, Z = Z + 1i*randn(n,m); end
end