function Opt = Initialize_parameters(model)

Opt.flag_gpu = 0;
Opt.model = model;  

Opt.nev = 10;

if strcmp(Opt.model,'AlWater')
    Opt.k_index = 5;
elseif strcmp(Opt.model,'SC')
    Opt.k_index = 10;
elseif strcmp(Opt.model,'SC_Quad') || strcmp(Opt.model,'SC_InnerFluid')
    Opt.k_index = 40;
end  

Opt.target_omega = 0;
Opt.eig_tol = 1e-10;
Opt.final_residual_tol = Opt.eig_tol;
Opt.max_outer_iterations = 200;  

Opt.subspace.init_dim = 12;
Opt.subspace.max_dim = 50;
Opt.subspace.restart_dim = 12;
Opt.subspace.expand_pairs = 2;
Opt.subspace.extra_candidates = 2;
Opt.add_cross = false;

Opt.linear.solver = 'minres';
Opt.linear.tol = 1e-2;     
Opt.linear.maxit = 10;     
Opt.linear.allow_fallback = false;

Opt.preconditioner.enable = true;
if strcmp(Opt.model,'SC') || strcmp(Opt.model,'AlWater')
    Opt.preconditioner.solid.droptol = 1e-3; 
    Opt.preconditioner.solid.lambda = 1e-4; 
    Opt.preconditioner.fluid.droptol = 1e-4;
    Opt.preconditioner.fluid.lambda = 1e-4; 
elseif strcmp(Opt.model,'SC_Quad')
    Opt.preconditioner.solid.droptol = 1e-3; 
    Opt.preconditioner.solid.lambda = 1e-5;  
    Opt.preconditioner.fluid.droptol = 1e-3;
    Opt.preconditioner.fluid.lambda = 0; 
end

Opt.preconditioner.solid.use_reorder = true;
Opt.preconditioner.fluid.use_reorder = true;

Opt.gamma.gammaAbsTol = 1e-12;
Opt.gamma.eliminateZero = true;       
Opt.gamma.useLift = true;             
Opt.gamma.nullTol = 1e-8;
Opt.gamma.couplingTol = 1e-8;
Opt.gamma.keepNullVectors = false;     
Opt.gamma.liftTau_s = Opt.preconditioner.solid.lambda;
Opt.gamma.liftTau_f = Opt.preconditioner.fluid.lambda;

Opt.init.solid = 'coarse';    
Opt.init.fluid = 'coarse';    
Opt.init.verbose = true;

Opt.Lanczos.eigensolver.tol = 1e-2;
Opt.Lanczos.eigensolver.itmax = 1e2;
Opt.Lanczos.LSsolver.tol = 1e-2;
Opt.Lanczos.LSsolver.itmax = 1e3;
Opt.Lanczos.subspace_dim = 20;
if strcmp(Opt.model,'SC') || strcmp(Opt.model,'AlWater')
    Opt.Lanczos.solid.droptol = 1e-3; 
    Opt.Lanczos.solid.diagcomp = 0;  
    Opt.Lanczos.fluid.droptol = 5e-3;
    Opt.Lanczos.fluid.diagcomp = 1e-5; 
elseif strcmp(Opt.model,'SC_Quad')
    Opt.Lanczos.solid.droptol = 1e-4; 
    Opt.Lanczos.solid.diagcomp = 0;  
    Opt.Lanczos.fluid.droptol = 1e-4;
    Opt.Lanczos.fluid.diagcomp = 0; 
end

Opt.stats.enable = true;         
Opt.stats.save_history = true;   
Opt.stats.track_matvec = true;   
Opt.stats.keep_inner_log = false;
Opt.stats.history_count = 10;
Opt.stats.detailed_timing = false; 
Opt.cacheBasisActions = true;

Opt.random_seed = 1;
Opt.reference_lambda = [];

Opt.eigs.tol = 1e-8;
Opt.eigs.inner_tol = 1e-8;
Opt.eigs.inner_maxit = 1000;

Opt.eigs.disp = 1;
Opt.eigs.profile_efficiency = true;

end
