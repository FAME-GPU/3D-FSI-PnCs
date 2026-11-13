function [v0, info] = init_eigs_v0_coupled_coarse(mesh_c, mesh_h, phys, wave_vector, T_s_h, T_f_h, opts)

if nargin < 7 || isempty(opts),opts=struct();end
if nargin < 5,T_s_h=[];end
if nargin < 6,T_f_h=[];end
opts=local_defaults_coupled(opts);
local_check_positive_integer(opts.nEig,'nEig');
if ~isvector(wave_vector)||numel(wave_vector)~=3
    error('wave_vector must be a 3-component Bloch vector.');
end
if ~isfield(phys,'meteria')||~isfield(phys,'lattice')
    error('phys must come from Physical_settings(model).');
end
all_clock=tic;
timing=struct();
t0=tic;
mesh_c=local_read_mesh(mesh_c);
mesh_h=local_read_mesh(mesh_h);
timing.mesh_load=toc(t0);

t0=tic;
mat_c=Assembly_mass_and_stiffness_matrix(mesh_c,phys.meteria);
[T_s_c,T_f_c]=Generate_Bloch_T_Matrix_fast( ...
    mesh_c,wave_vector,phys.lattice.lattice_vector);
K_s=local_hermitianize(T_s_c'*mat_c.K_s_tight*T_s_c);
K_f=local_hermitianize(T_f_c'*mat_c.K_f_tight*T_f_c);
M_s=local_hermitianize(T_s_c'*mat_c.M_s_tight*T_s_c);
M_f=local_hermitianize(T_f_c'*mat_c.M_f_tight*T_f_c);
C=T_f_c'*mat_c.C_tight*T_s_c;
clear mat_c;
timing.coarse_assembly=toc(t0);
ns=size(K_s,1);nf=size(K_f,1);n=ns+nf;
if n<=opts.nEig+1,error('nEig too large for coarse n=%d.',n);end
if opts.verbose
    fprintf('[CoupledCoarse] coarse DOFs: solid=%d, fluid=%d\n',ns,nf);
end

t0=tic;
P_s=Prec_setup_Ks_ict(K_s,opts.solidDroptol,opts.solidDiagcomp,true);
P_f=Prec_setup_Ks_ict(K_f,opts.fluidDroptol,opts.fluidDiagcomp,true);
timing.preconditioner=toc(t0);

stats=struct('op_calls',0,'s_calls',0,'f_calls',0, ...
    's_pcg_iter',0,'f_pcg_iter',0,'s_failed',0,'f_failed',0, ...
    's_pcg_time',0,'f_pcg_time',0);

ei=struct('isreal',false,'issym',false,'tol',opts.eigsTol, ...
    'maxit',opts.maxit,'disp',0);
ei.p=local_arnoldi_p(n,opts.nEig,opts.p);
t0=tic;
[V,D,flag_eigs]=eigs(@apply_AinvB,n,opts.nEig,'largestabs',ei);
timing.coarse_eigs=toc(t0);
[lambda_c,V] = local_positive_modes(1./diag(D),V,opts.realTol);
if numel(lambda_c)<opts.nEig
    warning('Only %d of the requested %d coarse eigenpairs are positive and essentially real.', ...
        numel(lambda_c),opts.nEig);
end
if isempty(lambda_c)
    error('No positive essentially-real coarse coupled eigenvalue found.');
end
if flag_eigs~=0
    warning('Coarse coupled eigs flag=%d; inspect original GEP residuals.',flag_eigs);
end
mode_idx=local_mode_indices(opts.modeIndices,numel(lambda_c));
selected_lambda=lambda_c(mode_idx);
phi_s_c=V(1:ns,mode_idx);
phi_f_c=V(ns+1:end,mode_idx);

t0=tic;
coarse_residual=zeros(numel(mode_idx),1);
for j=1:numel(mode_idx)
    us=phi_s_c(:,j);pf=phi_f_c(:,j);lam=selected_lambda(j);
    au=K_s*us-C'*pf;ap=K_f*pf;
    bu=M_s*us;bp=C*us+M_f*pf;
    coarse_residual(j)=norm([au-lam*bu;ap-lam*bp])/ ...
        (norm([au;ap])+abs(lam)*norm([bu;bp])+eps);
end
timing.coarse_residual=toc(t0);

t0=tic;
if isempty(T_s_h)||isempty(T_f_h)
    [T_s_new,T_f_new]=Generate_Bloch_T_Matrix_fast( ...
        mesh_h,wave_vector,phys.lattice.lattice_vector);
    if isempty(T_s_h),T_s_h=T_s_new;end
    if isempty(T_f_h),T_f_h=T_f_new;end
end
timing.fine_bloch=toc(t0);

popts=opts.prolongOpts;
popts.inputSpace='reduced';
popts.massOrth=false;
popts.computeSolid=true;
popts.computeFluid=true;
popts.verbose=opts.verbose && local_get(popts,'verbose',false);
t0=tic;
[U_h,P_h,prolong_info]=build_FSI_prolongation_init( ...
    mesh_c,mesh_h,T_s_c,T_f_c,T_s_h,T_f_h, ...
    phi_s_c,phi_f_c,[],[],popts);
timing.prolongation=toc(t0);

t0=tic;
Z_h=[U_h;P_h];
[v0,weights,column_norms]=local_mix_modes(Z_h,selected_lambda,opts);
timing.weighting=toc(t0);

info=struct();
info.method='coarse_full_coupled_eigs_prolongation';
info.options=opts;
info.coarse_dofs=[ns,nf];
info.fine_dofs=[size(T_s_h,2),size(T_f_h,2)];
info.coarse_eigs_flag=flag_eigs;
info.coarse_lambda_valid=lambda_c;
info.selected_mode_indices=mode_idx;
info.selected_lambda=selected_lambda;
info.selected_coarse_residual=coarse_residual;
info.weight_mode=opts.weightMode;
info.actual_weights=weights;
info.fine_mode_norms_before_normalization=column_norms;
info.inner=stats;
info.prolongation=prolong_info;
info.timings=timing;
info.total_time=toc(all_clock);
info.save_file=char(opts.saveFile);

if opts.verbose
    fprintf('[CoupledCoarse] v0=%d x 1, modes=%d, eigsFlag=%d, total=%.3f s, saved: %s\n', ...
        numel(v0),numel(mode_idx),flag_eigs,info.total_time,info.save_file);
end

    function y=apply_AinvB(x)
        
        y=zeros(size(x),'like',x);
        for kk=1:size(x,2)
            v=x(:,kk);
            u=v(1:ns);p=v(ns+1:end);
            rhs_p=C*u+M_f*p;
            tpcg=tic;
            [yp,ff,~,itf]=pcg(K_f,rhs_p,opts.innerTol,opts.innerMaxit,P_f.apply);
            stats.f_calls=stats.f_calls+1;
            stats.f_pcg_time=stats.f_pcg_time+toc(tpcg);
            stats.f_pcg_iter=stats.f_pcg_iter+local_iter(itf);
            stats.f_failed=stats.f_failed+(ff~=0);
            rhs_u=M_s*u+C'*yp;
            tpcg=tic;
            [yu,fs,~,its]=pcg(K_s,rhs_u,opts.innerTol,opts.innerMaxit,P_s.apply);
            stats.s_calls=stats.s_calls+1;
            stats.s_pcg_time=stats.s_pcg_time+toc(tpcg);
            stats.s_pcg_iter=stats.s_pcg_iter+local_iter(its);
            stats.s_failed=stats.s_failed+(fs~=0);
            y(:,kk)=[yu;yp];
            stats.op_calls=stats.op_calls+1;
        end
    end
end

function opts=local_defaults_sep(opts)
opts=local_defaults_common(opts);
opts=local_set(opts,'nSolid',8);
opts=local_set(opts,'nFluid',8);
opts=local_set(opts,'nCoupled',8);
opts=local_set(opts,'saveFile','eigs_init_sep_project.mat');
end

function opts=local_defaults_coupled(opts)
opts=local_defaults_common(opts);
opts=local_set(opts,'nEig',8);
opts=local_set(opts,'saveFile','eigs_init_coupled_coarse.mat');
end

function opts=local_defaults_common(opts)
opts=local_set(opts,'modeIndices',[]);
opts=local_set(opts,'weights',[]);
opts=local_set(opts,'weightMode','decay');
opts=local_set(opts,'eigsTol',1e-2);
opts=local_set(opts,'innerTol',1e-3);
opts=local_set(opts,'innerMaxit',200);
opts=local_set(opts,'maxit',300);
opts=local_set(opts,'p',[]);
opts=local_set(opts,'realTol',1e-6);
opts=local_set(opts,'solidDroptol',1e-3);
opts=local_set(opts,'solidDiagcomp',0);
opts=local_set(opts,'fluidDroptol',5e-3);
opts=local_set(opts,'fluidDiagcomp',1e-5);
opts=local_set(opts,'verbose',true);
opts=local_set(opts,'randomSeed',1);
opts=local_set(opts,'prolongOpts',struct());
end

function s=local_set(s,name,value)
if ~isfield(s,name)||isempty(s.(name)),s.(name)=value;end
end

function x=local_get(s,name,default)
if isfield(s,name)&&~isempty(s.(name)),x=s.(name);else,x=default;end
end

function mesh=local_read_mesh(mesh)
if isstring(mesh)||ischar(mesh)
    file=char(mesh);
    if ~isfile(file),error('Mesh file not found: %s',file);end
    s=load(file,'meshData');
    if ~isfield(s,'meshData')
        error('MAT file must contain variable meshData: %s',file);
    end
    mesh=s.meshData;
end
if isstruct(mesh)&&isfield(mesh,'meshData')&&~isfield(mesh,'nodes')
    mesh=mesh.meshData;
end
if ~isstruct(mesh)||~isfield(mesh,'nodes')||~isfield(mesh,'elements')
    error('mesh_c/mesh_h must be meshData or MAT path containing meshData.');
end
end

function H=local_hermitianize(H)
H=(H+H')*0.5;
end

function local_check_positive_integer(x,name)
if ~isscalar(x)||~isnumeric(x)||~isfinite(x)||x<1||x~=floor(x)
    error('%s must be a positive integer.',name);
end
end

function p=local_arnoldi_p(n,k,p_user)
if isempty(p_user)
    p=min(n,max(20,3*k));
else
    p=p_user;
end
if ~isscalar(p)||p<=k||p>n||p~=floor(p)
    error('Invalid eigs Arnoldi dimension p=%g: require k<p<=n.',p);
end
end

function [lambda,X]=local_positive_modes(lambda,X,rtol)
lambda=lambda(:);
good=isfinite(real(lambda))&isfinite(imag(lambda)) & real(lambda)>0 & ...
     abs(imag(lambda)) <= rtol .* max(1,abs(real(lambda)));
lambda=real(lambda(good));
X=X(:,good);
[lambda,ord]=sort(lambda,'ascend');
X=X(:,ord);
end

function Q=local_mass_orth(X,M)

G=local_hermitianize(X'*(M*X));
[R,p]=chol(full(G));
if p==0
    Q=X/R;
else
    [W,d]=eig(full(G),'vector');
    d=real(d(:));
    keep=d>max(d)*1e-11;
    if ~any(keep),error('Isolated M-Gram matrix is rank deficient.');end
    Q=(X*W(:,keep))*diag(1./sqrt(d(keep)));
    warning('M-Gram was singular/indefinite: retained %d directions.',nnz(keep));
end
end

function idx=local_mode_indices(idx,nvalid)
if isempty(idx),idx=1:nvalid;end
idx=idx(:).';
if isempty(idx)||any(~isfinite(idx))||any(idx~=round(idx))|| ...
        any(idx<1)||any(idx>nvalid)||numel(unique(idx))~=numel(idx)
    error('modeIndices must be distinct integers within 1:%d.',nvalid);
end
end

function [v0,weights,norms]=local_mix_modes(Z,lambda,opts)
if size(Z,2)~=numel(lambda)
    error('Number of columns must equal selected eigenvalue count.');
end
norms=sqrt(sum(abs(Z).^2,1)).';
if any(~isfinite(norms))||any(norms<=0)
    error('Prolongation returned zero/non-finite coupled mode.');
end

Z=bsxfun(@rdivide,Z,norms.');
k=size(Z,2);
if ~isempty(opts.weights)
    weights=opts.weights(:);
    if numel(weights)~=k||any(~isfinite(weights))
        error('opts.weights requires exactly %d finite coefficients.',k);
    end
else
    switch lower(char(opts.weightMode))
        case 'equal'
            weights=ones(k,1);
        case 'decay'
            weights=1./(1:k).';
        case 'first'
            weights=zeros(k,1);weights(1)=1;
        case 'inverse_lambda'
            weights=min(lambda)./lambda(:);
        case 'random'
            oldrng=rng;
            cleanup=onCleanup(@()rng(oldrng));
            rng(opts.randomSeed,'twister');
            weights=randn(k,1)+1i*randn(k,1);
        otherwise
            error('Unsupported weightMode: %s',char(opts.weightMode));
    end
end
if norm(weights)==0,error('Mode weights cannot all be zero.');end

weights=weights/norm(weights);
v0=Z*weights;
nv=norm(v0);
if ~isfinite(nv)||nv<1e-12
    error('Weighted modes nearly cancel; choose other weights.');
end
v0=full(v0/nv);
end

function k=local_iter(it)
if isempty(it),k=0;elseif isscalar(it),k=it;else,k=it(end);end
end

function local_save(v0,info,file)
file=char(file);
if isempty(file),error('opts.saveFile cannot be empty.');end
folder=fileparts(file);
if ~isempty(folder)&&~isfolder(folder),mkdir(folder);end
save(file,'v0','info','-v7.3');
end
