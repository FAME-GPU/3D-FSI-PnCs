function [Z_s, Z_f, ginfo] = build_gamma_nullspaces(K_s, K_f, M_s, M_f, C, opts)

g = opts.gamma;
lockSolid = g.lockSolid;
lockFluid = g.lockFluid;
nullTol = g.nullTol;
couplingTol = g.couplingTol;

ns = size(K_s,1);
nf = size(K_f,1);

Z_s = zeros(ns,0);
Z_f = zeros(nf,0);
if ~isreal(K_s) || ~isreal(M_s), Z_s = complex(Z_s); end
if ~isreal(K_f) || ~isreal(M_f), Z_f = complex(Z_f); end

solid_null_res = NaN;
fluid_null_res = NaN;

if lockSolid
    if isfield(g,'Zs') && ~isempty(g.Zs)
        Z_s_raw = g.Zs;
    else
        if mod(ns,3) ~= 0
            error(['Cannot construct the three rigid translations automatically because ', ...
                   'size(K_s,1) is not divisible by 3. Supply opts.gamma.Zs explicitly.']);
        end

        Z_s_raw = zeros(ns,3);
        if ~isreal(K_s), Z_s_raw = complex(Z_s_raw); end
        Z_s_raw(1:3:end,1) = 1;
        Z_s_raw(2:3:end,2) = 1;
        Z_s_raw(3:3:end,3) = 1;
    end

    Z_s = mass_orthonormalize(Z_s_raw,M_s,'solid Gamma nullspace');
    solid_null_res = relative_null_residual(K_s,Z_s);

    if solid_null_res > nullTol
        error(['lockSolid=true, but the constructed solid translations are not in ker(K_s): ', ...
               'relative residual %.3e > nullTol %.3e. ', ...
               'If K_s is nonsingular because of boundary conditions, set ', ...
               'opts.gamma.lockSolid=false. Otherwise select the correct DOF layout ', ...
               'or supply opts.gamma.Zs explicitly.'],solid_null_res,nullTol);
    end
end

if lockFluid
    if isfield(g,'Zf') && ~isempty(g.Zf)
        Z_f_raw = g.Zf;
    else
        Z_f_raw = ones(nf,1);
        if ~isreal(K_f), Z_f_raw = complex(Z_f_raw); end
    end

    Z_f = mass_orthonormalize(Z_f_raw,M_f,'fluid Gamma nullspace');
    fluid_null_res = relative_null_residual(K_f,Z_f);

    if fluid_null_res > nullTol
        error(['lockFluid=true, but the constant fluid vector is not in ker(K_f): ', ...
               'relative residual %.3e > nullTol %.3e. ', ...
               'If K_f is nonsingular because of boundary conditions, set ', ...
               'opts.gamma.lockFluid=false. Otherwise supply opts.gamma.Zf if needed.'], ...
              fluid_null_res,nullTol);
    end
end

if ~isempty(Z_s) && ~isempty(Z_f)
    C00 = Z_f'*(C*Z_s);
    coupling_scale = max(1,norm(C,'fro'));
    c00_rel = norm(C00,'fro')/coupling_scale;
    if c00_rel > couplingTol
        warning(['Gamma nullspace coupling Z_f''*C*Z_s is not negligible: %.3e > %.3e. ', ...
                 'The reduced solver will fall back to the full reduced QEP instead ', ...
                 'of the exact two-sided zero-space elimination.'],c00_rel,couplingTol);
    end
else
    C00 = zeros(size(Z_f,2),size(Z_s,2));
    c00_rel = 0;
end

ginfo = struct();
ginfo.enabled = true;
ginfo.lock_solid = lockSolid;
ginfo.lock_fluid = lockFluid;
ginfo.null_dim_s = size(Z_s,2);
ginfo.null_dim_f = size(Z_f,2);
ginfo.zero_mode_count_qep = ginfo.null_dim_s+ginfo.null_dim_f;
ginfo.solid_null_residual = solid_null_res;
ginfo.fluid_null_residual = fluid_null_res;
ginfo.C00 = C00;
ginfo.C00_relative = c00_rel;
ginfo.nullTol = nullTol;
ginfo.couplingTol = couplingTol;
end

function Z = mass_orthonormalize(Z,M,label)
G = Z'*(M*Z);
G = 0.5*(G+G');
try
    R = chol(full(G));
    Z = Z/R;
catch
    [U,D] = eig(full(G));
    d = real(diag(D));
    dmax = max(abs(d));
    if isempty(d) || dmax == 0 || any(d <= max(1,dmax)*1e-13)
        error('%s is rank deficient in the M-inner product.',label);
    end
    Z = Z*U*diag(1./sqrt(d));
end
end

function r = relative_null_residual(K,Z)
if isempty(Z)
    r = NaN;
    return;
end
num = norm(K*Z,'fro');
den = norm(K,'fro')*max(norm(Z,'fro'),1);
if den == 0
    r = num;
else
    r = num/den;
end
end