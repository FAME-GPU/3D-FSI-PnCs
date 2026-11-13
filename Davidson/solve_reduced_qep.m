function red = solve_reduced_qep(K_s, K_f, M_s, M_f, C, V_s, V_f, target_theta, opts, cache)

if nargin < 10
    cache = [];
end

ms = size(V_s,2);
mf = size(V_f,2);
m  = ms+mf;

red = struct();
red.theta = [];
red.qs = [];
red.qf = [];
red.gamma_elimination_used = false;
red.C00_relative = NaN;
red.path = '';
red.used_basis_cache = false;

if ~isempty(cache) && isstruct(cache) && ...
        isfield(cache,'KsVs') && size(cache.KsVs,2) == ms && ...
        isfield(cache,'KfVf') && size(cache.KfVf,2) == mf && ...
        isfield(cache,'CVs')  && size(cache.CVs,2)  == ms
    KsVs = cache.KsVs;
    KfVf = cache.KfVf;
    CVs  = cache.CVs;
    red.used_basis_cache = true;
else
    KsVs = K_s*V_s;
    KfVf = K_f*V_f;
    CVs  = C*V_s;
end

Ksh = V_s'*KsVs;
Kfh = V_f'*KfVf;
Ksh = 0.5*(Ksh+Ksh');
Kfh = 0.5*(Kfh+Kfh');
Chat = V_f'*CVs;

use_gamma = isfield(opts,'gamma') && isfield(opts.gamma,'enabled') && opts.gamma.enabled;
ns0 = 0;
nf0 = 0;
if use_gamma
    ns0 = get_field(opts.gamma,'nullDimS',0);
    nf0 = get_field(opts.gamma,'nullDimF',0);
end
eliminate_zero = use_gamma && (ns0+nf0 > 0) && get_field(opts.gamma,'eliminateZero',true);

if eliminate_zero
    couplingTol = get_field(opts.gamma,'couplingTol',1e-8);

    if ms <= ns0 || mf <= nf0
        red.path = 'gamma_no_working_space';
        return;
    end

    C00 = Chat(1:nf0,1:ns0);
    C0s = Chat(1:nf0,ns0+1:end);
    Cf0 = Chat(nf0+1:end,1:ns0);
    Cfs = Chat(nf0+1:end,ns0+1:end);

    if ns0 > 0 && nf0 > 0
        red.C00_relative = norm(C00,'fro')/max(1,norm(Chat,'fro'));
        coupling_ok = red.C00_relative <= couplingTol;
    else
        
        red.C00_relative = 0;
        coupling_ok = true;
    end

    if coupling_ok
        mws = ms-ns0;
        mwf = mf-nf0;

        Ksw = Ksh(ns0+1:end,ns0+1:end);
        Kfw = Kfh(nf0+1:end,nf0+1:end);

        Ksw_tilde = Ksw + C0s'*C0s;
        Kfw_tilde = Kfw + Cf0*Cf0';
        Ksw_tilde = 0.5*(Ksw_tilde+Ksw_tilde');
        Kfw_tilde = 0.5*(Kfw_tilde+Kfw_tilde');

        mw = mws+mwf;
        A0 = blkdiag(Ksw_tilde,Kfw_tilde);
        A1 = -[zeros(mws,mws), Cfs'; Cfs, zeros(mwf,mwf)];
        A2 = -eye(mw);

        [Xw,E] = polyeig(full(A0),full(A1),full(A2));
        if ns0 > 0 && nf0 > 0
            red.path = 'gamma_eliminated_both';
        elseif ns0 > 0
            red.path = 'gamma_eliminated_solid';
        else
            red.path = 'gamma_eliminated_fluid';
        end

        theta_all = E(:);
        good = admissible_roots(theta_all,opts);
        theta = real(theta_all(good));
        Xw = Xw(:,good);

        if isempty(theta)
            return;
        end

        q = zeros(m,numel(theta),'like',Xw);
        for j = 1:numel(theta)
            th = theta(j);
            bs = Xw(1:mws,j);
            bf = Xw(mws+1:end,j);
            as = -(Cf0'*bf)/th;
            af = -(C0s*bs)/th;

            qs = [as;bs];
            qf = [af;bf];
            qj = [qs;qf];
            nq = norm(qj);
            if nq > 0
                qj = qj/nq;
            end
            q(:,j) = qj;
        end

        [theta,q] = select_and_sort(theta,q,target_theta,opts);
        red.theta = theta(:);
        red.qs = q(1:ms,:);
        red.qf = q(ms+1:end,:);
        red.gamma_elimination_used = true;
        return;
    else
        warn_c00_once(red.C00_relative,couplingTol);
    end
end

A0 = blkdiag(Ksh,Kfh);
A1 = -[zeros(ms,ms), Chat'; Chat, zeros(mf,mf)];
A2 = -eye(m);

[X,E] = polyeig(full(A0),full(A1),full(A2));
if use_gamma
    if ns0+nf0 > 0
        red.path = 'gamma_full_fallback';
    else
        red.path = 'gamma_no_locked_modes';
    end
else
    red.path = 'standard_full';
end

theta_all = E(:);
good = admissible_roots(theta_all,opts);
theta = real(theta_all(good));
q = X(:,good);

if isempty(theta)
    return;
end

for j = 1:size(q,2)
    nq = norm(q(:,j));
    if nq > 0
        q(:,j) = q(:,j)/nq;
    end
end

[theta,q] = select_and_sort(theta,q,target_theta,opts);
red.theta = theta(:);
red.qs = q(1:ms,:);
red.qf = q(ms+1:end,:);
end

function good = admissible_roots(theta_all,opts)
minOmega = get_field(opts,'minOmega',0);
realTol = get_field(opts,'realTol',1e-9);
good = isfinite(theta_all) & abs(theta_all) > minOmega & ...
       abs(imag(theta_all)) <= realTol.*max(1,abs(theta_all)) & ...
       real(theta_all) > minOmega;
end

function [theta,q] = select_and_sort(theta,q,target_theta,opts)
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
theta = theta(ord);
q = q(:,ord);
end

function warn_c00_once(value,tol)
persistent already_warned
if isempty(already_warned) || ~already_warned
    warning(['Gamma exact two-sided zero-space elimination disabled because ', ...
             '||Z_f'' C Z_s|| is not negligible in the reduced basis ', ...
             '(relative %.3e > %.3e). Using the full reduced QEP with zero-root filtering.'], ...
            value,tol);
    already_warned = true;
end
end

function value = get_field(s,name,default_value)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default_value;
end
end
