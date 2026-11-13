function cands = evaluate_candidates(red, K_s, K_f, M_s, M_f, C, V_s, V_f, maxCheck, cache)

if nargin < 10
    cache = [];
end

ncheck = min(maxCheck,numel(red.theta));
empty = struct('theta',[],'q_s',[],'q_f',[],'u_s',[],'u_f',[], ...
               'r_s',[],'r_f',[],'res',[],'den',[]);
cands = repmat(empty,ncheck,1);

use_cache = has_valid_cache(cache,V_s,V_f);

for j = 1:ncheck
    th = red.theta(j);
    q_s = red.qs(:,j);
    q_f = red.qf(:,j);

    nrm2 = real(q_s'*q_s + q_f'*q_f);
    if nrm2 > 0 && isfinite(nrm2)
        scale = sqrt(nrm2);
        q_s = q_s/scale;
        q_f = q_f/scale;
    end

    u_s = V_s*q_s;
    u_f = V_f*q_f;

    if use_cache
        Ks_us = cache.KsVs*q_s;
        Ms_us = cache.MsVs*q_s;
        Cu    = cache.CVs*q_s;

        Kf_uf = cache.KfVf*q_f;
        Mf_uf = cache.MfVf*q_f;
        Ctu   = cache.CtVf*q_f;
    else
        
        Ks_us = K_s*u_s;
        Ms_us = M_s*u_s;
        Cu    = C*u_s;

        Kf_uf = K_f*u_f;
        Mf_uf = M_f*u_f;
        Ctu   = C'*u_f;
    end

    r_s = Ks_us - th^2*Ms_us - th*Ctu;
    r_f = Kf_uf - th^2*Mf_uf - th*Cu;

    den = norm(Ks_us) + norm(Kf_uf) + abs(th)*(norm(Ctu)+norm(Cu)) + ...
          abs(th)^2*(norm(Ms_us)+norm(Mf_uf));
    if den == 0 || ~isfinite(den)
        den = 1;
    end
    res = sqrt(norm(r_s)^2+norm(r_f)^2)/den;

    cands(j).theta = th;
    cands(j).q_s = q_s;
    cands(j).q_f = q_f;
    cands(j).u_s = u_s;
    cands(j).u_f = u_f;
    cands(j).r_s = r_s;
    cands(j).r_f = r_f;
    cands(j).res = res;
    cands(j).den = den;
end
end

function tf = has_valid_cache(cache,V_s,V_f)
tf = ~isempty(cache) && isstruct(cache) && ...
     isfield(cache,'KsVs') && size(cache.KsVs,2) == size(V_s,2) && ...
     isfield(cache,'MsVs') && size(cache.MsVs,2) == size(V_s,2) && ...
     isfield(cache,'CVs')  && size(cache.CVs,2)  == size(V_s,2) && ...
     isfield(cache,'KfVf') && size(cache.KfVf,2) == size(V_f,2) && ...
     isfield(cache,'MfVf') && size(cache.MfVf,2) == size(V_f,2) && ...
     isfield(cache,'CtVf') && size(cache.CtVf,2) == size(V_f,2);
end
