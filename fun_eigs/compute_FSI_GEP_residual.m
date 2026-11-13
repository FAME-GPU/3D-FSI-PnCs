function [residuals, Rmax, detail] = compute_FSI_GEP_residual(K_s, K_f, M_s, M_f, C, lambda, U_s, P)

lambda = lambda(:);
nv = numel(lambda);

if size(U_s,2) ~= nv || size(P,2) ~= nv
    error('Number of columns in U_s and P must match numel(lambda).');
end

residuals = NaN(nv,1);
rs_norm = NaN(nv,1);
rf_norm = NaN(nv,1);
denominator = NaN(nv,1);

for j = 1:nv
    lam = lambda(j);
    us = U_s(:,j);
    p  = P(:,j);

    Ks_us = K_s * us;
    Ms_us = M_s * us;
    Kf_p  = K_f * p;
    Mf_p  = M_f * p;
    Cp    = C' * p;
    Cus   = C * us;

    r_s = Ks_us - Cp - lam * Ms_us;
    r_f = Kf_p - lam * (Cus + Mf_p);

    den = norm(Ks_us) + norm(Cp) + abs(lam)*norm(Ms_us) + ...
          norm(Kf_p) + abs(lam)*(norm(Cus) + norm(Mf_p));
    if ~isfinite(den) || den <= 0
        den = 1;
    end

    rs_norm(j) = norm(r_s);
    rf_norm(j) = norm(r_f);
    denominator(j) = den;
    residuals(j) = hypot(rs_norm(j), rf_norm(j)) / den;
end

if isempty(residuals)
    Rmax = NaN;
else
    finite_res = residuals(isfinite(residuals));
    if isempty(finite_res)
        Rmax = NaN;
    else
        Rmax = max(finite_res);
    end
end

detail = struct();
detail.solid_residual_norm = rs_norm;
detail.fluid_residual_norm = rf_norm;
detail.denominator = denominator;
end
