function Prec = Prec_setup_fluid_ict(Kf, Mf, droptol, lambda, use_reorder)

    if nargin < 4
        use_reorder = true;
    end

    KMs = Kf + lambda * Mf;
    KMs = (KMs + KMs') / 2;
    n  = size(KMs, 1);
        
    d0 = real(diag(KMs));
    if min(d0) <= 0
        error('Ks has nonpositive diagonal entries: min(diag(Ks)) = %.4e', min(d0));
    end

    if use_reorder
        p = symamd(KMs);
    else
        p = (1:n).';
    end

    Kp = KMs(p, p);
    Kp = (Kp + Kp') / 2;

    d = real(diag(Kp));
    s = 1 ./ sqrt(d);
    S = spdiags(s, 0, n, n);

    Kh = S * Kp * S;
    Kh = (Kh + Kh') / 2;

    setup.type     = 'ict';
    setup.droptol  = droptol;
    
    setup.michol   = 'off';

    L = ichol(Kh, setup);

    Prec.L = L;
    Prec.p = p;
    Prec.s = s;
    Prec.apply = @(r) apply_Ks_prec_safe(r, L, p, s);
end

function y = apply_Ks_prec_safe(r, L, p, s)

    y = zeros(size(r), 'like', r);

    rp = r(p);
    rhs = s .* rp;

    z = L \ rhs;
    z = L' \ z;

    yp = s .* z;
    y(p) = yp;

end