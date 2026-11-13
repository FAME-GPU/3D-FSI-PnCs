function [ev, ew, iter] = Lanczos(fncA, fncB, fncBinv, n, k, tol, maxiter, p, v0, flag_gpu)

if flag_gpu
    v0 = gpuArray(v0);
end
v0 = v0 / sqrt(abs(v0' * fncB(v0)));

[ew, ev, iter, res] = Lanczos_iteration(fncA, fncB, fncBinv, n, k, tol, maxiter, p, v0, flag_gpu);

end

function [d, U, kk, res] = Lanczos_iteration(fncA, fncB, fncBinv, n, k, tol, maxiter, p, v, flag_gpu)

eps = 2.2204e-16;

Q = zeros(n, p);
if flag_gpu
    Q = gpuArray(Q);
end

k0 = k;  
Alpha = [];
Beta = [];
alpha = 0;
beta = 0;
d = [];
c = [];
justRestarted = false;
sizeV = 1;
kk = 0;

for mm = 1 : maxiter
    
    for jj = sizeV : p
        
        Q(:, jj) = v;
        r = fncA(Q(:, jj));
        alpha = real(r' * Q(:, jj));
        
        if jj == 1
            pr = fncBinv(r);
            r = pr - alpha * Q(:, jj);
        elseif justRestarted
            pr = fncBinv(r);
            r = GramSchmitGEP(Q, pr, fncB, jj);
            justRestarted = false;
        else
            pr = fncBinv(r);
            r = pr - alpha * Q(:, jj) - beta * Q(:, jj - 1);
            r = r - fncB(Q(:, jj - 1))' * r * Q(:, jj - 1);
        end
        
        r = r - fncB(Q(:, jj))' * r * Q(:, jj);
        
        r = GramSchmitGEP(Q, r, fncB, jj);
        
        rr = fncB(r);
        beta = sqrt(r' * rr);
        v = r / beta;
        
        if flag_gpu
            alphac = gather(alpha);
            betac = gather(beta);
        else
            alphac = alpha;
            betac = beta;
        end
                
        Alpha = [Alpha; alphac]; 
        Beta = [Beta; betac];
        
        kk = kk + 1;
        
    end
    
    H1 = diag(Alpha) + diag(Beta(1:end-1), 1) + diag(Beta(1:end-1), -1);
    H = blkdiag(diag(d), H1);
    if ~isempty(d)
        H( 1:k , k+1 ) = c';
        H( k+1 , 1: k ) = c;
    end
    Alpha = [];
    Beta = [];
    
    [U, d] = eig(H, 'vector');
    
    res = abs(beta * U(end, :));
    
    [~, ind] = sort(real(d), 'descend');  
    d = d(ind);
    res = res(ind);
    
    isNotConverged = ~(res(1:k0)' < tol*max(eps^(2/3), abs(d(1:k0))));
    nconv = nnz(~isNotConverged);
    
    if nconv >= k0 || mm == maxiter
        
        break;
    else
        
        k = k0 + min(nconv, floor((p - k0) / 2));
        if k == 1 && p > 3
            k = floor(p / 2);
        end
    end
    
    ind = ind(1 : k);
    U = U(:, ind);
    
    if flag_gpu
        Q(:, 1 : k) = Q * gpuArray(U);
    else
        Q(:, 1 : k) = Q * U;
    end
    d = d(1 : k);
    c = betac * U(end, :);
    justRestarted = true;
    sizeV = k + 1;
    
end

U = U(:, ind(1 : k0));
d = d(1 : k0);

if flag_gpu
    U = Q * gpuArray(U);
else
    U = Q * U;
end
    
end

function r = GramSchmitGEP(V, r, fncB, jj)

for k = 1 : jj
    tmp = fncB(V(:, k))' * r;
    r = r - tmp * V(:, k);

end

end
