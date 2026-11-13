function Q = m_orthonormalize(X, M, tol, maxReorth)

    if nargin < 3 || isempty(tol)
        tol = 1e-12;
    end
    if nargin < 4 || isempty(maxReorth)
        maxReorth = 5;
    end

    if isempty(X)
        Q = X;
        return;
    end

    Q = zeros(size(X,1), 0, 'like', X);

    for j = 1:size(X,2)
        v = X(:,j);

        for rep = 1:maxReorth
            if ~isempty(Q)
                h = Q' * (M * v);
                v = v - Q * h;
            end
        end

        Mv = M * v;
        nv = sqrt(abs(v' * Mv));
        if isfinite(nv) && nv > tol
            Q(:,end+1) = v / nv; 
        end
    end
end