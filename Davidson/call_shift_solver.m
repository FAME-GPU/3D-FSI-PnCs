function [z, solve_info] = call_shift_solver(fh, theta, b, A, K, M, opts)
    nf = nargin(fh);
    if nf == 1
        z = fh(b);
    elseif nf == 2
        [z, solve_info] = fh(theta, b);
    elseif nf == 3
        z = fh(theta, b, A);
    elseif nf == 4
        z = fh(theta, b, A, K);
    elseif nf == 5
        z = fh(theta, b, A, K, M);
    else
        z = fh(theta, b, A, K, M, opts);
    end
end