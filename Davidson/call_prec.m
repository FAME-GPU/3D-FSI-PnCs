function z = call_prec(fh, theta, b)
    nf = nargin(fh);
    if nf == 1
        z = fh(b);
    else
        z = fh(theta, b);
    end
end