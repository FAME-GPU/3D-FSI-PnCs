function A = hermitianize(A)
    A = 0.5 * (A + A');
end