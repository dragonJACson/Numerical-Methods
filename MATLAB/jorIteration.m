function x = jorIteration(A, b, x0, omega, eps, max_n)
    D = diag(diag(A));
    L = tril(A, -1);
    U = tril(A, 1);
    B = D \ ((1 - omega) * D - omega * (L + U));
    F = omega * D \ b;
    x = B * x0 + F;
    n = 1;
    while norm(x - x0, inf) >= eps
        x0 = x;
        x = B * x0 + F;
        n = n + 1;
        if n >= max_n
            return;
        end
    end
end
