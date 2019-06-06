function [X0, k] = congrad(A, b)
    n = length(b);
    X0 = zeros(n, 1);
    r0 = b - A * X0;
    P0 = r0;
    k = 0;
    while norm(r0) >= 10^(-12)
        alpha = r0' * r0 / (P0' * A * P0);
        X1 = X0 + alpha * P0;
        r1 = r0 - alpha * A * P0;
        beta = r1' * r1 / (r0' * r0);
        P1 = r1 + beta * P0;
        r0 = r1;
        X0 = X1;
        P0 = P1;
        k = k + 1;
    end
end
