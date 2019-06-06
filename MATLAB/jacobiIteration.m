function [x, n] = jacobiIteration(A, b, x0, eps, max_n)
    D = diag(diag(A));
    L = tril(A, -1);
    U = triu(A, 1);
    B = -D \ (L + U);
    F = D \ b;
    x = B * x0 + F;
    n = 1;
    while norm(x - x0) >= eps
        x0 = x;
        x = B * x0 + F;
        n = n + 1;
        if n >= max_n
            fprintf('迭代次数过多');
            return;
        end
    end
    fprintf('迭代了 %f 次。\n', n);
end

